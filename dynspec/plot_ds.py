import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys
import argparse

from timing import *


def create_axs(figsize=(6, 10)):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(3, 1, figure=fig)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1:, 0], sharex=ax0)
    return ax0, ax1

def main():
    parser = argparse.ArgumentParser(description="Plot dynamic spectrum")

    parser.add_argument('ds_file', help="The name of the input .pkl file")
    parser.add_argument('--output_image', help="The filename to use for the average DS image. If not given, will plt.show().")
    parser.add_argument('--fscrunch_factor', type=int, help="Average this many channels together for the ouput DS image.")
    parser.add_argument('--tscrunch_factor', type=int, help="Average this many time bins together for the ouput DS image.")
    parser.add_argument('--add_predicted_toas', action='store_true', help="Add lines on plot to show when pulses are expected to arrive.")
    parser.add_argument('--dm', type=float, help="Dedisperse the spectrum using this DM (default is use the DM in timing.py)")
    parser.add_argument('--use_bin_units', action='store_true', help="Use bin/chan numbers for time/freq axes instead of physical units")
    parser.add_argument('--xlim', nargs=2, type=float, help="Zoom into this x-axis range")
    parser.add_argument('--pol', choices=['I', 'Q', 'U', 'V'], default='I', help="Which polarisation to plot (default = 'I')")
    parser.add_argument('--do_not_apply_pb_corr', action='store_true', help="Explicitly do not apply the primary beam correction, even if such a correction is available")
    parser.add_argument('--time_format', default="mjd", help="Display time in x-axis label using this (astropy) format (default='mjd')")
    #parser.add_argument('--rm', type=float, help="Apply de-Faraday rotation using this rotation measure")

    args = parser.parse_args()

    ephemeris = get_J1755_ephemeris()
    src_dm = ephemeris['DM']

    # Setup output plot specs
    axs = create_axs()

    #print(f"Opening {args.ds_file}...")
    dat = np.load(args.ds_file, allow_pickle=True)
    t = Time(dat['TIMES']/86400.0, scale='utc', format='mjd')
    dt = t[1] - t[0]
    f = dat['FREQS'] * u.Hz
    df = f[1] - f[0]
    f += df/2 # <--- a hack because the labels should be in the middle of the channels instead of the lower edge

    pb_corr = not args.do_not_apply_pb_corr
    S, fmask, _ = Stokes_ds(dat, pol=args.pol, pb_corr=pb_corr)

    if args.dm == 0.0:
        Sdd = S
        src_dm = 0.0 * u.pc / u.cm**3
    else:
        # Dedispersion won't work with nans
        S[np.isnan(S)] = 0.0

        if args.dm is not None:
            src_dm = args.dm * u.pc / u.cm**3

        Sdd = dedisperse_ds(S, src_dm, f, np.mean(f), dt)

    if 'DM' in dat.keys():
        total_dm = dat['DM']*u.pc/u.cm**3 + src_dm
    else:
        total_dm = src_dm

    # Reapply channel flags
    if len(fmask) > 0:
        Sdd[:,fmask] = np.nan

    # Form the profile
    lightcurve = np.nanmean(Sdd, axis=1)

    # Do a bit of extra averaging in frequency
    if args.fscrunch_factor is not None:
        Sdd_scr = np.nanmean(np.reshape(Sdd, (Sdd.shape[0], -1, args.fscrunch_factor)), axis=-1)
        fscr = np.nanmean(np.reshape(f, (-1, args.fscrunch_factor)), axis=-1)
    else:
        Sdd_scr = Sdd
        fscr = f

    # If add_predicted_toas is requested, add them!
    predicted_topo_toas = []
    t_base = (t - t[0]).to('s')
    if args.add_predicted_toas:
        src_period = ephemeris['period']
        src_pepoch = ephemeris['PEPOCH']
        src_coord = ephemeris['coord']

        # Subtract a baseline from the lightcurve to estimate noise
        A = np.vstack([t_base.value, np.ones_like(t_base.value)]).T
        m, b = np.linalg.lstsq(A, lightcurve, rcond=None)[0]
        debased_lightcurve = lightcurve - (m*t_base.value + b)
        noise = np.std(debased_lightcurve)

        # Take the beginning and end of the observation and convert them to barycentric
        t_range = Time([dat['TIMES'][0]/86400, dat['TIMES'][-1]/86400],
                       scale='utc', format='mjd',
                       location=EarthLocation.of_site(dat['TELESCOPE']))

        bary_range = t_range + t_range.light_travel_time(src_coord, ephemeris='jpl') - calc_dmdelay(src_dm, np.mean(f), np.inf*u.MHz)

        # Convert them to integer pulse numbers that can be used for indexing
        pulse_range = np.ceil(((bary_range - src_pepoch) / src_period).decompose()).astype(int)
        for pulse in range(*pulse_range):
            predicted_bary_toa = pulse*src_period + src_pepoch
            print(predicted_bary_toa.mjd, 3*noise)
            predicted_topo_toa = predicted_bary_toa - predicted_bary_toa.light_travel_time(src_coord, ephemeris='jpl', location=EarthLocation.of_site(dat['TELESCOPE'])) + calc_dmdelay(src_dm, np.mean(f), np.inf*u.MHz)
            predicted_topo_toas.append(predicted_topo_toa)

    # Time scrunching, if requested
    if args.tscrunch_factor is not None and args.tscrunch_factor > 1:
        nbins, nchan = Sdd_scr.shape
        try:
            Sdd_scr = np.nanmean(np.reshape(Sdd_scr, (-1, args.tscrunch_factor, nchan)), axis=1)
            t_base = np.nanmean(np.reshape(t_base, (-1, args.tscrunch_factor)), axis=1)
            lightcurve = np.nanmean(np.reshape(lightcurve, (-1, args.tscrunch_factor)), axis=1)
        except:
            print(f"warning: Couldn't scrunch bins. Possible reason: {args.tscrunch_factor} doesn't divide {nbins} evenly")

    # Plots
    if args.use_bin_units:
        axs[0].plot(lightcurve*1e3) # convert to mJy
        axs[1].pcolormesh(Sdd_scr.T)
        for predicted_topo_toa in predicted_topo_toas:
            axs[0].axvline((predicted_topo_toa.mjd - t[0].mjd)/dt.to('d').value, c='r', ls='--', alpha=0.3)

        axs[-1].set_xlabel("Time (bins)")
        axs[1].set_ylabel("Frequency (channels)")
    else:
        axs[0].plot(t_base, lightcurve*1e3)
        axs[1].pcolormesh(t_base.to('s').value, fscr.to('MHz').value, Sdd_scr.T)
        for predicted_topo_toa in predicted_topo_toas:
            axs[0].axvline((predicted_topo_toa - t[0]).to('s').value, c='r', ls='--', alpha=0.3)

        axs[-1].set_xlabel(f"Time (s) after {args.time_format.upper()} {t[0].to_value(args.time_format)}")
        axs[1].set_ylabel("Frequency (MHz)")
    axs[0].set_ylabel("Flux density (mJy/beam)")

    axs[0].set_title(f"Stokes {args.pol}, dedispersed to {total_dm}")

    if args.xlim is not None:
        axs[0].set_xlim(args.xlim)

    if args.output_image is not None:
        plt.tight_layout()
        plt.savefig(args.output_image)
    else:
        plt.show()

if __name__ == '__main__':
    main()
