import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys
import argparse

from timing import *
from forward_model import exponnorm

def cost_function(params, DM, τ_1GHz, t, f, ds, counts):
    width_s, toa_offset_s, brightness_in_lowest_channel, spectral_idx = params
    width = width_s * u.s
    toa_offset = toa_offset_s * u.s

    # Scattering timescale
    T, F = np.meshgrid(t, f)
    τ = τ_1GHz * (f/(1*u.GHz)).decompose()**(-4)
    TAU = τ[:,np.newaxis]

    # Apply DM
    T -= calc_dmdelay(DM, f, np.inf*u.MHz)[:,np.newaxis]

    # Form the pulses
    pulses = exponnorm(T, tao_offset, width, TAU)

    # Apply the spectral index
    pulses *= brightness_in_lowest_channel * (f[:,np.newaxis]/f[0])**spectral_idx

    # Get the weighted residuals squared
    residuals = (pulses - ds)**2 * counts

    return np.sum(residuals)


def main():
    parser = argparse.ArgumentParser(description="Stack dynamic spectra together according to a folding ephemeris")

    parser.add_argument('--dedisperse', action='store_true', help="Dedisperse (using DM in timing.py)")

    parser.add_argument('bin_size', type=float, help="The size of a time bin, in seconds")
    parser.add_argument('--fscrunch_factor', type=int, help="How many frequency channels to average together in the final output plot")
    parser.add_argument('--draw_dm_sweep', type=float, nargs='*', help="Time(s) (at centre frequency) at which to draw DM sweep curves in the synamic spectra panel")

    parser.add_argument('output_pulsestack_image', help="The filename to use for the pulsestack image")
    parser.add_argument('output_average_ds_image', help="The filename to use for the average DS image")

    parser.add_argument('ds_files', nargs='*', help="The names of the input .pkl files, without the .pkl extension")

    args = parser.parse_args()

    ephemeris = get_J1755_ephemeris()
    src_coord = ephemeris['coord']
    src_period = ephemeris['period']
    src_pepoch = ephemeris['PEPOCH']
    if args.dedisperse:
        src_dm = ephemeris['DM']

    # Binning parameters
    bin_size = args.bin_size * u.s
    bin_size_phase = (bin_size / src_period).decompose()

    fscrunch_factor = args.fscrunch_factor

    # Output filenames
    pulsestack_png = args.output_pulsestack_image
    average_ds_png = args.output_average_ds_image

    min_phase_bin = np.inf
    max_phase_bin = -np.inf

    pkls = args.ds_files
    dats = []
    phase_bins_list = []

    # Load data, and convert the times to phases and phase bins
    for pkl in pkls:
        print(f"Opening {pkl}...")
        dat = np.load(pkl, allow_pickle=True)
        dat['OBSID'] = pkl[:-4].replace('_', ' ')

        location = EarthLocation.of_site(dat['TELESCOPE'])
        times = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=location)
        barytimes = barycentre(times, src_coord)

        _, phases = fold(barytimes, src_period, src_pepoch)
        phase_bins = np.round(phases / bin_size_phase).astype(int)
        phase_bins_list.append(phase_bins)

        if min(phase_bins) < min_phase_bin:
            min_phase_bin = int(np.round(min(phase_bins)))
        if max(phase_bins) > max_phase_bin:
            max_phase_bin = int(np.round(max(phase_bins)))

        dat['PHASE_BINS'] = phase_bins
        dats.append(dat)

    f = dats[0]['FREQS'] * u.Hz
    
    # Work out the dimensions of the final output array
    nphase_bins = max_phase_bin - min_phase_bin + 1
    nchans = len(dats[0]['FREQS'])

    output_ds = np.zeros((nphase_bins, nchans))
    counts = np.zeros((nphase_bins, nchans)) # For tracking how many items contribute to the same bin

    # Go through the dynamic spectra and add them to the final output array. Also make a pulsestack
    fig, ax = plt.subplots(figsize=(8,12))
    i = 0
    yticks = []
    ylabels = []
    for dat in dats:

        I = StokesI_ds(dat)
        phase_bins = dat['PHASE_BINS']

        ax.plot(phase_bins*bin_size, np.nanmean(I, axis=1)+i, 'k', alpha=0.5)
        yticks.append(i)
        ylabels.append(dat['OBSID'])
        i += 1

        # All of the dynamic spectra I'm considering in this run have the same sample time, which is also what the bin size has been set to. This means that each output bin will only ever have one pixel from each dynamic spectrum added to it. If this assumption breaks, then the following method of adding each dynamic spectrum into the output dynamic spectrum won't work.
        start = phase_bins[0] - min_phase_bin
        end = start + I.shape[0]
        output_ds[start:end] = np.nansum([output_ds[start:end], I], axis=0)
        counts[start:end] += ~np.isnan(counts[start:end])

        #plt.pcolormesh(I.T)
        #plt.show()
        #print(phase_bins)

    mean_ds = output_ds / counts
    t = bin_size.to('s') * (np.arange(nphase_bins) + min_phase_bin)

    # ---- TESTING cost_function() ----
    #print(cost_function(, dm, 

    # ----------- END TEST ------------

    # If requested, dedisperse the composite dynamic spectrum
    if args.dedisperse:
        f_ref = np.nanmean(f)
        mean_ds = dedisperse_ds(mean_ds, src_dm, f, f_ref, bin_size)
        counts = dedisperse_ds(counts, src_dm, f, f_ref, bin_size)
        t -= calc_dmdelay(src_dm, f_ref, np.inf * u.MHz)

    # Form the profile
    profile = np.nanmean(mean_ds, axis=1)

    # Do a bit of extra averaging in frequency
    if fscrunch_factor is not None:
        mean_ds = np.nanmean(np.reshape(mean_ds, (mean_ds.shape[0], -1, fscrunch_factor)), axis=-1)
        counts = np.nanmean(np.reshape(counts, (counts.shape[0], -1, fscrunch_factor)), axis=-1)
        fscr = np.nanmean(np.reshape(f, (-1, fscrunch_factor)), axis=-1)
    else:
        fscr = f

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ObsID")
    plt.tight_layout()
    plt.savefig(pulsestack_png)
    plt.close(fig)

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(5,8))

    axs[0].plot(t, profile)
    axs[1].pcolormesh(t.to('s').value, fscr.to('MHz'), mean_ds.T, vmin=-0.1, vmax=0.6)

    # Draw a curve representing the DM
    if args.dedisperse and args.draw_dm_sweep is not None:
        first = True
        dmdelay = -calc_dmdelay(src_dm, f, np.nanmean(f))
        # The minus sign ^^^ is because the dynamic spectrum will already be dedispersed
        for time in args.draw_dm_sweep:
            label = f"Zero DM curve" if first else None
            axs[1].plot(dmdelay + time, f, 'w', lw=2, label=label)
            first = False

    cax = axs[2].pcolormesh(t.to('s').value, fscr.to('MHz').value, counts.T, vmin=0)
    cbar = fig.colorbar(cax, ax=axs[2], orientation='horizontal')
    cbar.set_label("Number of dynamic spectra\nthat contribute to each point")
    axs[-1].set_xlabel("Time (s)")
    axs[0].set_ylabel("Flux density (Jy/beam)")
    axs[1].set_ylabel("Frequency (MHz)")
    if args.dedisperse and args.draw_dm_sweep is not None:
        axs[1].legend()
    if args.dedisperse:
        axs[0].set_title(f"Dedispersed to DM = {src_dm}")
    axs[2].set_ylabel("Frequency (MHz)")

    plt.tight_layout()

    # Change the shape of the "counts" panel to make the x-label ("Time (s)")
    # not be obscured by the color bar
    pos = axs[2].get_position() # [x0, y0, width, height]
    delta = 0.02 # inches?
    pos.y0 += delta
    axs[2].set_position(pos)

    plt.savefig(average_ds_png)

if __name__ == '__main__':
    main()
