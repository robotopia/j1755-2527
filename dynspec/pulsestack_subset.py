import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import colormaps, colors, cm
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys
import argparse

from timing import *

def main():
    parser = argparse.ArgumentParser(description="Plot a pulsestack from the given (pickle) dynamic spectra files")

    parser.add_argument('ds_files', nargs='*', help="The names of the input .pkl files, without the .pkl extension")
    parser.add_argument('--output_pulsestack_image', help="The filename to use for the pulsestack image")
    parser.add_argument('--nobary', action='store_true', help="Do not apply a barycentric correction")
    parser.add_argument('--colormap', default='rainbow_r', help="Apply this colormap for frequencies (default='rainbow_r')")
    parser.add_argument('--vmin', type=float, default=150.0, help="Define lowest frequency (MHz) on colormap (default=150)")
    parser.add_argument('--vmax', type=float, default=1500.0, help="Define highest frequency (MHz) on colormap (default=1500)")
    parser.add_argument('--xlim', type=float, nargs=2, help="Define limits of x-axis")
    parser.add_argument('--grid', action='store_true', help="Toggle on grid on plot")
    parser.add_argument('--ncols', type=int, default=1, help="Number of columns in output pulsestack")

    args = parser.parse_args()

    ephemeris = get_J1755_ephemeris()

    fig, axs = plt.subplots(ncols=args.ncols, figsize=(10,10), squeeze=False)
    yticks = [[] for i in range(args.ncols)]
    ylabels = [[] for i in range(args.ncols)]
    #cmap = colormaps[args.colormap]
    #cmap_norm = lambda x : ((x - args.vmin*u.MHz)/(args.vmax*u.MHz - args.vmin*u.MHz)).decompose()

    pkls = args.ds_files

    # Load data, and convert the times to phases and phase bins
    ys = [0 for i in range(args.ncols)]
    pulses_per_col = int(np.ceil(len(pkls) / args.ncols))
    for i in range(len(pkls)):
        pkl = pkls[i]
        print(f"Opening {pkl}...")
        col = i // pulses_per_col
        dat = np.load(pkl, allow_pickle=True)
        dat['OBSID'] = pkl[:-4].replace('_', ' ')

        I, fmask, _ = Stokes_ds(dat, pol='I')

        f = dat['FREQS'] * u.Hz
        #color = cmap(cmap_norm(np.mean(f)))

        location = EarthLocation.of_site(dat['TELESCOPE'])
        times = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=location)
        dt = (times[1] - times[0]).to('s')
        if dt < 3*u.s:
            times = times[1::2] - dt/2
            dt *= 2
            t_len = len(times)
            I = np.nanmean(np.reshape(I[:t_len*2,:], (-1, 2, I.shape[1])), axis=1)

        # Dedispersion requires no nans
        I[np.isnan(I)] = 0.0

        # Dedisperse
        f_ref = np.mean(f)
        Idd = dedisperse_ds(I, ephemeris['DM'], f, f_ref, dt)
        times -= calc_dmdelay(ephemeris['DM'], f_ref, np.inf*u.Hz)

        location = EarthLocation.of_site(dat['TELESCOPE'])
        if args.nobary == False: # i.e. if we DO want to barycentre (which is the default)
            times = barycentre(times, ephemeris['coord'])

        pulses, phases = fold(times, ephemeris['period'], ephemeris['PEPOCH'])

        # Reapply channel flags
        if len(fmask) > 0:
            Idd[:,fmask] = np.nan

        # Form lightcurve
        lightcurve = np.nanmean(Idd, axis=1)

        # Fit a scattered pulse
        τ = (ephemeris['tau_sc']/ephemeris['period'] * (f_ref/u.GHz)**(-4)).decompose()
        def model(ph, A, μ, m, c, σ):
            return scattered_pulse_model(ph, A, μ, m, c, σ, τ)
        peak_idx = np.nanargmax(lightcurve)
        p0 = [lightcurve[peak_idx], phases[peak_idx], 0, 0, 20/ephemeris['period'].to('s').value]
        #def model(ph, m, c, A, mu, sigma):
        #    return m*ph + c + A*np.exp(-0.5*(ph-mu)**2/sigma**2)
        #p0 = [0, 0, np.max(lightcurve), 0, 10/ephemeris['period'].to('s').value]
        popt, _ = curve_fit(model, phases, lightcurve, p0=p0)
        A, μ, m, c, σ = popt
        #print(f'{pkl}: {m = }, {c = }')
        lightcurve -= m*phases + c

        if len(set(pulses)) != 1:
            phase_0_idx = np.argmin(np.abs(phases))
            pulse_at_phase_0 = pulses[phase_0_idx]
            mask = (pulses == pulse_at_phase_0)
            phases = phases[mask]
            lightcurve = lightcurve[mask]

        ph_fine = np.linspace(phases[0], phases[-1], 1000)

        y = pulse_model(ph_fine, A, μ, 0, 0, σ)
        y1 = scattered_pulse_model(ph_fine, A, μ, 0, 0, σ, τ)

        ynorm = y / np.nanmax(y)
        y1norm = y1 / np.nanmax(y1)

        axs[0][col].plot(phases*ephemeris['period'], 0.9*lightcurve/np.max(y1) + ys[col], 'k', lw=1, alpha=0.7)
        #axs[0][col].plot(phases*ephemeris['period'], 0.9*lightcurve/np.max(lightcurve) + ys[col], 'k', lw=0.5, alpha=0.5)
        axs[0][col].plot(ph_fine*ephemeris['period'], 0.9*ynorm + ys[col], 'b', lw=2, alpha=0.3)
        axs[0][col].plot(ph_fine*ephemeris['period'], 0.9*y1norm + ys[col], 'r', lw=2, alpha=0.3)

        yticks[col].append(ys[col])
        telescope = dat['TELESCOPE'] if dat['TELESCOPE'] != "Siding Spring Observatory" else "ATCA"
        ylabels[col].append(f"{times[0].iso[:10]}\n{telescope}\n{f_ref.to('MHz'):.0f}")#\n(Pulse #{int(np.round(np.median(pulses)))})")
        ys[col] += 1

    for col in range(args.ncols):
        axs[0][col].set_yticks(yticks[col])
        axs[0][col].set_yticklabels(ylabels[col], fontsize=16)
        axs[0][col].set_xlabel("Time from predicted ToA (s)", fontsize=20)
        if args.xlim:
            axs[0][col].set_xlim(args.xlim)
        axs[0][col].set_ylim([-0.2, ys[col]+0.1])
        axs[0][col].axvline(0, ls='--', alpha=0.5, c='g', zorder=-100)
        axs[0][col].tick_params(axis='x', labelsize=16)
    #axs[0][0].set_ylabel("Observation start time (UTC)", fontsize=20)

    if args.grid:
        plt.grid(axis='x')

    #fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=args.vmin, vmax=args.vmax), cmap=cmap),
    #             ax=ax, label='Frequency (MHz)')
    plt.tight_layout()

    if args.output_pulsestack_image is not None:
        plt.savefig(args.output_pulsestack_image)
    else:
        plt.show()

if __name__ == '__main__':
    main()
