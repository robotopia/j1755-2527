import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps, colors, cm
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys
import argparse

from timing import fold, barycentre, calc_dmdelay, StokesI_ds, dedisperse_ds, get_J1755_ephemeris

def main():
    parser = argparse.ArgumentParser(description="Plot a pulsestack from the given (pickle) dynamic spectra files")

    parser.add_argument('ds_files', nargs='*', help="The names of the input .pkl files, without the .pkl extension")
    parser.add_argument('--output_pulsestack_image', help="The filename to use for the pulsestack image")
    parser.add_argument('--nobary', action='store_true', help="Do not apply a barycentric correction")
    parser.add_argument('--colormap', default='rainbow_r', help="Apply this colormap for frequencies (default='rainbow_r')")
    parser.add_argument('--vmin', type=float, default=150.0, help="Define lowest frequency (MHz) on colormap (default=150)")
    parser.add_argument('--vmax', type=float, default=1500.0, help="Define highest frequency (MHz) on colormap (default=1500)")
    parser.add_argument('--grid', action='store_true', help="Toggle on grid on plot")

    args = parser.parse_args()

    ephemeris = get_J1755_ephemeris()

    fig, ax = plt.subplots(figsize=(9,12))
    yticks = []
    ylabels = []
    cmap = colormaps[args.colormap]
    cmap_norm = lambda x : ((x - args.vmin*u.MHz)/(args.vmax*u.MHz - args.vmin*u.MHz)).decompose()

    pkls = args.ds_files

    # Load data, and convert the times to phases and phase bins
    for i in range(len(pkls)):
        pkl = pkls[i]
        print(f"Opening {pkl}...")
        dat = np.load(pkl, allow_pickle=True)
        dat['OBSID'] = pkl[:-4].replace('_', ' ')

        I = StokesI_ds(dat)

        f = dat['FREQS'] * u.Hz
        color = cmap(cmap_norm(np.mean(f)))

        location = EarthLocation.of_site(dat['TELESCOPE'])
        times = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=location)
        dt = (times[1] - times[0]).to('s')

        # Dedisperse
        f_ref = np.mean(f)
        Idd = dedisperse_ds(I, ephemeris['DM'], f, f_ref, dt)
        times -= calc_dmdelay(ephemeris['DM'], f_ref, np.inf*u.Hz)

        location = EarthLocation.of_site(dat['TELESCOPE'])
        if args.nobary == False: # i.e. if we DO want to barycentre (which is the default)
            times = barycentre(times, ephemeris['coord'])

        _, phases = fold(times, ephemeris['period'], ephemeris['PEPOCH'])

        ax.plot(phases*ephemeris['period'], np.nanmean(I, axis=1)+i, 'k', color=color)
        yticks.append(i)
        ylabels.append(dat['OBSID'])
        i += 1
    
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ObsID")

    if args.grid:
        plt.grid(axis='x')

    fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=args.vmin, vmax=args.vmax), cmap=cmap),
                 ax=ax, label='Frequency (MHz)')
    plt.tight_layout()

    if args.output_pulsestack_image is not None:
        plt.savefig(args.output_pulsestack_image)
    else:
        plt.show()

if __name__ == '__main__':
    main()
