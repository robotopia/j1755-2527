import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.time import Time
from astropy.coordinates import EarthLocation
from timing import *
from brokenaxes import brokenaxes

def main():
    parser = argparse.ArgumentParser(description="Construct a LaTeX table from ToA data (default=toas.ecsv) and print to stdout")
    parser.add_argument('--toas_file', default="toas.ecsv", help="Choose a different input file from the default 'toas.ecsv")
    parser.add_argument('--output_plot', default="toa_details.pdf", help="Set output plot filename (default='toa_details.pdf')")

    args = parser.parse_args()

    # Load the table
    table = QTable.read(args.toas_file)

    # Calculate timing residuals
    ephem = get_J1755_ephemeris()
    freq = table['freq']
    locations = [EarthLocation.of_site(table['telescope'][i]) for i in range(len(table['ToA']))]
    toas = Time([Time(table['ToA'][i], scale='utc', format='mjd', location=locations[i]) for i in range(len(table['ToA']))])
    bary_toas = barycentre(toas, ephem['coord'])
    bary_toas_dd = bary_toas - calc_dmdelay(ephem['DM'], freq, np.inf)
    pulses, phases = fold(bary_toas_dd, ephem['period'], ephem['PEPOCH'])

    # Set the 2023 MeerKAT fluences to None --- we don't trust them!
    mask = (table['telescope'] == 'MeerKAT') & (np.array([d.isocalendar().year for d in toas.to_datetime()]) == 2023)
    for row_number in np.where(mask)[0]:
        table[row_number]['fluence'] = np.nan

    # Make plots
    nrows = 3
    ncols = 3
    width_ratios = [1, 3, 5]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row',
                            gridspec_kw={'width_ratios': [0.1, 0.1, 0.8], 'wspace': 0.05, 'hspace': 0},
                            figsize=(10,10))
    toa_row, fluence_row, width_row = 0, 1, 2

    point_fmt = {'fmt': 'o', 'ls': 'none', 'capsize': 4, 'markersize': 2}
    for col in range(ncols):
        axs[fluence_row, col].errorbar(table['ToA'], table['fluence'], yerr=table['fluence_err'], **point_fmt)
        axs[width_row, col].errorbar(table['ToA'], table['width'], yerr=table['width_err'], **point_fmt)
        axs[toa_row, col].errorbar(table['ToA'], phases*ephem['period'].to('s'), yerr=table['ToA_err'].to('s'), **point_fmt)
    axs[fluence_row, 0].set_ylabel(f"Fluence ({table['fluence'].unit})")
    axs[width_row, 0].set_ylabel(f"$\\sigma$ ({table['width'].unit})")
    axs[toa_row, 0].set_ylabel(f"Timing residual (s)")
    axs[fluence_row, 0].set_yscale('log')

    # Add diagonal "break marks" between subplots
    # Hide inner spines and ticks, and draw break marks
    d = 0.015  # size of break marker
    for i in range(nrows):
        for j in range(ncols-1):
            axs[i, j].spines['right'].set_visible(False)
            axs[i, j].tick_params(labelright=False)

            break_mark_args = dict(transform=axs[i,j].transAxes, color='k', clip_on=False)
            break_mark_args.update(transform=axs[i,j].transAxes)
            axs[i, j].plot((1-2*d/width_ratios[j], 1+2*d/width_ratios[j]), (1-d, 1+d), **break_mark_args)
            axs[i, j].plot((1-2*d/width_ratios[j], 1+2*d/width_ratios[j]), (-d, d), **break_mark_args)
        for j in range(1,ncols):
            axs[i, j].spines['left'].set_visible(False)
            axs[i, j].tick_params(labelleft=False)

            break_mark_args = dict(transform=axs[i,j].transAxes, color='k', clip_on=False)
            break_mark_args.update(transform=axs[i,j].transAxes)
            axs[i, j].plot((-2*d/width_ratios[j], 2*d/width_ratios[j]), (1-d, 1+d), **break_mark_args)
            axs[i, j].plot((-2*d/width_ratios[j], 2*d/width_ratios[j]), (-d, d), **break_mark_args)
            if j != ncols - 1:
                axs[i,j].tick_params(axis='y', which='both', left=False, labelleft=False)
                #axs[i, j].set_yticks([])
            else:
                axs[i, j].yaxis.tick_right()

    # Add shared x label
    fig.text(0.5, 0.04, 'Time (x)', ha='center')

    axs[0, 0].set_xlim([59964, 59967])
    axs[0, 1].set_xlim([60092.85, 60093.1])
    axs[0, 2].set_xlim([60475, 60610])

    plt.tight_layout()
    plt.savefig(args.output_plot)
    #plt.show()

if __name__ == '__main__':
    main()
