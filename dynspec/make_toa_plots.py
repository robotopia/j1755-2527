import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
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
    width_ratios = [1, 1, 8]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row',
                            gridspec_kw={'width_ratios': width_ratios, 'wspace': 0.03, 'hspace': 0},
                            figsize=(10,10))
    toa_row, fluence_row, width_row = 0, 1, 2

    # Colormap setup
    cmap = plt.cm.rainbow_r
    norm = PowerNorm(0.3, freq.min().value, 1400)

    for col in range(ncols):
        for i in range(len(table['ToA'])):
            color = cmap(norm(freq[i].value))
            point_fmt = {'fmt': 'o', 'ls': 'none', 'capsize': 4, 'markersize': 2, 'color': color, 'ecolor': color}
            axs[fluence_row, col].errorbar(table['ToA'][i], table['fluence'][i], yerr=table['fluence_err'][i], **point_fmt)
            axs[width_row, col].errorbar(table['ToA'][i], table['width'][i], yerr=table['width_err'][i], **point_fmt)
            axs[toa_row, col].errorbar(table['ToA'][i], phases[i]*ephem['period'].to('s'), yerr=table['ToA_err'][i].to('s'), **point_fmt)
        axs[toa_row, col].axhline(0, ls='--', color='k', alpha=0.2)

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
    fig.text(0.5, 0.04, 'Time (MJD)', ha='center')

    # Add secondary y-axis for residuals
    def forward(y): return y/ephem['period'].to('s').value
    def inverse(y): return y*ephem['period'].to('s').value
    axs[0, -1].tick_params(axis='y', which='both', right=False, labelright=False)
    secax = axs[0, -1].secondary_yaxis('right', functions=(forward, inverse))
    secax.set_ylabel("Rotation phase")

    # Set xlims for each time chunk
    axs[0, 0].set_xlim([59964, 59967])
    axs[0, 1].set_xlim([60092.85, 60093.1])
    axs[0, 2].set_xlim([60475, 60610])

    # Add colorbar manually
    from matplotlib.cm import ScalarMappable
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axs[:, -1])
    cbar.set_label(f'Frequency ({freq.unit})')

    plt.tight_layout()
    plt.savefig(args.output_plot)
    #plt.show()

if __name__ == '__main__':
    main()
