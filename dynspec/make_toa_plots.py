import argparse
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
from matplotlib.gridspec import GridSpec

from astropy.table import QTable
from astropy.time import Time
from astropy.coordinates import EarthLocation
from timing import *
from brokenaxes import brokenaxes

def main():
    parser = argparse.ArgumentParser(description="Construct a LaTeX table from ToA data (default=toas.ecsv) and print to stdout")
    parser.add_argument('--toas_file', default="toas.ecsv", help="Choose a different input file from the default 'toas.ecsv'")
    parser.add_argument('--output_plot', help="Set output plot filename (default: display using plt.show())")

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

    # Make plots
    nrows = 4
    ncols = 3
    width_ratios = [1, 1, 7, 2]
    fig = plt.figure(figsize=(9,6), constrained_layout=True)
    gs = GridSpec(nrows, ncols+1, figure=fig, width_ratios=width_ratios, wspace=0.03, hspace=0)
    axs = np.empty((nrows, ncols), dtype=object)
    for i in range(nrows):
        for j in range(ncols):
            axs[i, j] = fig.add_subplot(gs[i, j],
                                        sharey=axs[i, 0] if j > 0 else None)
    #fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row',
    #                        gridspec_kw={'width_ratios': width_ratios,},# 'wspace': 0.03, 'hspace': 0},
    #                        figsize=(9,5), constrained_layout=True)
    #fig.subplots_adjust(right=0.8)
    toa_row, width_row, peak_flux_row, scaled_peak_flux_row = 0, 1, 2, 3

    # Colormap setup
    cmap = plt.cm.rainbow_r
    norm = PowerNorm(0.4, 150, 1400)
    point_types = {
        170:  {'label': 'MWA (170 MHz)',      'color': '#F0E422', 'marker': '.'},
        185:  {'label': 'MWA (185 MHz)',      'color': '#009E73', 'marker': '.'},
        200:  {'label': 'MWA (200 MHz)',      'color': '#0072B2', 'marker': '.'},
        813:  {'label': 'MeerKAT (813 MHz)',  'color': '#CC79A7', 'marker': 'x'},
        888:  {'label': 'ASKAP (888 MHz)',    'color': '#D55E00', 'marker': 'o', 'markerfacecolor': 'none'},
        1284: {'label': 'MeerKAT (1284 MHz)', 'color': '#000000', 'marker': 'x'},
    }

    for col in range(ncols):
        for i in range(len(table['ToA'])):
            #color = cmap(norm(freq[i].value))
            rounded_freq = int(np.round(freq[i].value))
            color = point_types[rounded_freq]['color']
            fmt = point_types[rounded_freq]['marker']
            markerfacecolor = point_types[rounded_freq].get('markerfacecolor')
            point_fmt = {
                'fmt': fmt,
                'ls': 'none',
                'capsize': 4,
                #'markersize': 2,
                'color': color,
                'ecolor': color,
                'markerfacecolor': markerfacecolor,
            }
            #axs[fluence_row, col].errorbar(
            #    table['ToA'][i],
            #    table['fluence'][i],# / (freq[i]/u.GHz).decompose()**(-3.1),
            #    yerr=table['fluence_err'][i],# / (freq[i]/u.GHz).decompose()**(-3.1),
            #    **point_fmt
            #)
            axs[scaled_peak_flux_row, col].errorbar(
                table['ToA'][i],
                table['fitted_peak_flux_density'][i] / (freq[i]/u.GHz).decompose()**(-2.0),
                yerr=table['fitted_peak_flux_density_err'][i] / (freq[i]/u.GHz).decompose()**(-2.0),
                **point_fmt
            )
            axs[width_row, col].errorbar(
                table['ToA'][i],
                table['width'][i],
                yerr=table['width_err'][i],
                **point_fmt
            )
            axs[peak_flux_row, col].errorbar(
                table['ToA'][i],
                table['fitted_peak_flux_density'][i],# / (freq[i]/u.GHz).decompose()**(-3.1),
                yerr=table['fitted_peak_flux_density_err'][i],# / (freq[i]/u.GHz).decompose()**(-3.1),
                **point_fmt
            )
            axs[toa_row, col].errorbar(
                table['ToA'][i],
                phases[i]*ephem['period'].to('s'),
                yerr=table['ToA_err'][i].to('s'),
                **point_fmt
            )
        axs[toa_row, col].axhline(0, ls='--', color='k', alpha=0.2)

    custom_lines = [Line2D([0], [0], **v) for v in point_types.values()]
    fig.legend(handles=custom_lines, loc='center right', bbox_to_anchor=(1.0, 0.19))

    #axs[fluence_row, 0].set_ylabel(f"Fluence ({table['fluence'].unit})")
    axs[scaled_peak_flux_row, 0].set_ylabel(f"$S_{{\\rm peak,1 GHz}}$ ({table['fitted_peak_flux_density'].unit})")
    axs[width_row, 0].set_ylabel(f"$\\sigma$ ({table['width'].unit})")
    axs[peak_flux_row, 0].set_ylabel(f"$S_{{\\rm peak}}$ ({table['fitted_peak_flux_density'].unit})")
    axs[toa_row, 0].set_ylabel(f"Timing residual (s)")
    #axs[fluence_row, 0].set_yscale('log')
    axs[scaled_peak_flux_row, 0].set_yscale('log')
    axs[peak_flux_row, 0].set_yscale('log')

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
    #fig.text(0.5, 0.04, 'Time (MJD)', ha='center')
    #fig.subplots_adjust(bottom=0.15)
    fig.supxlabel('Time (MJD)                       ', y=0.0, ha='center', fontsize=10)

    # Add secondary y-axis for residuals
    def forward(y): return y/ephem['period'].to('s').value
    def inverse(y): return y*ephem['period'].to('s').value
    axs[0, -1].tick_params(axis='y', which='both', right=False, labelright=False)
    secax = axs[0, -1].secondary_yaxis('right', functions=(forward, inverse))
    secax.set_ylabel("Rotation phase")

    # Set xlims for each time chunk
    for i in range(nrows):
        axs[i, 0].set_xlim([59964, 59967])
        axs[i, 1].set_xlim([60035, 60105])
        axs[i, 2].set_xlim([60475, 60610])

    # Set custom xticks for clarity
    for i in range(nrows-1):
        for j in range(ncols):
            axs[i, j].set_xticks([])
    axs[-1, 0].set_xticks([59964, 59966])
    axs[-1, 1].set_xticks([60040, 60070, 60100])
    axs[-1, 0].ticklabel_format(axis='x', useOffset=False, style='plain')
    axs[-1, 1].ticklabel_format(axis='x', useOffset=False, style='plain')

    # And rotate them a bit so that they don't collide
    axs[-1, 1].tick_params(axis='x', labelrotation=35)
    axs[-1, 0].tick_params(axis='x', labelrotation=35)

    # Add colorbar manually
    #from matplotlib.cm import ScalarMappable
    #sm = ScalarMappable(cmap=cmap, norm=norm)
    #sm.set_array([])
    #cbar = fig.colorbar(sm, ax=axs[:, -1])
    #cbar.set_label(f'Frequency ({freq.unit})')

    #formatter = ScalarFormatter(useOffset=True, useMathText=False)
    #formatter.set_powerlimits((-1000, 1000))  # Prevent scientific notation
    #axs[-1,0].xaxis.set_major_formatter(formatter)

    if args.output_plot is not None:
        #plt.tight_layout()
        plt.savefig(args.output_plot)
    else:
        plt.show()

if __name__ == '__main__':
    main()
