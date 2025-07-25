import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerTuple
import astropy.units as u

point_types = {
    'Binaries hosting WDs': {
        'color': 'r',
        'ecolor': 'r',
        'fmt': 's',
        'fillstyle': 'none',
    },
    'Unknown': {
        'color': 'k',
        'ecolor': 'k',
        'fmt': 'x',
        'fillstyle': 'none',
    },
    'AR Sco-like': {
        'color': 'b',
        'ecolor': 'b',
        'fmt': 'o',
        'fillstyle': 'full',
    },
}

lpts = [
    {
        'name': 'ASKAP J1755-2527',
        'type': 'Unknown',
        'min_width': 40 * u.s, # Be sure to note: these are W10, unscattered
        'max_width': 172 * u.s,
        'period': 4186.32827 * u.s,
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
            'fontweight': 'bold',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        },
    },
    {
        'name': 'ASKAP J1832-0911',
        'type': 'Unknown',
        'min_width': 2 * u.min,
        'max_width': 4 * u.min,
        'period': 2656.24 * u.s,
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        }
    },
    {
        'name': 'ASKAP J1839-0756',
        'type': 'Unknown',
        'min_width': 320 * u.s,
        'max_width': 710 * u.s,
        'period': 23221.29 * u.s,
        'text_kwargs': {
            'ha': 'right',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        },
        'voffset': 1,
    },
    {
        'name': 'ASKAP J1935+2148',
        'type': 'Unknown',
        'min_width': 10 * u.s,
        'max_width': 50 * u.s,
        'period': 3225.313 * u.s,
        'text_kwargs': {
            'ha': 'right',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        }
    },
    {
        'name': 'CHIME J0630+25',
        'type': 'Unknown',
        'min_width': 0.2 * u.s,
        'max_width': 0.91 * u.s,
        'period': 421.35542 * u.s,
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        }
    },
    {
        'name': 'GCRT J1745-3009',
        'type': 'Unknown',
        'min_width': 600 * u.s,
        'max_width': 600 * u.s,
        'period': 4620.720 * u.s,  # Remember to cite Spreeuw et al. (2009) for this number.
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 0,
        }
    },
    {
        'name': 'GLEAM-X J0704-37',
        'type': 'Binaries hosting WDs',
        'min_width': 30 * u.s,
        'max_width': 80 * u.s,
        'period': 10496.6 * u.s,
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        }
    },
    {
        'name': 'CHIME/ILT J1634+44',
        'type': 'Binaries hosting WDs',
        'min_width': 0.07 * u.s,
        'max_width': 0.8 * u.s,
        'period': 841.245895 * u.s,
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        }
    },
    {
        'name': 'GLEAM-X J1627-5235',
        'type': 'Unknown',
        'min_width': 30 * u.s,
        'max_width': 60 * u.s,
        'period': 1091.169 * u.s,
        'text_kwargs': {
            'ha': 'right',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        },
        'voffset': 1,
    },
    {
        'name': 'GPM J1839-10',
        'type': 'Unknown',
        'min_width': 30 * u.s,
        'max_width': 300 * u.s,
        'period': 1318.1957 * u.s,
        'text_kwargs': {
            'ha': 'right',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        }
    },
    {
        'name': 'ILT J1101+5521',
        'type': 'Binaries hosting WDs',
        'min_width': 30 * u.s,
        'max_width': 90 * u.s,
        'period': 7531.17 * u.s,
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 2,
        }
    },
    {
        'name': 'J1912-4410',
        'type': 'AR Sco-like',
        'min_width': 16 * u.s, # TODO: Get Csanad to give me better numbers from data
        'max_width': 16 * u.s,
        'period': 319.34903 * u.s,
        'text_kwargs': {
            'ha': 'left',
            'va': 'center',
        },
        'errorbar_kwargs': {
            'ls': 'none',
            'capsize': 0,
        }
    },
    #{
    #    'name': 'AR Sco',
    #    'type': 'AR Sco-like',
    #    'min_width': 0 * u.s, # TODO: Get numbers for this!!
    #    'max_width': 1.97/2 * u.min,
    #    'period': 1.97 * u.min,
    #    'text_kwargs': {
    #        'ha': 'left',
    #        'va': 'center',
    #    },
    #    'errorbar_kwargs': {
    #        'ls': 'none',
    #        'capsize': 2,
    #    }
    #},
]

period_units = 'min'

plt.figure(figsize=(8, 3))

for lpt in lpts:
    min_duty_cycle = (lpt['min_width']/lpt['period']).decompose()
    max_duty_cycle = (lpt['max_width']/lpt['period']).decompose()
    duty_cycle = 0.5*(max_duty_cycle + min_duty_cycle)
    duty_cycle_err = 0.5*(max_duty_cycle - min_duty_cycle)

    plt.errorbar(
        lpt['period'].to(period_units).value,
        duty_cycle,
        yerr=duty_cycle_err,
        **lpt['errorbar_kwargs'],
        **point_types[lpt['type']],
    )

    lpt_label = " " + lpt['name'] if lpt['text_kwargs']['ha'] == 'left' else lpt['name'] + " "
    if 'voffset' in lpt.keys():
        lpt_label = '\n'*lpt['voffset'] + lpt_label
    plt.annotate(lpt_label,
                 (lpt['period'].to(period_units).value, duty_cycle + 0.5*duty_cycle_err),
                 **lpt['text_kwargs'], alpha=0.7)

# Construct custom legend
custom_lines = []
for point_type, kwargs in point_types.items():
    custom_lines.append(mlines.Line2D(
        [], [],
        color=kwargs['color'],
        marker=kwargs['fmt'],
        fillstyle=kwargs['fillstyle'],
        #linestyle=('--' if point_type == 'Binaries hosting WDs' else '-'),
        label=point_type,
    ))
plt.legend(handles=custom_lines, loc='lower right')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(f"Period ({period_units})")
plt.ylabel("Pulse duty cycle")

plt.tight_layout()
plt.savefig('lpt_comparison.pdf')
