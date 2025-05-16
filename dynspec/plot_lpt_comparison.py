import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

lpts = [
    {
        'name': 'ASKAP J1832-0911',
        'min_width': 0 * u.s,  # TODO
        'max_width': 0 * u.s,  # TODO
        'period': 2656.24 * u.s,
    },
    {
        'name': 'ASKAP J1839-0756',
        'min_width': 320 * u.s,
        'max_width': 710 * u.s,
        'period': 23221.29 * u.s,
    },
    {
        'name': 'ASKAP J1935+2148',
        'min_width': 10 * u.s,
        'max_width': 50 * u.s,
        'period': 3225.313 * u.s,
    },
    {
        'name': 'CHIME J0630+25',
        'min_width': 0.2 * u.s,
        'max_width': 0.91 * u.s,
        'period': 421.35542 * u.s,
    },
    {
        'name': 'GCRT J1745-3009',
        'min_width': 600 * u.s,
        'max_width': 600 * u.s,
        'period': 4620.720 * u.s,  # Remember to cite Spreeuw et al. (2009) for this number.
    },
    {
        'name': 'GLEAM-X J1627-5235',
        'min_width': 30 * u.s,
        'max_width': 60 * u.s,
        'period': 1091.169 * u.s,
    },
    {
        'name': 'GPM J1839-10',
        'min_width': 30 * u.s,
        'max_width': 300 * u.s,
        'period': 1318.1957 * u.s,
    },
    {
        'name': 'ILT J1101+5521',
        'min_width': 30 * u.s,
        'max_width': 90 * u.s,
        'period': 7531.17 * u.s,
    },
]

periods = [lpt['period'].to('s').value for lpt in lpts] * u.s
min_duty_cycles = np.array([(lpt['min_width']/lpt['period']).decompose() for lpt in lpts])
max_duty_cycles = np.array([(lpt['max_width']/lpt['period']).decompose() for lpt in lpts])
duty_cycles = 0.5*(max_duty_cycles + min_duty_cycles)
duty_cycles_err = 0.5*(max_duty_cycles - min_duty_cycles)

plt.errorbar(periods.to('min').value, duty_cycles, yerr=duty_cycles_err, capsize=2, fmt='none', ls='none')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Period (min)")
plt.ylabel("Duty cycle")
plt.show()
