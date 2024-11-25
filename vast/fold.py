import csv
import numpy as np
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt

PEPOCH = Time(59965.03767627493, scale='utc', format='mjd')
P = 4186.32874813198 * u.s

# TODO: Add DM and Barycentre corrections

with open('VAST-Pipeline-Source-ID-34260695.csv') as csvfile:
    myreader = csv.reader(csvfile)
    # Skip the header
    next(myreader)

    rows = [row for row in myreader]

    #print([row[2][:22] for row in myreader])
    start_times = Time([row[2].split('+')[0] for row in rows], scale='utc', format='isot')
    is_detection = [row[16] == "true" for row in rows]

start_pulses, start_phases = np.divmod(((start_times - PEPOCH) / P).decompose(), 1)
obs_duration = 11*u.min
end_phases   = start_phases + (obs_duration / P).decompose()

# Fix the start_pulses so that things that start at phases > 0.5 get added 1 to it
start_pulses += start_phases >= 0.5
start_pulses = start_pulses.astype(int)

centre_phases = 0.5*(start_phases + end_phases)
xerr = 0.5*(end_phases - start_phases)

pulse_width = 2*u.min # Approximate only
pulse_width_phase = (pulse_width / P).decompose()

ids = np.arange(len(centre_phases))

fig, ax = plt.subplots()
ax.errorbar(centre_phases, ids, xerr=xerr, fmt='none', capsize=5)
ax.errorbar(centre_phases - 1, ids, xerr=xerr, fmt='none', capsize=5)
y_min, y_max = ax.get_ylim()
ax.fill_betweenx(y=[y_min-1, y_max+1], x1=-pulse_width_phase/2, x2=pulse_width_phase/2, color='lightcoral', alpha=0.7)
ax.set_xlim([-0.5, 0.5])
ax.set_ylim([y_min, y_max])
ax.set_yticks(ids)
ax.set_yticklabels([str(pulse) for pulse in start_pulses])
ax.set_xlabel('Pulse phase')
ax.set_ylabel('Pulse number')

ax2 = ax.twinx()
ax2.set_yticks(ids)
ax2.set_ylim([y_min, y_max])
ax2.set_yticklabels([t.isot for t in start_times])

plt.tight_layout()
plt.savefig('fold.png')

