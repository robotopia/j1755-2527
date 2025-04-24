import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from timing import *

ephem = get_J1755_ephemeris()
DM = ephem['DM'].value

data = [
    {'pkl': 'D0042_observations_185MHz_stacked.pkl', 'fscrunch_factor': 12,
     'color': 'red', 'label': 'D0042 (185 MHz)'},
    {'pkl': 'D0042_observations_200MHz_stacked.pkl', 'fscrunch_factor': 12,
     'color': 'orange', 'label': 'D0042 (200 MHz)'},
    {'pkl': 'GPM_observations_stacked.pkl',          'fscrunch_factor': 12,
     'color': 'yellow', 'label': 'GPM 2024'},
    {'pkl': '2023_meerkat_stacked.pkl',              'fscrunch_factor': 8 ,
     'color': 'green', 'label': 'MeerKAT (2023)'},
    {'pkl': '2024_meerkat_stacked.pkl',              'fscrunch_factor': 33,
     'color': 'blue', 'label': 'MeerKAT (2024)'},
]

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8,12), sharex=True)
rects = []

for datum in data:
    pkl = datum['pkl']
    fscrunch_factor = datum['fscrunch_factor']
    dat = np.load(pkl, allow_pickle=True)
    I = np.nanmean(np.reshape(dat['I'], (dat['I'].shape[0], -1, fscrunch_factor)), axis=-1)
    f = np.nanmean(np.reshape(dat['f'], (-1, fscrunch_factor)), axis=-1)/1e6
    t = dat['t']
    ax = axs[0] if 'meerkat' in pkl else axs[1]
    ax.pcolormesh(t, f, I.T, alpha=0.8)
    # Draw a rectangle around the stacked spectrum
    dt = t[1] - t[0]
    df = f[1] - f[0]
    x = t[0] - 0.5*dt
    y = f[0] - 0.5*df
    width = t[-1] - t[0] + dt
    height = f[-1] - f[0] + df
    rect = Rectangle((x, y), width, height, linewidth=2, edgecolor=datum['color'], facecolor='none', alpha=0.8, zorder=100, label=datum['label'])
    ax.add_patch(rect)
    rects.append(rect)

# Plot didespersion sweep
ffine = np.linspace(570, 1710, 1000)
axs[0].plot(4148.808*DM/ffine**2, ffine, 'w')

ffine = np.linspace(169.96, 215.5, 1000)
axs[1].plot(4148.808*DM/ffine**2, ffine, 'w')

axs[0].set_xlim([-250, 400])
axs[0].set_ylim([555, 1725])
axs[1].set_ylim([169, 216.2])
axs[0].set_ylabel("Frequency (MHz)")
axs[1].set_ylabel("Frequency (MHz)")

plt.xlabel("Time (s)")

axs[0].legend(rects, [rect.get_label() for rect in rects])

plt.tight_layout()
plt.savefig('stacked_spectra.png')
