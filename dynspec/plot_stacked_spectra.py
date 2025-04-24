import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from timing import *

ephem = get_J1755_ephemeris()
DM = ephem['DM'].value
tau_sc_1GHz = ephem['tau_sc'].value

def calc_dm_delay(DM, f):
    '''
    DM in pc/cm³
    f in MHz
    result is in seconds
    '''
    return 4148.808 * DM / f**2

data = [
    {'pkl': 'D0042_observations_185MHz_stacked.pkl', 'fscrunch_factor': 12,
     'tscrunch_factor': 1, 'color': 'red', 'label': 'D0042 (185 MHz)'},
    {'pkl': 'D0042_observations_200MHz_stacked.pkl', 'fscrunch_factor': 12,
     'tscrunch_factor': 1, 'color': 'orange', 'label': 'D0042 (200 MHz)'},
    {'pkl': 'GPM_observations_stacked.pkl',          'fscrunch_factor': 12,
     'tscrunch_factor': 1, 'color': 'yellow', 'label': 'GPM 2024'},
    {'pkl': '2023_meerkat_stacked.pkl',              'fscrunch_factor': 16,
     'tscrunch_factor': 2, 'color': 'green', 'label': 'MeerKAT (2023)'},
    {'pkl': '2024_meerkat_stacked.pkl',              'fscrunch_factor': 33,
     'tscrunch_factor': 1, 'color': 'blue', 'label': 'MeerKAT (2024)'},
]

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(5,8))
rects = []

for datum in data:
    pkl = datum['pkl']
    fscrunch_factor = datum['fscrunch_factor']
    tscrunch_factor = datum['tscrunch_factor']
    dat = np.load(pkl, allow_pickle=True)
    t_len = (dat['I'].shape[0]//tscrunch_factor)*tscrunch_factor
    I = np.nanmean(np.reshape(dat['I'][:t_len,:], (dat['I'].shape[0]//tscrunch_factor, tscrunch_factor, -1, fscrunch_factor)), axis=(1,3))
    f = np.nanmean(np.reshape(dat['f'], (-1, fscrunch_factor)), axis=-1)/1e6
    t = np.nanmean(np.reshape(dat['t'][:t_len], (-1, tscrunch_factor)), axis=-1)
    ax = axs[0] if 'meerkat' in pkl else axs[1]
    ax.pcolormesh(t, f, I.T, alpha=0.8, rasterized=True)
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

# Take the MeerKAT 2024 as a template, and produce scattered pulses at different frequencies
dat = np.load('2024_meerkat_stacked.pkl', allow_pickle=True)
I = dat['I']
t = dat['t']
dt = t[1] - t[0]
f = dat['f']/1e6 # to MHz
dm_delays = calc_dm_delay(DM, f) # in seconds
dm_delays_bins = np.round(dm_delays / dt).astype(int)
Idd = np.array([np.roll(row, -shift) for row, shift in zip(I.T, dm_delays_bins)]).T

t = t[10:-10]
lc = np.nanmean(Idd, axis=1)[10:-10] + 0.011 + 0.00001*t
fref = np.nanmean(f)

def plot_scattered_pulse(ax, f, height):
    dm_delay = calc_dm_delay(DM, f)
    if f == fref:
        lc_sc = lc
    else:
        tau_sc = tau_sc_1GHz * (f/1e3)**-4
        kernel_t = (np.arange(201) - 100)*dt
        kernel = np.exp(-kernel_t / tau_sc)
        kernel[kernel_t < 0] = 0.0
        lc_sc = np.convolve(lc, kernel, mode='same')
    lc_f = lc_sc/np.max(lc_sc)*height + f - height/2
    ax.plot(t + dm_delay, lc_f, c='white', lw=2, zorder=200)

plot_scattered_pulse(axs[0], fref, 100)
for f in [210, 195, 175]:
    plot_scattered_pulse(axs[1], f, 7)

# Plot didespersion sweep
ffine = np.linspace(570, 1710, 1000)
axs[0].plot(calc_dm_delay(DM, ffine), ffine, c='cyan', ls='--')
axs[0].plot(calc_dm_delay(DM, ffine) - 125, ffine, c='cyan', ls='--')

ffine = np.linspace(169.96, 215.5, 1000)
axs[1].plot(calc_dm_delay(DM, ffine), ffine, c='cyan', ls='--')
axs[1].plot(calc_dm_delay(DM, ffine) - 125, ffine, c='cyan', ls='--')

for ax in axs:
    ax.set_xlim([-230, 330]) # Preferable to sharex if you want the tickmarks repeated on both axes

axs[0].set_ylim([555, 1725])
axs[1].set_ylim([169, 216.2])
axs[0].set_ylabel("Frequency (MHz)")
axs[1].set_ylabel("Frequency (MHz)")

plt.xlabel("Time (s)")

axs[0].legend(rects, [rect.get_label() for rect in rects], loc='upper right')

plt.tight_layout()
plt.savefig('stacked_spectra.pdf')
