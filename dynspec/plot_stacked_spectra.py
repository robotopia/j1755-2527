import numpy as np
import matplotlib.pyplot as plt

data = [
    {'pkl': 'D0042_observations_185MHz_stacked.pkl', 'fscrunch_factor': 12},
    {'pkl': 'D0042_observations_200MHz_stacked.pkl', 'fscrunch_factor': 12},
    {'pkl': 'GPM_observations_stacked.pkl',          'fscrunch_factor': 12},
    {'pkl': '2023_meerkat_stacked.pkl',              'fscrunch_factor': 32},
]

fmin = np.inf
fmax = -np.inf

for datum in data:
    pkl = datum['pkl']
    fscrunch_factor = datum['fscrunch_factor']
    dat = np.load(pkl, allow_pickle=True)
    I = np.nanmean(np.reshape(dat['I'], (dat['I'].shape[0], -1, fscrunch_factor)), axis=-1)
    f = np.nanmean(np.reshape(dat['f'], (-1, fscrunch_factor)), axis=-1)/1e6
    t = dat['t']
    plt.pcolormesh(t, f, I.T, alpha=0.8)
    if min(f) < fmin:
        fmin = min(f)
    if max(f) > fmax:
        fmax = max(f)

# Plot didespersion sweep
ffine = np.linspace(fmin, fmax, 1000)
DM = 710
plt.plot(4148.808*DM/ffine**2, ffine, 'r')
plt.yscale('log')
plt.show()
