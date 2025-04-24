import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import pickle

on = np.loadtxt('ts_on_source.txt')
off = np.loadtxt('ts_off_source.txt')
on_obsprior = np.loadtxt('ts_on_source_obsprior.txt')
off_obsprior = np.loadtxt('ts_off_source_obsprior.txt')

mean_on = np.mean(on)
mean_off = np.mean(off)
mean_on_obsprior = np.mean(on_obsprior)
mean_off_obsprior = np.mean(off_obsprior)

# From Angie:
# On source: pixel coordinates 1546, 1258, value of 1.07391
# Off source: pixel coordinates 1549, 1563, value of 0.00331
# On source (obs prior): pixel coordinates 1501, 1255, value of 0.030
# Off source (obs prior): pixel coordinates 1506, 1206, value of -0.120

mean_on_should_be = 1.07391
mean_off_should_be = 0.00331

mean_on_obsprior_should_be = 0.030
mean_off_obsprior_should_be = -0.120

mean_on_adjust = mean_on_should_be - mean_on
mean_off_adjust = mean_off_should_be - mean_off
mean_on_obsprior_adjust = mean_on_obsprior_should_be - mean_on_obsprior
mean_off_obsprior_adjust = mean_off_obsprior_should_be - mean_off_obsprior

on += mean_on_adjust
off += mean_off_adjust
on_obsprior += mean_on_obsprior_adjust
off_obsprior += mean_off_obsprior_adjust

b = 20 # Number of time samples to bin together
sample_time = 0.5 # seconds

binned_on = np.mean(np.reshape(np.append(on_obsprior, on), (-1, b)), axis=-1)
binned_off = np.mean(np.reshape(np.append(off_obsprior, off), (-1, b)), axis=-1)

yerr = np.std(off)/np.sqrt(b)

t_rel = np.arange(len(binned_on)) * sample_time * b
t = Time(t_rel + 1358385896 + b*sample_time/2, scale='utc', format='gps')
np.savetxt('ts_on_source_binned.txt', np.stack((t.mjd, binned_on)).T, fmt="%.13f %.18e")

plt.errorbar(t_rel, binned_on, yerr=yerr, label="On")
plt.errorbar(t_rel, binned_off, yerr=yerr, label="Off");

plt.xlabel("Time (s)")
plt.ylabel("Flux density (Jy/beam)")

plt.legend()
plt.tight_layout()

plt.savefig(f'binned_to_{b*sample_time:.0f}s.png')

# Save to pickle file
ds = np.empty((len(binned_on), 1, 4))
ds[:,0,0] = binned_on # Stokes I
ds[:,:,1:] = np.nan # Stokes Q, U, V
dat = {
    'DS': ds,
    'DS_MED': ds,
    'DS_STD': ds,
    'POLS': ['I', 'Q', 'U', 'V'],
    'TIMES': t.mjd*86400,
    'FREQS': np.array([161920000]),
    'TELESCOPE': 'MWA',
}
pkl = f'../dynspec/{t[0].gps:.0f}_ips.pkl'
with open(pkl, 'wb') as f:
    pickle.dump(dat, f)
