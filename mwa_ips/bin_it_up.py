import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

on = np.loadtxt('ts_on_source.txt')
off = np.loadtxt('ts_off_source.txt')

mean_on = np.mean(on)
mean_off = np.mean(off)

# From Angie:
# On source: pixel coordinates 1546, 1258, value of 1.07391
# off source: pixel coordinates 1549, 1563, value of 0.00331

mean_on_should_be = 1.07391
mean_off_should_be = 0.00331

mean_on_adjust = mean_on_should_be - mean_on
mean_off_adjust = mean_off_should_be - mean_off

on += mean_on_adjust
off += mean_off_adjust

b = 20 # Number of time samples to bin together
sample_time = 0.5 # seconds

binned_on = np.mean(np.reshape(on, (-1, b)), axis=-1)
binned_off = np.mean(np.reshape(off, (-1, b)), axis=-1)

yerr = np.std(off)/np.sqrt(b)

t_rel = np.arange(b) * sample_time * b
t = Time(t_rel + 1358386096 + b*sample_time/2, scale='utc', format='gps')
np.savetxt('ts_on_source_binned.txt', np.stack((t.mjd, binned_on)).T, fmt="%.13f %.18e")

plt.errorbar(t_rel, binned_on, yerr=yerr, label="On")
plt.errorbar(t_rel, binned_off, yerr=yerr, label="Off");

plt.xlabel("Time (s)")
plt.ylabel("Flux density (Jy/beam)")

plt.legend()
plt.tight_layout()

plt.savefig(f'binned_to_{b*sample_time:.0f}s.png')
