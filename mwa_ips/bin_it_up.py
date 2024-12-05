import numpy as np
import matplotlib.pyplot as plt

on = np.loadtxt('ts_on_source.txt')
off = np.loadtxt('ts_off_source.txt')

b = 20 # Number of time samples to bin together
sample_time = 0.5 # seconds

binned_on = np.mean(np.reshape(on, (-1, b)), axis=-1)
binned_off = np.mean(np.reshape(off, (-1, b)), axis=-1)

yerr = np.std(off)/np.sqrt(b)

plt.errorbar(np.arange(b), binned_on, yerr=yerr, label="On")
plt.errorbar(np.arange(b), binned_off, yerr=yerr, label="Off");

plt.xlabel("Time bin")
plt.ylabel("Flux density (Jy/beam)")

plt.legend()
plt.tight_layout()

plt.savefig(f'binned_to_{b*sample_time:.0f}s.png')
