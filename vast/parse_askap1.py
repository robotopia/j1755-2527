import numpy as np
import matplotlib.pyplot as plt

pulse_diff = np.loadtxt('askap1-mpath.txt')
pulse = np.cumsum(pulse_diff) / 1e3  # Because the ULP site expects Jy, not mJy

np.savetxt('askap1-I.txt', pulse)

plt.plot(np.arange(len(pulse))*10, pulse)
plt.xlabel("Time from start of observations (s)")
plt.ylabel("Flux density (Jy)")
plt.savefig('askap1-I.png')
