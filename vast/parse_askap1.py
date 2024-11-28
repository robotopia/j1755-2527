import numpy as np
import matplotlib.pyplot as plt

pulse_diff = np.loadtxt('askap1-mpath.txt')
pulse = np.cumsum(pulse_diff)

np.savetxt('askap1-I.txt', pulse)

plt.plot(pulse)
plt.savefig('askap1-I.png')
