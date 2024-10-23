import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt('peak_flux.txt')
toa_MJD, freq_MHz, bw_MHz, peak_flux_Jy = data.T

def powerlaw(nu_GHz, S1GHz, alpha):
    return S1GHz*nu_GHz**alpha

p0 = [0.1, -3.1]
popt, pcov = curve_fit(powerlaw, freq_MHz/1e3, peak_flux_Jy, p0=p0)

print(f"{popt = }")
print(f"{pcov = }")

X = np.logspace(-1, 0, 10) # frequencies in GHz
plt.plot(freq_MHz/1e3, peak_flux_Jy, 'k.', label="ASKAP and MWA points")
plt.plot(X, powerlaw(X, *popt), '--', label="Least squares fit")
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()

