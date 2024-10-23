import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
import glob

def gauss(t, A, mu, sigma):
    return A*np.exp(-0.5*((t - mu)/sigma)**2) - 0.1

lightcurve_files = glob.glob("1410777264-I_and_1410777560-I_lightcurve_dm-????.txt")

widths = []
widths_err = []
dms = []
for lightcurve_file in lightcurve_files:
    data = np.loadtxt(lightcurve_file)

    t = data[5:,0]
    lightcurve = data[5:,1]

    p0 = (1, np.median(t), 50.0)

    popt, pcov = curve_fit(gauss, t, lightcurve, p0=p0)

    dms.append(int(lightcurve_file[44:48]))
    widths.append(popt[2])
    widths_err.append(np.sqrt(pcov[2,2]))

    #plt.clf()
    #plt.plot(t, lightcurve, label="Data")
    #plt.plot(t, gauss(t, *popt), 'k--', alpha=0.5, label="Fit")
    #plt.plot(t, gauss(t, *p0), 'r--', alpha=0.5, label="Initial guess")
    #plt.legend()
    #plt.show()

plt.clf()
plt.errorbar(dms, widths, yerr=widths_err, fmt='k.')
plt.xlabel("DM (pc/cm^3)")
plt.ylabel("Fitted pulse width (s)")
plt.savefig("1410777264-I_and_1410777560-I_lightcurve_width_vs_dm.png")
