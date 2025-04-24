import numpy as np
from scipy.optimize import curve_fit
import sys
import os
sys.path.append(os.path.abspath(".."))
from timing import *
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter('always')

dat = np.load('../1413381294_meerkat.pkl', allow_pickle=True)

t = dat['TIMES'] - dat['TIMES'][0]  # in seconds from start of obs
f = dat['FREQS'] / 1e6              # in MHz
I, _, _ = Stokes_ds(dat)

# Do some f-scrunching
factor = 11
f = np.nanmean(np.reshape(f, (-1, factor)), axis=1)
I = np.nanmean(np.reshape(I, (I.shape[0], -1, factor)), axis=2)

F, T = np.meshgrid(f, t)

Tflat = T.ravel()
Fflat = F.ravel()
Iflat = I.ravel()

mask = ~np.isnan(Iflat)

Tclean = Tflat[mask]
Fclean = Fflat[mask]
Iclean = Iflat[mask]

pixels = np.vstack((Tclean, Fclean))
all_pixels = np.vstack((Tflat, Fflat))

def calc_dmdelay(DM, f, f_ref):
    return 4148.808 * DM * (1/f**2 - 1/f_ref**2)

def model(pixel, comp1_1GHz, comp1_mu, comp1_sigma, comp1_alpha, comp2_1GHz, comp2_mu, comp2_sigma, comp2_alpha, baseline_m, baseline_c_1GHz, baseline_alpha, DM):

    t, f = pixel
    t_dd = t - 4148.808 * DM / f**2

    baseline = (baseline_m*t_dd + baseline_c_1GHz)*(f/1e3)**baseline_alpha
    comp1    = comp1_1GHz * (f/1e3)**comp1_alpha * np.exp(-0.5*(t_dd - comp1_mu)**2 / comp1_sigma**2)
    comp2    = comp2_1GHz * (f/1e3)**comp2_alpha * np.exp(-0.5*(t_dd - comp2_mu)**2 / comp2_sigma**2)

    return baseline + comp1 + comp2

p0 = [0.001, 220, 14, -3.7, 0.020, 340, 10, -3.7, 0, -0.01, 0, 710]
bounds = [
    (0.0001, 150, 1, -np.inf, 0.001, 280, 1, -np.inf, -np.inf, -np.inf, -np.inf, 500),
    (0.02,   280, 15, np.inf, 0.5,   400, 18, np.inf,  np.inf,  np.inf,  np.inf, 1500),
]
popt, pcov = curve_fit(model, pixels, Iclean, p0=p0, bounds=bounds)
errs = np.sqrt(np.diag(pcov))

print("Component 1:")
print(f"  peak     = ({popt[0]:.3f} ± {errs[0]:.3f}) (f/1GHz)^({popt[3]:.2f} ± {errs[3]:.2f})")
print(f"  position = {popt[1]:.0f} ± {errs[1]:.0f}")
print(f"  width    = {popt[2]:.1f} ± {errs[2]:.1f}")
print("Component 2:")
print(f"  peak     = ({popt[4]:.3f} ± {errs[4]:.3f}) (f/1GHz)^({popt[7]:.2f} ± {errs[7]:.2f})")
print(f"  position = {popt[5]:.0f} ± {errs[5]:.0f}")
print(f"  width    = {popt[6]:.1f} ± {errs[6]:.1f}")
print("Baseline:")
print(f"  ({popt[8]:.3f} * t + {popt[9]:.3f})*(f/1e3)**{popt[10]:.3f}")
print("DM:")
print(f"  DM = {popt[11]:.0f} ± {errs[11]:.0f} pc/cm³")

Imodel = model(all_pixels, *popt).reshape(I.shape)

print(np.sum((Iclean - model(pixels, *p0))**2))

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12,10))
axs[0,0].plot(t, np.nanmean(I, axis=1))
axs[0,1].plot(t, np.nanmean(Imodel, axis=1))
axs[0,2].plot(t, np.nanmean(I - Imodel, axis=1))

ylim = axs[0,0].get_ylim()
axs[0,1].set_ylim(ylim)
axs[0,2].set_ylim(ylim)

pc = axs[1,0].pcolormesh(t, f, I.T)
vmin, vmax = pc.get_clim()
axs[1,1].pcolormesh(t, f, Imodel.T, vmin=vmin, vmax=vmax)
axs[1,2].pcolormesh(t, f, (I - Imodel).T, vmin=vmin, vmax=vmax)
      
axs[0,0].set_title("Data")
axs[0,1].set_title("Model")
axs[0,2].set_title("Residuals")
      
axs[0,0].set_ylabel("Frequency (MHz)")
axs[1,0].set_xlabel("Time (s)")
axs[1,1].set_xlabel("Time (s)")
axs[1,2].set_xlabel("Time (s)")

plt.tight_layout()
plt.savefig('fit_component_spectra.png')
