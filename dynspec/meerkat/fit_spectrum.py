#!/usr/bin/env python

import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mticker
from matplotlib.ticker import MultipleLocator
import warnings
#from vis_ds import *
import glob
from scipy.optimize import curve_fit
from astropy import units as u

from astropy.coordinates import Angle, SkyCoord
from astropy import units as au
from astropy import constants as const

# Make sure to do export PYTHONPATH=$HOME/Programs/GPMTransient/dynspec/
import dedisperse_dynspec as dd

ref_nu = 1.e9

def pl(nu, norm, alpha):
    spec_nu = nu / ref_nu
    return norm * spec_nu **alpha

val = np.load("./J1755-25_1729341386_scan18.pkl", allow_pickle=True)
mjds = Time((val['TIMES']*u.s).to(u.day), format='mjd')
freqs = val['FREQS']
It = 1000*np.real((val['DS'][:,:,0]+val['DS'][:,:,3]))

# Remove small negative background
m = np.nanmean(It[0:75,:], axis=0)
bkg = np.tile(m, (It.shape[0],1))
It -= bkg

# Clip out final part of array as it has RFI
It = It[0:250, :]
mjds = mjds[0:250]

d = dd.Dynspec(dynspec=It, sample_time=2, freqlo=565.25, bw=0.1328125, time_offset=1413381294.4718554, transpose=True)
d.freq_ref = 1000
DM = 1200
d.dedisperse(DM)
d.t = d.get_time_at_infinite_frequency()

# Transpose and average to produce a light curve

lc = np.nanmean(d.dynspec.T, axis=1)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(mjds.gps, lc)
ax.axvline(1413381619.5701714, label="predicted TOA", color='k', ls=":")
ax.axvline(1413381612.5588295, label="dedispersed TOA", color='k', ls="-")
ax.set_xlabel("Time / GPS")
ax.set_ylabel("Average light curve / Jy beam$^{-1}$")
ax.legend()
fig.savefig("dedispersed_light_curve.png", bbox_inches="tight")

# Use light curve to generate some weights
cutoff = 10 #mJy
weights_1D = np.copy(lc)
weights_1D[weights_1D < cutoff] = 0.0
weights_2D = np.tile(weights_1D, (It.shape[1],1)).T

# Check weights resemble what we want
fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(121)
ax1.imshow(It.T, aspect='auto', origin="lower")
ax2 = fig.add_subplot(122)
ax2.imshow(weights_2D.T, aspect='auto', origin="lower")
fig.savefig("testing_weights.png", bbox_inches="tight")

# Produce average spectrum using the weights
# https://stackoverflow.com/questions/21113384/python-numpy-weighted-average-with-nans
dsI = np.ma.MaskedArray(It, mask=np.isnan(It))
weighted_spectrum = np.ma.average(dsI, weights=weights_2D, axis=0)

# Weighted spectrum
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(freqs, weighted_spectrum)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("Frequency / Hz")
ax.set_ylabel("Flux density / Jy")
fig.savefig("weighted_spectrum.png", bbox_inches="tight")

errs = np.nanstd(dsI[0:75,:], axis=0)
# But we are not just using single pixels, we are performing a weighted sum over the light curve
# To first order, we are using N measurements, where N is the number of non-zero-weight pixels in the light curve
# Errors roughly decrease as sqrt(N) for independent samples
weights_1D_num = np.copy(weights_1D)
weights_1D_num[weights_1D_num > 0] = 1
N = np.sum(weights_1D_num)
print(N)
errs = errs/np.sqrt(N)

# Don't fit the bottom end of the band because it has bandpass issues
clip = 80 
fit_freqs = freqs[clip:]
fit_flux = weighted_spectrum.data[clip:]
fit_errs = errs.data[clip:]

# Set any negative or zero fluxes to NaN
fit_errs[fit_flux <= 0] = np.nan

# Remove NaNs from the data that goes into the fit because they make scipy unhappy
fit_freqs = fit_freqs[~np.isnan(fit_flux)]
fit_errs = fit_errs[~np.isnan(fit_flux)]
fit_flux = fit_flux[~np.isnan(fit_flux)]
fit_freqs = fit_freqs[~np.isnan(fit_errs)]
fit_flux = fit_flux[~np.isnan(fit_errs)]
fit_errs = fit_errs[~np.isnan(fit_errs)]

#print(np.unique(np.isnan(fit_freqs)))
#print(np.unique(np.isnan(fit_flux)))
#print(np.unique(np.isnan(fit_errs)))

#print(len(fit_freqs))
#print(len(fit_flux))
#print(len(fit_errs))

#print(fit_freqs)
#print(fit_flux)
#print(fit_errs)

fit_p0 = [100, -3]
fit_res = curve_fit(
    pl,
    fit_freqs,
    fit_flux,
    fit_p0,
    sigma=fit_errs,
    absolute_sigma=True
)

amp = fit_res[0][0]
alpha = fit_res[0][1]

plot_freqs = np.linspace(np.nanmin(freqs), np.nanmax(freqs), 1000)

# Fitted weighted spectrum
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.axvline(freqs[clip]/1.e9, color='blue', ls=":")
ax.errorbar(freqs/1e9, weighted_spectrum, yerr=errs, lw=0, elinewidth=0.1, ls=None, color='black', alpha=0.5)
ax.scatter(freqs/1e9, weighted_spectrum, label="data", marker='.', color='black', alpha=0.3, s=10)
ax.plot(plot_freqs/1.e9, pl(plot_freqs, amp, alpha), zorder=10, label=f"$\\alpha=${alpha:3.1f}", alpha=0.8, color='red')
#ax.scatter(1.4, hi_flux, marker='*', color='k', label=f"HI flux density = {1.e3*hi_flux:2.1f}mJy")
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("Frequency / GHz")
ax.set_ylabel("Flux density / mJy")
ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.legend()
fig.savefig("fitted_weighted_spectrum.png", bbox_inches="tight")

hi_flux = pl(1.4e9, amp, alpha)
s_flux = pl(3e9, amp, alpha)
print("Flux density at 1.4 GHz = {0:3.1f} mJy".format(hi_flux))
print("Flux density at 3 GHz = {0:3.1f} mJy".format(s_flux))

