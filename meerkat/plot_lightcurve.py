# Script for plotting the lightcurve obtained from correlator imaging the MeerKAT fields
# (credit Dougal, Natasha). The input data file (1685306788_sdp_l0_1024ch_J1755-2527.pkl)
# is too big to keep in this repo. As of this writing, I'm putting a copy on Acacia:
# mwasci:smcsweeney/j1755-2527/meerkat_2023/1685306788_sdp_l0_1024ch_J1755-2527.pkl

import numpy as np
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.constants import c
import matplotlib.pyplot as plt

dat = np.load('1685306788_sdp_l0_1024ch_J1755-2527.pkl', allow_pickle=True)
I = np.real(0.5*(dat['DS'][:,:,0] + dat['DS'][:,:,3])) * u.Jy
Q = np.real(0.5*(dat['DS'][:,:,0] - dat['DS'][:,:,3])) * u.Jy
U = np.real(0.5*(dat['DS'][:,:,1] + dat['DS'][:,:,2])) * u.Jy
V = np.imag(0.5*(dat['DS'][:,:,1] - dat['DS'][:,:,2])) * u.Jy
t = Time(dat['TIMES']/86400, scale='utc', format='mjd')

# Defaraday rotate
f = dat['FREQS'] * u.Hz
λ = c/f
L = Q + U*1j
RM = 961 * u.rad / u.m**2
φ = (RM * λ**2).to('rad')
Lfr = L * np.exp(1j*φ[np.newaxis,:] / u.rad)
Lfr2 = L * np.exp(-1j*φ[np.newaxis,:] / u.rad)
plt.plot(f.to("GHz"), np.mod(φ.value, 2*np.pi) - np.pi)
plt.xlabel("Frequency (GHz)")
plt.ylabel("RM phase (rad)")
plt.savefig("rmphase.png")

PEPOCH = Time(59965.03767627493, scale='utc', format='mjd')
P = 4186.32874813198 * u.s

coord = SkyCoord('17:55:34.87 -25:27:49.1', unit=(u.hourangle, u.deg), frame='fk4')

meerkat = EarthLocation.of_site("MeerKAT")
bary_predicted_pulses = np.arange(2639, 2645)*P + PEPOCH
topo_predicted_pulses = bary_predicted_pulses - bary_predicted_pulses.light_travel_time(coord, ephemeris='jpl', location=meerkat)

fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(15,12))
for ax in axs:
    ax.vlines(topo_predicted_pulses.mjd, -2, 6, ls='dashed', colors=['r'], label="Ephemeris-predicted ToAs")
    ax.set_ylabel("Flux density (mJy/beam)")
#plt.vlines(bary_predicted_pulses.mjd, 0, 0.01, ls='dashed', colors=['g'])
axs[0].plot(t.mjd, np.nanmean(I, axis=1).to('mJy'), label="Total intensity")
axs[1].plot(t.mjd, np.abs(np.nanmean(Lfr, axis=1).to('mJy')), label="Linear (RM = ∓961 rad/m²)")
axs[1].plot(t.mjd, np.abs(np.nanmean(Lfr2, axis=1).to('mJy')), label="Linear (RM = ±961 rad/m²)")
axs[2].plot(t.mjd, np.nanmean(V, axis=1).to('mJy'), label="Circular")
axs[-1].set_xlabel("Time (MJD)")
for ax in axs:
    ax.legend(loc="upper right")
plt.tight_layout()
#plt.show()
plt.savefig("lightcurve.png")
