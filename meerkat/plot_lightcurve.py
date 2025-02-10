# Script for plotting the lightcurve obtained from correlator imaging the MeerKAT fields
# (credit Dougal, Natasha). The input data file (1685306788_sdp_l0_1024ch_J1755-2527.pkl)
# is too big to keep in this repo. As of this writing, I'm putting a copy on Acacia:
# mwasci:smcsweeney/j1755-2527/meerkat_2023/1685306788_sdp_l0_1024ch_J1755-2527.pkl

import numpy as np
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
import matplotlib.pyplot as plt

dat = np.load('1685306788_sdp_l0_1024ch_J1755-2527.pkl', allow_pickle=True)
I = np.real(0.5*(dat['DS'][:,:,0] + dat['DS'][:,:,3])) * u.Jy
t = Time(dat['TIMES']/86400, scale='utc', format='mjd')
f = dat['FREQS']

PEPOCH = Time(59965.03767627493, scale='utc', format='mjd')
P = 4186.32874813198 * u.s

coord = SkyCoord('17:55:34.87 -25:27:49.1', unit=(u.hourangle, u.deg), frame='fk4')

meerkat = EarthLocation.of_site("MeerKAT")
bary_predicted_pulses = np.arange(2639, 2645)*P + PEPOCH
topo_predicted_pulses = bary_predicted_pulses - bary_predicted_pulses.light_travel_time(coord, ephemeris='jpl', location=meerkat)

plt.figure(figsize=(15,4))
plt.vlines(topo_predicted_pulses.mjd, -2, 6, ls='dashed', colors=['r'], label="Ephemeris-predicted ToAs")
#plt.vlines(bary_predicted_pulses.mjd, 0, 0.01, ls='dashed', colors=['g'])
plt.plot(t.mjd, np.nanmean(I, axis=1).to('mJy'))
plt.xlabel("Time (MJD)")
plt.ylabel("Flux density (mJy/beam)")
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig("lightcurve.png")
