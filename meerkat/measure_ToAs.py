# Measure ToAs from the lightcurve obtained from correlator imaging the MeerKAT fields
# (credit Dougal, Natasha). The input data file (1685306788_sdp_l0_1024ch_J1755-2527.pkl)
# is too big to keep in this repo. As of this writing, I'm putting a copy on Acacia:
# mwasci:smcsweeney/j1755-2527/meerkat_2023/1685306788_sdp_l0_1024ch_J1755-2527.pkl

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps, colors, cm
from scipy.optimize import curve_fit
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys

sys.path.append('../dynspec')

from timing import *

# Model for MeerKAT pulses
def gaussian(x, A, μ, σ, baseline, baseline_slope):
    return A*np.exp(-0.5*(x-μ)**2/σ**2) + baseline_slope*(x-μ) + baseline

pkl = '1685306788_sdp_l0_1024ch_J1755-2527.pkl'
print(f"Opening {pkl}...")

dat = np.load(pkl, allow_pickle=True)
I = StokesI_ds(dat)
freqs = dat['FREQS'] * u.Hz
f_ctr = np.mean(freqs)
bw = freqs[-1] - freqs[0] + freqs[1] - freqs[0]

ephemeris = get_J1755_ephemeris()
P = ephemeris['period']
PEPOCH = ephemeris['PEPOCH']
coord = ephemeris['coord']
DM = ephemeris['DM']

dmdelay = calc_dmdelay(DM, f_ctr, np.inf*u.Hz)

meerkat = EarthLocation.of_site("MeerKAT")
bary_predicted_toas = np.arange(2639, 2643)*P + PEPOCH
topo_predicted_toas = bary_predicted_toas - bary_predicted_toas.light_travel_time(coord, ephemeris='jpl', location=meerkat) + dmdelay

# One of the predicted ToAs is off to the side of an observation,
# where there is no data (yes, the pulse is truncated). Just add a
# a small offset so that it initialises correctly. It doesn't matter
# in terms of the fitting itself, so long as it converges.
topo_predicted_toas[1] += 5*u.s

times = Time(dat['TIMES']/86400.0, scale='utc', format='mjd')

# For each ToA, just pull out the snippet of data within a few minutes either side of a predicted pulse
t_radius = 3 * u.min

for toa in topo_predicted_toas:
    # Produce a dedispersed lightcurve that we can use to fit
    mask = np.abs((times - toa).to('s').value) <= t_radius.to('s').value
    t = times[mask]
    t0 = (t - t[0]).to('s').value # Use "seconds since start" for actual fitting
    dt = t[1] - t[0]
    I_cutout = I[mask,:]
    I_dd = dedisperse_ds(I_cutout, DM, freqs, f_ctr, dt)
    lightcurve = np.nanmean(I_dd, axis=1)

    # Do the fitting
    p0 = (np.max(lightcurve) - np.min(lightcurve), # A
          (toa - t[0]).to('s').value, # μ (seconds since beginning)
          10, # σ (seconds)
          np.min(lightcurve),
          0.0)
    bounds = [(0.0, -np.inf, 0.0, -np.inf, -np.inf,),
              (np.inf, np.inf, np.inf, np.inf, np.inf)]

    popt, pcov = curve_fit(gaussian, t0, lightcurve, p0=p0, bounds=bounds)
    fitted_toa = popt[1]*u.s + t[0]
    err = pcov[1,1]*u.s

    fig, axs = plt.subplots(nrows=2, sharex=True)
    axs[0].plot(t.mjd, lightcurve)
    axs[0].plot(t.mjd, gaussian(t0, *popt), alpha=0.5)
    axs[1].pcolormesh(t.mjd, freqs.to('MHz').value, I_dd.T)
    axs[0].axvline(fitted_toa.mjd, c='r', ls='--')
    print(f'{fitted_toa.mjd = }\n{f_ctr.to("MHz") = }\n{bw.to("MHz") = }\n{err.to("d") = }\n')

    axs[0].set_title(f"Fitted ToA = {fitted_toa.mjd:.20f}")
    plt.savefig(f'{fitted_toa.mjd}.png')
    plt.close(fig)


