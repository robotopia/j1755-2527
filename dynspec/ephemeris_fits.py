from astropy.table import QTable
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import numpy as np
from scipy.optimize import curve_fit
from timing import *

table = QTable.read('toas.ecsv', format='ascii.ecsv')
ephemeris = get_J1755_ephemeris()

def period_model(toas_mjd, period_d, pepoch_mjd):
    _, phases = np.divmod((toas_mjd - pepoch_mjd)/period_d + 0.5, 1)
    return phases - 0.5

toas_mjd = []
for row in table:

    location = EarthLocation.of_site(row['telescope'])
    freq = row['freq']
    toa = Time(row['ToA'], scale='utc', format='mjd', location=location)
    toa_bary_dd = barycentre(toa, ephemeris['coord']) - calc_dmdelay(ephemeris['DM'], freq, np.inf*u.MHz)
    toas_mjd.append(toa_bary_dd.mjd)

toas_mjd = np.array(toas_mjd)
toas_mjd_err = table['ToA_err']

p0 = [ephemeris['period'].to('d').value,
      ephemeris['PEPOCH'].mjd]

popt, pcov = curve_fit(period_model, toas_mjd, np.zeros(toas_mjd.shape), p0=p0)
errs = np.sqrt(np.diag(pcov))
print("Period-only model:")
print(f'P            = {popt[0]*86400:.5f} ± {errs[0]*86400:.5f} s')
print(f'PEPOCH (MJD) = {popt[1]:.6f} ± {errs[1]:.6f}')
print()

def period_and_pdot_model(toas_mjd, period_d, pepoch_mjd, fdot):
    t = toas_mjd - pepoch_mjd
    _, phases = np.divmod(fdot*t**2 + t/period_d + 0.5, 1)
    return phases - 0.5

p0.append(0.0)
popt, pcov = curve_fit(period_and_pdot_model, toas_mjd, np.zeros(toas_mjd.shape), p0=p0)
errs = np.sqrt(np.diag(pcov))

print("Period and Pdot model:")
print(f'P            = {popt[0]*86400} ± {errs[0]*86400} s')
print(f'Fdot         = {popt[2]} ± {errs[2]} / day²')
print(f'PEPOCH (MJD) = {popt[1]} ± {errs[1]}')
print()

