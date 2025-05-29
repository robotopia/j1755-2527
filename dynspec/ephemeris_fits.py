from astropy.table import QTable
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import numpy as np
from scipy.optimize import curve_fit
from timing import *

table = QTable.read('toas.ecsv', format='ascii.ecsv')
ephemeris = get_J1755_ephemeris()
freqs = table['freq']

def period_model(toas_mjd, period_d, pepoch_mjd):
    _, phases = np.divmod((toas_mjd - pepoch_mjd)/period_d + 0.5, 1)
    return phases - 0.5

toas_mjd = []
for row in table:

    telescope = row['telescope'] if row['telescope'] != 'ATCA' else 'Siding Spring Observatory'
    location = EarthLocation.of_site(telescope)
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
print(f'  P            = {popt[0]*86400:.5f} ± {errs[0]*86400:.5f} s')
print(f'  PEPOCH (MJD) = {popt[1]:.6f} ± {errs[1]:.6f}')
print()

def period_and_pdot_model(toas_mjd, period_d, pepoch_mjd, fdot):
    t = toas_mjd - pepoch_mjd - calc_dmdelay(ephemeris['DM'], freqs, np.inf*u.MHz).to('d').value
    _, phases = np.divmod(fdot*t**2 + t/period_d + 0.5, 1)
    return phases - 0.5

p0.append(0.0)
popt, pcov = curve_fit(period_and_pdot_model, toas_mjd, np.zeros(toas_mjd.shape), p0=p0)
errs = np.sqrt(np.diag(pcov))

P = popt[0] * u.d
P_err = errs[0] * u.d
F = 1/P
Fdot = popt[2] / u.d**2
Fdot_err = errs[2] / u.d**2
Pdot = -Fdot/F**2
Pdot_err = np.abs(Pdot * np.sqrt((Fdot_err/Fdot)**2 + 2*(P_err/P)**2))

print(f"Period and Pdot model (assumes DM of {ephemeris['DM']}):")
print(f'  P            = {popt[0]*86400:.5f} ± {errs[0]*86400:.5f} s')
print(f'  Fdot         = {Fdot.to("Hz2").value} ± {Fdot_err.to("Hz2").value} Hz²')
print(f'  Pdot         = {Pdot.decompose().value} ± {Pdot_err.decompose().value}')
print(f'  PEPOCH (MJD) = {popt[1]:.6f} ± {errs[1]:.6f}')
print()

def period_and_dm_model(toas_mjd_and_freqs, period_d, pepoch_mjd, dm):
    toas_mjd, freqs = toas_mjd_and_freqs.T
    t = toas_mjd - pepoch_mjd - calc_dmdelay(dm*u.pc/u.cm**3 - ephemeris['DM'], freqs*u.MHz, np.inf*u.MHz).to('d').value
    _, phases = np.divmod(t/period_d + 0.5, 1)
    return phases - 0.5

p0[2] = 710
toas_mjd_and_freqs = np.array([toas_mjd, freqs.to('MHz').value]).T
popt, pcov = curve_fit(period_and_dm_model, toas_mjd_and_freqs, np.zeros(toas_mjd.shape), p0=p0)
errs = np.sqrt(np.diag(pcov))

print("Period and DM model:")
print(f'  P            = {popt[0]*86400:.5f} ± {errs[0]*86400:.5f} s')
print(f'  PEPOCH (MJD) = {popt[1]:.6f} ± {errs[1]:.6f}')
print(f'  DM           = {popt[2]:.0f} ± {errs[2]:.0f} pc/cm³')
print()

