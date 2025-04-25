from astropy.table import QTable
import numpy as np
from scipy.optimize import curve_fit
from timing import *

table = QTable.read('toas.ecsv', format='ascii.ecsv')

def period_model(toas_mjd, period_d, pepoch_mjd):
    pulses, phases = np.divmod((toas_mjd - pepoch_mjd)/period_d, 1)
    rounded_pulses = np.round(pulses)
    predicted_toas = rounded_pulses * period_d + pepoch_mjd
    return predicted_toas

toas_mjd = table['ToA']
toas_mjd_err = table['ToA_err']

ephemeris = get_J1755_ephemeris()

p0 = [ephemeris['period'].to('d').value,
      ephemeris['PEPOCH'].mjd]

popt, pcov = curve_fit(period_model, toas_mjd, toas_mjd, p0=p0)
errs = np.sqrt(np.diag(pcov))
print("Period-only model:")
print(f'P            = {popt[0]*86400} ± {errs[0]*86400} s')
print(f'PEPOCH (MJD) = {popt[1]} ± {errs[1]}')
print()

def period_and_pdot_model(toas_mjd, period_d, pepoch_mjd, fdot):
    t = toas_mjd - pepoch_mjd
    pulses, phases = np.divmod(fdot*t**2 + t/period_d, 1)
    rounded_pulses = np.round(pulses)
    predicted_toas = rounded_pulses * period_d + pepoch_mjd
    return predicted_toas

p0.append(0.0)
popt, pcov = curve_fit(period_and_pdot_model, toas_mjd, toas_mjd, p0=p0)
errs = np.sqrt(np.diag(pcov))

print("Period and Pdot model:")
print(f'P            = {popt[0]*86400} ± {errs[0]*86400} s')
print(f'Fdot         = {popt[2]} ± {errs[2]} / day²')
print(f'PEPOCH (MJD) = {popt[1]} ± {errs[1]}')
print()

