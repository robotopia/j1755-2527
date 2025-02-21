import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys

def get_J1755_ephemeris():
    return {
        'PEPOCH': Time(59965.03819710826, scale='utc', format='mjd'),
        'period': 4186.32874813198 * u.s,
        'DM': 710.0 * u.pc / u.cm**3,
        'tau_sc': 0.0*u.s,
        'coord': SkyCoord("17:55:34.87 -25:27:49.1", unit=(u.hourangle, u.deg), frame='icrs'),
    }

def fold(times, period, pepoch):

    pulses_phases = ((times - pepoch) / period).decompose()
    pulses, phases = np.divmod(pulses_phases + 0.5, 1)
    phases -= 0.5

    return pulses, phases

def barycentre(times, src_coord):
    bc_correction = times.light_travel_time(src_coord, ephemeris='jpl')
    return times + bc_correction

def calc_dmdelay(DM, f, f_ref):
    D = 4148.808 * u.MHz**2 / u.pc * u.cm**3 * u.s
    return D * DM * (1/f**2 - 1/f_ref**2)

def StokesI_ds(dat):
    if dat['POLS'][0] == 'XX' and dat['POLS'][3] == 'YY':
        return 0.5*np.real(dat['DS'][:,:,0] + dat['DS'][:,:,3])
    else:
        raise NotImplementedError(f"I haven't been taught how to deal with POLS = {ds['POLS']} yet")

def dedisperse_ds(ds, dm, f, f_ref, bin_size):
    dmdelays = calc_dmdelay(dm, f, f_ref)

    # It doesn't matter if the edges are wrapped, so use simple FFT to do
    # the DM delay adjustment (essentially, a "roll" with fractional bins)
    FFT = np.fft.rfft(ds, axis=0)
    FFT_freqs = np.fft.rfftfreq(ds.shape[0], d=bin_size)

    # Make a mesh grid of delays and frequencies with the same size as the FFT matrix
    DMDELAYS, FFT_FREQS = np.meshgrid(dmdelays, FFT_freqs)

    # Apply the phase ramp
    FFT *= np.exp(2j*np.pi*(FFT_FREQS*DMDELAYS).decompose())

    # Invert the FFT
    return np.fft.irfft(FFT, n=ds.shape[0], axis=0)

