import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf, erfc, erfcx
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys

π = np.pi

def get_J1755_ephemeris():
    return {
        'PEPOCH': Time(59965.037829, scale='utc', format='mjd'),
        'period': 4186.328785 * u.s,
        'DM': 771.5 * u.pc / u.cm**3,
        'tau_sc': 7.0e-2*u.s,
        'coord': SkyCoord("17:55:34.87 -25:27:49.1", unit=(u.hourangle, u.deg), frame='icrs'),
    }

def gaussian_cdf(t, A, μ, σ):
    '''
    ... not *quite* CDF of a Gaussian, since this is scaled up by A.
    '''
    return A*0.5*(1 + erf((t - μ)/np.sqrt(2)/σ))

def discrete_gaussian(t, A, μ, σ):
    #return A*np.exp(-0.5*(x - μ)**2 / σ**2)
    # Integrating over a bin width lets us preserve fluence when fitting

    dt = t[1] - t[0]
    return (gaussian_cdf(t + dt/2, A, μ, σ) - gaussian_cdf(t - dt/2, A, μ, σ)) / dt

def emg_cdf(t, A, μ, σ, τ):
    '''
    ... not *quite* CDF of an exponentially modified Gaussian, since this is scaled up by A.
    '''
    λ = 1/τ
    arg1 = (μ + λ*σ**2 - t)/(np.sqrt(2)*σ)
    arg2 = λ/2 * (2*μ + λ*σ**2 - 2*t)

    # The following two are mathematically equivalent, but computationally different in
    # how they are subject to problems of overflow
    #return A*(gaussian_cdf(t, 1, μ, σ) - 0.5*np.exp(arg2) * erfc(arg1))
    return A*(gaussian_cdf(t, 1, μ, σ) - 0.5*np.exp(arg2 - arg1**2) * erfcx(arg1))

def discrete_emg(t, A, μ, σ, τ):
    '''
    Exponentially modified gaussian
    '''
    dt = t[1] - t[0]
    pdf_value = (emg_cdf(t + dt/2, A, μ, σ, τ) - emg_cdf(t - dt/2, A, μ, σ, τ)) / dt
    #print(f'{A = }, {μ = }, {σ = }, {τ = }')#, {pdf_value = }')
    return pdf_value


def pulse_model(t_s, A, μ_s, m=0, c=0, σ_s=30):
    baseline = m*t_s + c
    return discrete_gaussian(t_s, A, μ_s, σ_s) + baseline

def scattered_pulse_model(t_s, A, μ_s, m, c, σ_s, τ_s):
    baseline = m*t_s + c
    return discrete_emg(t_s, A, μ_s, σ_s, τ_s) + baseline

def pulse_fluence(A, μ_s, m=0, c=0, σ_s=30):
    return A*σ_s*np.sqrt(2*np.pi)

def fold(times, period, pepoch):

    pulses_phases = ((times - pepoch) / period).decompose()
    pulses, phases = np.divmod(pulses_phases + 0.5, 1)
    phases -= 0.5

    return pulses, phases

def barycentre(times, src_coord):
    bc_correction = times.light_travel_time(src_coord, ephemeris='jpl')
    return times + bc_correction

def topocentre(times, src_coord):
    bc_correction = times.light_travel_time(src_coord, ephemeris='jpl')
    return times - bc_correction

def calc_dmdelay(DM, f, f_ref):
    D = 4148.808 * u.MHz**2 / u.pc * u.cm**3 * u.s
    return D * DM * (1/f**2 - 1/f_ref**2)

def Stokes_ds(dat, pol='I', pb_corr=True):
    if pol not in ['I', 'Q', 'U', 'V']:
        raise NotImplementedError(f"Expecting pol to be one of [I, Q, U, V]")

    if pb_corr:

        if 'PB_CORR' not in dat.keys():
            print("WARNING: no 'PB_CORR' field found in pickle file; no PB correction applied")
            dat['PB_CORR'] = np.ones((len(dat['FREQS']), len(dat['POLS'])))
        else:
            for i in range(len(dat['POLS'])):
                this_pol = dat['POLS'][i]

            if np.all(np.isnan(dat['PB_CORR'][:,i])):
                print(f"WARNING: PB corrections for {this_pol} are all nans. No correction applied to this pol.")
                dat['PB_CORR'][:,i] = np.ones(dat['FREQS'].shape)

    else:
        # This creates a PB correction field (full of ones) in memory, but doesn't save to disk
        dat['PB_CORR'] = np.ones((len(dat['FREQS']), len(dat['POLS'])))

    if dat['POLS'] == ['XX', 'XY', 'YX', 'YY']:
        if pol == 'I':
            S = 0.5*np.real(dat['DS'][:,:,0]/dat['PB_CORR'][:,0] + dat['DS'][:,:,3]/dat['PB_CORR'][:,3])
        elif pol == 'Q':
            S = 0.5*np.real(dat['DS'][:,:,0]/dat['PB_CORR'][:,0] - dat['DS'][:,:,3]/dat['PB_CORR'][:,3])
        elif pol == 'U':
            S = 0.5*np.real(dat['DS'][:,:,1]/dat['PB_CORR'][:,1] + dat['DS'][:,:,2]/dat['PB_CORR'][:,2])
        elif pol == 'V':
            S = 0.5*np.imag(dat['DS'][:,:,1]/dat['PB_CORR'][:,1] - dat['DS'][:,:,2]/dat['PB_CORR'][:,2])
    elif dat['POLS'] == ['I', 'Q', 'U', 'V']:
        if pol == 'I':
            S = np.real(dat['DS'][:,:,0]/dat['PB_CORR'][:,0])
        elif pol == 'Q':
            S = np.real(dat['DS'][:,:,1]/dat['PB_CORR'][:,1])
        elif pol == 'U':
            S = np.real(dat['DS'][:,:,2]/dat['PB_CORR'][:,2])
        elif pol == 'V':
            S = np.real(dat['DS'][:,:,3]/dat['PB_CORR'][:,3])
    else:
        raise NotImplementedError(f"I haven't been taught how to deal with POLS = {ds['POLS']} yet")

    # Get and apply channel and time flags
    fmask, tmask = get_flags(dat)

    if len(fmask) > 0:
        S[:,fmask] = np.nan

    if len(tmask) > 0:
        S[tmask,:] = np.nan

    return S, fmask, tmask

def get_flags(dat):
    # Get channel flags
    if 'CHAN_FLAGS' in dat.keys():
        fmask = set(dat['CHAN_FLAGS'])
    else:
        fmask = set()

    if len(dat['FREQS']) > 1: # A hack to save my single-frequency dynspecs
        fmask |= set(np.where(np.all(np.isnan(dat['DS']), axis=0))[0])
    fmask = list(fmask)

    # Get time flags
    if 'TIME_FLAGS' in dat.keys():
        tmask = dat['TIME_FLAGS']
    else:
        tmask = []

    return fmask, tmask

def dedisperse_ds(ds, dm, f, f_ref, bin_size):
    dmdelays = calc_dmdelay(dm, f, f_ref)
    iscomplex = np.iscomplexobj(ds)

    # It doesn't matter if the edges are wrapped, so use simple FFT to do
    # the DM delay adjustment (essentially, a "roll" with fractional bins)
    FFT = np.fft.fft(ds, axis=0)
    FFT_freqs = np.fft.fftfreq(ds.shape[0], d=bin_size)

    # Make a mesh grid of delays and frequencies with the same size as the FFT matrix
    DMDELAYS, FFT_FREQS = np.meshgrid(dmdelays, FFT_freqs)

    # Apply the phase ramp
    FFT *= np.exp(2j*np.pi*FFT_FREQS*DMDELAYS)

    # Invert the FFT
    IFFT = np.fft.ifft(FFT, n=ds.shape[0], axis=0)

    # If the original dynamic spectrum was real-valued, just get back the real part
    if not iscomplex:
        IFFT = np.real(IFFT)

    return IFFT

