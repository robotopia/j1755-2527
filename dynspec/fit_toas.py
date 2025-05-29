import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.table import QTable
import sys
import argparse

from timing import *

ephemeris = get_J1755_ephemeris()
src_coord = ephemeris['coord']
src_period = ephemeris['period']
#src_Porb = ephemeris['Porb']
src_pepoch = ephemeris['PEPOCH']
src_dm = ephemeris['DM']
src_tau_sc_1GHz = ephemeris['tau_sc']

def fit_toa(dat, fit_scattering=False, output_plot=None, obsid=None):

    t = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=EarthLocation.of_site(dat['TELESCOPE']))
    dt = t[1] - t[0]
    f = dat['FREQS'] * u.Hz
    try:
        df = f[1] - f[0]
    except:
        df = dat['BW'] * u.Hz
    nchans = len(f)

    f += 0.5*df # A hack because frequencies are labelled as lower edge of channels instead of middle of channels where they should be
    try:
        I, fmask, tmask = Stokes_ds(dat, pol='I')
    except:
        print(f"Couldn't open dynamic spectrum for {obsid}. Perhaps 'DS' wasn't correctly made?")
        return

    # Before dedispersing and dealing with nans, get dynamic spectrum and lightcurve stats
    noise = np.nanmedian(np.nanstd(I, axis=1)/np.sqrt(I.shape[1]))

    # Dedispersing requires no nans
    I[np.isnan(I)] = 0.0

    # Dedisperse
    f_ref = np.nanmean(f)
    Idd = dedisperse_ds(I, src_dm, f, f_ref, dt)

    # Reapply channel flags 
    if len(fmask) > 0:
        Idd[:,fmask] = np.nan

    # Form lightcurve
    lightcurve = np.nanmean(Idd, axis=1)

    # Debase: fit a baseline and remove it
    t_base = (t - t[0]).to('s')
    A = np.vstack([t_base.value, np.ones_like(t_base.value)]).T
    m, b = np.linalg.lstsq(A, lightcurve, rcond=None)[0]
    debased_lightcurve = lightcurve - (m*t_base.value + b)

    # Find the thing that looks most like our (unscattered) model pulse.
    # The assumption here is that the number of bins is enough to capture
    # the pulse. Some other clever strategy might be needed for different
    # resolutions / data sizes.
    model_pulse = pulse_model(np.arange(11)*dt.to('s').value, 1, 5*dt.to('s').value, 0.0, 0.0)
    matched_filter_results = np.convolve(debased_lightcurve, model_pulse, mode='same')
    peak_idx = np.argmax(matched_filter_results)
    A0 = matched_filter_results[peak_idx] / np.sum(model_pulse)

    # Fit a pulse to it
    p0 = (100*A0, t_base[peak_idx].value, 0.0, 0.0, 30) # Not 100% sure why we need 100*A0, but it helps!
    bounds = [(0,      -np.inf, -np.inf, -np.inf, 0.1*dt.to('s').value),
              (np.inf,  np.inf,  np.inf,  np.inf, np.inf              )]

    if fit_scattering:
        def scattered_pulse_model_at_freq(t_s, A, μ_s, m, c, σ_s):
            τ_s = src_tau_sc_1GHz.to('s').value * (f_ref/u.GHz).decompose()**(-4)
            return scattered_pulse_model(t_s, A, μ_s, m, c, σ_s, τ_s)
        model = scattered_pulse_model_at_freq
    else:
        model = pulse_model
    #print(f"{obsid = }")
    #print(f"{peak_idx = }")
    #print(f"{p0 = }")
    #print(f"{bounds = }")
    popt, pcov = curve_fit(model, t_base.to('s').value, lightcurve, p0=p0, sigma=noise if noise != 0.0 else None, bounds=bounds)
    A, μ, m, c, σ = popt
    A_err, μ_err, m_err, c_err, σ_err = np.sqrt(np.diag(pcov))
    Aσ_err = pcov[0,4] # Needed for peak flux error calculation

    # Pull out the ToA
    toa = t[0] + μ*u.s
    toa_err = (μ_err*u.s).to('d')

    # For plotting, do a bit of frequency scrunching
    fscrunch_factor = 1
    fscr = np.mean(np.reshape(f, (-1, fscrunch_factor)), axis=-1)
    Iddscr = np.mean(np.reshape(Idd, (-1, nchans//fscrunch_factor, fscrunch_factor)), axis=-1)

    #... and use a finer time resolution for the fitted pulse
    t_fine_base = np.arange(t_base[0].to('s').value, t_base[-1].to('s').value, dt.to('s').value/100) * u.s
    t_fine = t[0] + t_fine_base

    # The A parameter *is* the pulse fluence
    fluence = A * u.Jy * u.s
    fluence_err = A_err * u.Jy * u.s

    # Estimate peak S/N by getting std of lightcurve with fitted pulse subtracted
    noise = np.std(lightcurve - model(t_base.to('s').value, *popt)) * u.Jy
    peak_signal = A/σ/np.sqrt(2*np.pi) * u.Jy
    peak_signal_err = peak_signal * np.sqrt((A_err/A)**2 + (σ_err/σ)**2 + 2*Aσ_err/(A*σ))
    peak_snr = (peak_signal / noise).decompose().value
    #print(f'{peak_signal = }')
    #print(f'cf. = {A/σ/np.sqrt(2*np.pi)}')

    if output_plot is not None:

        # Make plot
        fig = plt.figure(figsize=(6, 10))
        gs = gridspec.GridSpec(5, 1, figure=fig)
        ax0 = fig.add_subplot(gs[:2, 0])
        ax1 = fig.add_subplot(gs[2:, 0], sharex=ax0)
        axs = [ax0, ax1]

        baseline = m*(t - t[0]).to('s').value + c
        axs[0].plot((t - t[0]).to('s'),
                    lightcurve - baseline,
                    label='data')
        #axs[0].plot((t_fine - t[0]).to('s'),
        #            model(t_fine_base.value, *p0),
        #            label='initial guess')
        axs[0].plot((t_fine - t[0]).to('s'),
                    model(t_fine_base.value, A, μ, m=0, c=0, σ_s=σ),
                    'r', alpha=0.4, label='model')
        if fit_scattering:
            modelled_pulse = pulse_model(t_fine_base.value, A, μ, m=0, c=0, σ_s=σ)
            axs[0].plot((t_fine - t[0]).to('s'), modelled_pulse,
                        'g--', alpha=0.4, label='model (no scattering)')

        axs[1].pcolormesh((t - t[0]).to('s').value, fscr.to('MHz').value, Iddscr.T)

        axs[0].set_ylabel("Flux density (Jy/beam)")
        axs[1].set_ylabel("Frequency (MHz)")
        axs[1].set_xlabel(f"Time (s) since MJD {t[0].mjd:.5f}")

        axs[0].legend()

        axs[0].set_title(f"{obsid or ''}\nToA = {toa.mjd:.5f} ± {toa_err:.5f}")

        plt.tight_layout()
        plt.savefig(output_plot)
        plt.close(fig)

    return f_ref, toa, toa_err, fluence, fluence_err, peak_signal, peak_signal_err, peak_snr, σ*u.s, σ_err*u.s


def main():
    parser = argparse.ArgumentParser(description="Fit ToAs to pulses")

    parser.add_argument('ds_files', nargs='*', help="The names of the input .pkl files")
    parser.add_argument('--fit_scattering', action='store_true', help="Include scattering in the fit")
    parser.add_argument('--output_table', default='toas.ecsv', help="Output filename for Astropy table of ToAs (default=toas.ecsv)")

    args = parser.parse_args()

    # Summarise everything into a table
    table = QTable(names=['ObsID', 'freq', 'ToA', 'ToA_err', 'telescope',
                          'fluence', 'fluence_err', 'fitted_peak_flux_density', 'fitted_peak_flux_density_err',
                          'peak_snr', 'width', 'width_err'],
                   dtype=[str, float, float, float, str, float, float, float, float, float, float, float],
                   units=[None, u.MHz, None, u.d, None, u.Jy*u.s, u.Jy*u.s, u.Jy, u.Jy, None, u.s, u.s])

    for pkl in args.ds_files:
        print(f"{pkl}...")
        obsid = pkl[:-4]

        # Load data
        dat = np.load(pkl, allow_pickle=True)

        output_plot = f"{pkl[:-4]}_toa.png"

        f_ref, toa, toa_err, fluence, fluence_err, peak_signal, peak_signal_err, peak_snr, σ, σ_err = fit_toa(dat, fit_scattering=args.fit_scattering, output_plot=output_plot, obsid=obsid)

        table.add_row([
            obsid,
            f_ref,
            toa.mjd,
            toa_err,
            dat['TELESCOPE'] if dat['TELESCOPE'] != "Siding Spring Observatory" else "ATCA",
            fluence,
            fluence_err,
            peak_signal,
            peak_signal_err,
            peak_snr,
            σ,
            σ_err,
        ])

    table.write(args.output_table, format="ascii.ecsv", overwrite=True)

if __name__ == '__main__':
    main()
