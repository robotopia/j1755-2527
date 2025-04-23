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

def fit_toa(dat, fit_scattering=False, output_plot=None, obsid=None):

    t = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=EarthLocation.of_site(dat['TELESCOPE']))
    dt = t[1] - t[0]
    f = dat['FREQS'] * u.Hz
    df = f[1] - f[0]
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
    if fit_scattering:
        p0 = (A0, t_base[peak_idx].value, 0.1*dt.to('s').value, 0.0, 0.0, 1.5)
        bounds = [(0, -np.inf, 1e-6*dt.to('s').value, -np.inf, -np.inf, 0.1*dt.to('s').value), (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)]
        model = scattered_pulse_model
    else:
        p0 = (A0, t_base[peak_idx].value, 0.0, 0.0)
        bounds = [(0, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)]
        model = pulse_model
    #print(f"{obsid = }")
    #print(f"{peak_idx = }")
    #print(f"{p0 = }")
    #print(f"{bounds = }")
    popt, pcov = curve_fit(model, t_base.to('s').value, lightcurve, p0=p0, sigma=noise, bounds=bounds)

    # Pull out the ToA
    toa = t[0] + popt[1]*u.s
    toa_err = (np.sqrt(pcov[1,1])*u.s).to('d')

    # Estimate peak S/N by getting std of lightcurve with fitted pulse subtracted
    noise = np.std(lightcurve - model(t_base.to('s').value, *popt)) * u.Jy
    peak_signal = popt[0] * u.Jy # *** MAY NEED TO BE CHANGED IF PULSE MODEL CHANGES ***
    peak_snr = (peak_signal / noise).decompose().value

    # For plotting, do a bit of frequency scrunching
    fscrunch_factor = 1
    fscr = np.mean(np.reshape(f, (-1, fscrunch_factor)), axis=-1)
    Iddscr = np.mean(np.reshape(Idd, (-1, nchans//fscrunch_factor, fscrunch_factor)), axis=-1)

    #... and use a finer time resolution for the fitted pulse
    t_fine_base = np.arange(t_base[0].to('s').value, t_base[-1].to('s').value, dt.to('s').value/100) * u.s
    t_fine = t[0] + t_fine_base

    # Calculate pulse fluence from model
    if fit_scattering:
        A, μ, _, m, c, σ = popt
        fluence = pulse_fluence(A, μ, σ) * u.Jy * u.s
    else:
        fluence = pulse_fluence(*popt) * u.Jy * u.s
    fluence_err = fluence*np.sqrt(pcov[0,0])/popt[0]  # <--- If pulse model changes, change this accordingly

    if output_plot is not None:
        if fit_scattering:

            # Make plot
            fig = plt.figure(figsize=(6, 10))
            gs = gridspec.GridSpec(5, 1, figure=fig)
            ax0 = fig.add_subplot(gs[0, 0])
            ax1 = fig.add_subplot(gs[1:3, 0], sharex=ax0)
            ax2 = fig.add_subplot(gs[3:, 0])
            axs = [ax0, ax1, ax2]

            axs[0].plot((t - t[0]).to('s'), lightcurve, label='data')
            A, μ, τ, m, c, σ = popt
            axs[0].plot((t_fine - t[0]).to('s'), pulse_model(t_fine_base.value, A, μ, m=m, c=c, σ_s=σ), 'g--', alpha=0.4, label='model (no scattering)')

            # Plot error ellipse
            mini_pcov = np.array([[pcov[5,5], pcov[5,2]], [pcov[2,5], pcov[2,2]]])
            eigenvalues, eigenvectors = np.linalg.eigh(mini_pcov)
            mean = [σ, τ]
            axs[2].scatter(*mean, color='red', marker='x')
            for sigma in [1, 2, 3]:
                width, height = sigma * np.sqrt(eigenvalues)
                angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
                ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle, edgecolor='black', facecolor='none', lw=2, alpha=1-0.25*sigma)
                axs[2].add_patch(ellipse)#, label=f'{sigma}σ')
            axs[2].set_xlabel("σ of unscattered pulse (s)")
            axs[2].set_ylabel(f"Scattering timescale at {f_ref.to('MHz')} (s)")
            #axs[2].set_xlim(mean[0] - width if mean[0] - width > 0 else 0, mean[0] + width)
            #axs[2].set_ylim(mean[1] - height if mean[1] - height > 0 else 0, mean[1] + height)
            axs[2].set_xlim([0.0, 4.0])
            axs[2].set_ylim([0.0, 4.0])
            axs[2].axhline(2.7477e-4*((f_ref/u.GHz).decompose())**-4, ls='--', c='g', label='NE2001')
            axs[2].axhline(3.1225e-5*((f_ref/u.GHz).decompose())**-4, ls='--', c='orange', label='YMW16')
            axs[2].legend()
        else:
            fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 9))
            axs[0].plot((t - t[0]).to('s'), lightcurve, label='data')

        axs[0].plot((t_fine - t[0]).to('s'), model(t_fine_base.value, *popt), 'r', alpha=0.4, label='model')
        axs[1].pcolormesh((t - t[0]).to('s').value, fscr.to('MHz').value, Iddscr.T)

        axs[0].set_ylabel("Flux density (Jy/beam)")
        axs[1].set_ylabel("Frequency (MHz)")
        axs[1].set_xlabel("Time (MJD)")

        axs[0].legend()

        axs[0].set_title(f"{obsid or ''}\nToA = {toa.mjd} ± {toa_err}")

        plt.tight_layout()
        plt.savefig(output_plot)
        plt.close(fig)

    return f_ref, toa, toa_err, fluence, fluence_err, peak_signal, peak_snr


def main():
    parser = argparse.ArgumentParser(description="Fit ToAs to pulses")

    parser.add_argument('ds_files', nargs='*', help="The names of the input .pkl files")
    parser.add_argument('--fit_scattering', action='store_true', help="Include scattering in the fit")
    parser.add_argument('--output_table', default='toas.ecsv', help="Output filename for Astropy table of ToAs (default=toas.ecsv)")

    args = parser.parse_args()

    # Summarise everything into a table
    table = QTable(names=['ObsID', 'freq', 'ToA', 'ToA_err', 'telescope', 'fluence', 'fluence_err', 'fitted_peak_flux_density', 'peak_snr'],
                   dtype=[str, float, float, float, str, float, float, float, float],
                   units=[None, u.MHz, None, u.d, None, u.Jy*u.s, u.Jy*u.s, u.Jy, None])

    for pkl in args.ds_files:
        print(f"{pkl}...")
        obsid = pkl[:-4]

        # Load data
        dat = np.load(pkl, allow_pickle=True)

        output_plot = f"{pkl[:-4]}_toa.png"

        f_ref, toa, toa_err, fluence, fluence_err, peak_signal, peak_snr = fit_toa(dat, fit_scattering=args.fit_scattering, output_plot=output_plot, obsid=obsid)

        table.add_row([obsid, f_ref, toa.mjd, toa_err, dat['TELESCOPE'], fluence, fluence_err, peak_signal, peak_snr])

    table.write(args.output_table, format="ascii.ecsv", overwrite=True)

if __name__ == '__main__':
    main()
