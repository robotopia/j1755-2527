import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.constants import c
import sys
import argparse

from timing import *

import warnings

def safe_nanmean(arr, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        np.seterr(all="ignore")
        result = np.nanmean(arr, **kwargs)
        np.seterr(all="warn")
    return result

# Model
def L_model(f_MHz, RM, L_1GHz, α, PA_rad):
    # RM in rad/m^2
    # L_1GHz is linearly polarised flux density at 1 GHz
    # α is spectral index of linearly polarised component
    # PA_rad is the intrinsic polarisation angle, in radians
    L_f = L_1GHz*(f_MHz/1000)**α # Abs value of linear polarisation at given frequencies
    λ = 299792458/(f_MHz*1e6) # Wavelengths in m
    ψ = RM * λ**2 + PA_rad
    return L_f*np.exp(2j*ψ)

# But curve_fit doesn't deal with models with complex outputs, so here's a wrapper
def real_L_model(f_MHz, RM, L_1GHz, α, PA_rad):
    L = L_model(f_MHz, RM, L_1GHz, α, PA_rad)
    return np.concatenate([L.real, L.imag])

def parse_args():
    parser = argparse.ArgumentParser(description="Measure RM")

    parser.add_argument('--ds_files', nargs='*', help="The name of the input .pkl file")
    parser.add_argument('--output_image', help="The filename to use for the average DS image. If not given, will plt.show().")
    parser.add_argument('--fscrunch_factor', type=int, help="Average this many channels together for the output plot")
    parser.add_argument('--bins', type=int, nargs='*', help="Use the specified bins for the RM fit")

    return parser.parse_args()

def get_Q_and_U(dat, bins):
    Q = safe_nanmean(Stokes_ds(dat, pol='Q')[0][bins,:], axis=0)
    U = safe_nanmean(Stokes_ds(dat, pol='U')[0][bins,:], axis=0)

    return Q, U, f

def fit_RM_to_QU(Q, U, f, p0=None, bounds=None):
    mask = ~np.logical_or(np.isnan(Q), np.isnan(U))

    xdata = f[mask].to('MHz').value
    ydata = np.concatenate([Q[mask], U[mask]])

    return curve_fit(real_L_model, xdata, ydata, p0=p0, bounds=bounds)

def plot_RM(ax, f, Q, U, popt=None, pcov=None, fscrunch_factor=None):

    if fscrunch_factor is not None:
        favg = safe_nanmean(np.reshape(f, (-1, fscrunch_factor)), axis=-1)
        Qavg = safe_nanmean(np.reshape(Q, (-1, fscrunch_factor)), axis=-1)
        Uavg = safe_nanmean(np.reshape(U, (-1, fscrunch_factor)), axis=-1)
    else:
        favg, Qavg, Uavg = f, Q, U

    ax.plot(favg.to('MHz'), Qavg, label='Q data', alpha=0.5)
    ax.plot(favg.to('MHz'), Uavg, label='U data', alpha=0.5)

    if popt is not None and pcov is not None:
        # Unpack popt and pcov
        RM, L_1GHz, α, PA_rad = popt
        RM_err, L_1GHz_err, α_err, PA_rad_err = np.sqrt(np.diag(pcov))
        PA_rad = (np.divmod(PA_rad / np.pi + 0.5, 1)[-1] - 0.5) * np.pi

        f_fine = np.linspace(f[0].value, f[-1].value, 1000) * f.unit
        L = L_model(f_fine.to('MHz').value, *popt)

        ax.plot(f_fine.to('MHz'), np.real(L), label='Q fit')
        ax.plot(f_fine.to('MHz'), np.imag(L), label='U fit')
        ax.plot(f_fine.to('MHz'), np.abs(L), label='|L| fit', alpha=0.5, ls='--')

    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Flux density (Jy)")

    ax.set_title(f"L = $({L_1GHz*1e3:.1f} \\pm {L_1GHz_err*1e3:.1f})$ mJy ${{}} \\times \\left(\\frac{{\\nu}}{{1 \\, GHz}}\\right)^{{{α:.1f} \\pm {α_err:.1f}}}$\nRM ${{}} = {RM:.1f} \\pm {RM_err:.1f}$ rad/m²\nPA ${{}} = ({np.rad2deg(PA_rad):.0f} \\pm {np.rad2deg(PA_rad_err):.0f})^\\circ$")

    plt.legend()


def main():
    args = parse_args()
    fig = plt.figure(figsize=(4*len(args.ds_files),5))
    gs = gridspec.GridSpec(4, len(args.ds_files), figure=fig, hspace=0.0)

    for i in range(len(args.ds_files)):
        ds_file = args.ds_files[i]
        dat = np.load(ds_file, allow_pickle=True)

        ax_RM = fig.add_subplot(gs[0, i])
        ax_PA = fig.add_subplot(gs[1, i], sharex=ax_RM)
        ax_lc = fig.add_subplot(gs[2:, i], sharex=ax_RM)

        I, fmask, tmask = Stokes_ds(dat, pol='I')
        Q, fmask, tmask = Stokes_ds(dat, pol='Q')
        U, fmask, tmask = Stokes_ds(dat, pol='U')
        V, fmask, tmask = Stokes_ds(dat, pol='V')

        f = dat['FREQS'] * u.Hz
        df = f[1] - f[0]
        f += df/2 # TODO: Check this!!!

        t = Time(dat['TIMES']/86400, scale='utc', format='mjd')
        dt = t[1] - t[0]

        ephem = get_J1755_ephemeris()
        I[np.isnan(I)] = 0.0
        Q[np.isnan(Q)] = 0.0
        U[np.isnan(U)] = 0.0
        V[np.isnan(V)] = 0.0
        I = dedisperse_ds(I, ephem['DM'], f, safe_nanmean(f), dt)
        Q = dedisperse_ds(Q, ephem['DM'], f, safe_nanmean(f), dt)
        U = dedisperse_ds(U, ephem['DM'], f, safe_nanmean(f), dt)
        V = dedisperse_ds(V, ephem['DM'], f, safe_nanmean(f), dt)
        I[:, fmask] = np.nan
        Q[:, fmask] = np.nan
        U[:, fmask] = np.nan
        V[:, fmask] = np.nan

        for phase_bin in range(len(t)):
            # Get Q and U data from file
            Qcol = Q[phase_bin, :]
            Ucol = U[phase_bin, :]

            # Curve fit RM to data
            p0 = [960, 0.03, -2.4, 0]
            bounds = [(-np.inf, 0, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)]

            try:
                popt, pcov = fit_RM_to_QU(Qcol, Ucol, f, p0=p0, bounds=bounds)
            except:
                continue
            errs = np.sqrt(np.diag(pcov))

            RM, RM_err = popt[0], errs[0]
            L_1GHz, L_1GHz_err = popt[1], errs[1]
            alpha, alpha_err = popt[2], errs[2]
            PA, PA_err = np.divmod(popt[3] + np.pi/2, np.pi)[1] - np.pi/2, errs[3]

            ax_RM.errorbar([t[phase_bin].mjd], [RM], yerr=[RM_err], capsize=2, color='k', fmt='.')
            ax_PA.errorbar([t[phase_bin].mjd], [np.rad2deg(PA)], yerr=[PA_err], capsize=2, color='k', fmt='.')

        # Choose an RM
        RM = 961 * u.rad/u.m**2
        λ = c/f
        ψ = RM*λ**2
        L = Q + U*1j
        Lrot = L*np.exp(-2j*ψ/u.rad)
        Qrot = np.real(Lrot)
        Urot = np.imag(Lrot)
        ax_lc.plot(t.mjd, np.nanmean(I, axis=-1), color='k')
        ax_lc.plot(t.mjd, np.sqrt(safe_nanmean(Qrot, axis=-1)**2 + safe_nanmean(Urot)**2), color='r')

        #plot_RM(ax, f, Qcol, Ucol, popt=popt, pcov=pcov, fscrunch_factor=args.fscrunch_factor)

        ax_RM.set_ylim([955, 967])
        ax_PA.set_ylim([-90, 90])

    if args.output_image is not None:
        plt.tight_layout()
        plt.savefig(args.output_image)
    else:
        plt.show()

if __name__ == '__main__':
    main()
