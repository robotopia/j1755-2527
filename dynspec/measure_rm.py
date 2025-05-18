import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.constants import c
import sys
import argparse

from timing import *

import warnings

data = [
    {
        'ds_file': '1358297519_askap.pkl',
        'flip_RM': False,
        'ToA': Time(59965.042934515615, scale='utc', format='mjd'),
        'xlim': [-149, 149],
        'title': "ASKAP, MJD 59965, 888 MHz",
        'do_parallactic_angle_correction': False,
    },
    {
        'ds_file': '1404832334_askap.pkl',
        'flip_RM': False,
        'ToA': Time(60503.63475701395, scale='utc', format='mjd'),
        'xlim': [-149, 149],
        'title': "ASKAP, MJD 60503, 888 MHz",
        'do_parallactic_angle_correction': False,
    },
    {
        'ds_file': '1413381294_meerkat.pkl',
        'flip_RM': True,
        'baseline': [-0.05087631833089632, -0.010672092754814125],
        'ToA': Time(60602.58357327181, scale='utc', format='mjd'),
        'xlim': [-199, 99],
        'title': "MeerKAT, MJD 60602, 813 MHz",
        'do_parallactic_angle_correction': True,
    },
]

# Hardcoded "baseline", derived from pulsestack_all.py

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
    fig = plt.figure(figsize=(4*len(data),7))
    gs = gridspec.GridSpec(5, len(data), figure=fig, hspace=0.0, wspace=0.0)

    threshold = 1.4 # Alpha measurements with errors less than this are included on plots

    axs_RM = []
    axs_PA = []
    axs_lc = []
    axs_al = []

    for i in range(len(data)):
        ds_file = data[i]['ds_file']
        print(f'Opening {ds_file}...')
        dat = np.load(ds_file, allow_pickle=True)

        axs_RM.append(fig.add_subplot(gs[0, i]))
        axs_al.append(fig.add_subplot(gs[1, i]))
        axs_PA.append(fig.add_subplot(gs[2, i]))
        axs_lc.append(fig.add_subplot(gs[3:,i]))

        ax_RM = axs_RM[-1]
        ax_PA = axs_PA[-1]
        ax_lc = axs_lc[-1]
        ax_al = axs_al[-1]

        I, fmask, tmask = Stokes_ds(dat, pol='I')
        Q, fmask, tmask = Stokes_ds(dat, pol='Q')
        U, fmask, tmask = Stokes_ds(dat, pol='U')
        V, fmask, tmask = Stokes_ds(dat, pol='V')

        f = dat['FREQS'] * u.Hz
        df = f[1] - f[0]
        f += df/2 # TODO: Check this!!!

        t = Time(dat['TIMES']/86400, scale='utc', format='mjd')
        dt = t[1] - t[0]
        t_base = (t - data[i]['ToA']).to('s')

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

        # MeerKAT RM seems to be negative of the others. Compensate by flipping, e.g., Q
        if data[i]['flip_RM'] == True:
            Q *= -1

        # MeerKAT also needs to have parallactic angle correction applied to it
        parallactic_angle = 0*u.deg # Default behaviour = assume correction has already been applied
        if data[i]['do_parallactic_angle_correction'] == True:
            location = EarthLocation.of_site(dat['TELESCOPE'])
            lst = t[0].sidereal_time('apparent', longitude=location.lon)
            H = lst - ephem['coord'].ra
            φ = location.lat
            δ = ephem['coord'].dec
            parallactic_angle = np.arctan2(np.sin(H), np.tan(φ)*np.cos(δ) - np.sin(δ)*np.cos(H))
            print(f"Applying parallactic angle of {parallactic_angle.to('deg')}")

        RMs = []
        alphas = []
        alpha_errs = []

        for phase_bin in range(len(t)):
            # Get Q and U data from file
            Qcol = Q[phase_bin, :]
            Ucol = U[phase_bin, :]

            # Curve fit RM to data
            p0 = [961, 0.03, -2.4, 0]
            bounds = [(-np.inf, 0, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf)]

            try:
                popt, pcov = fit_RM_to_QU(Qcol, Ucol, f, p0=p0, bounds=bounds)
            except:
                RMs.append(np.nan)
                alphas.append(np.nan)
                alpha_errs.append(np.nan)
                continue
            errs = np.sqrt(np.diag(pcov))

            RM, RM_err = popt[0], errs[0]
            L_1GHz, L_1GHz_err = popt[1], errs[1]
            alpha, alpha_err = popt[2], errs[2]
            PA, PA_err = np.divmod(popt[3] + np.pi/2, np.pi)[1] - np.pi/2, errs[3]

            RMs.append(RM)
            alphas.append(alpha)
            alpha_errs.append(alpha_err)
            if alpha_err < threshold:
                ax_RM.errorbar([t_base[phase_bin].value], [RM], yerr=[RM_err], capsize=2, color='k', fmt='.')
                ax_al.errorbar([t_base[phase_bin].value], [alpha], yerr=[alpha_err], capsize=2, color='k', fmt='.')

        # Plot only decent PA points
        mask = (np.array(alpha_errs) < threshold) & (np.array(alphas) > -10)
        alphas = np.array(alphas)[mask]
        print(f'avg(α_L) = {np.nanmean(alphas)} ± {np.nanstd(alphas)}')

        # Choose an RM. Weight the individual RMs by profile
        weights = safe_nanmean(I, axis=-1)[mask]
        RM = safe_nanmean(np.array(RMs)[mask] * weights)/np.mean(weights) * u.rad/u.m**2
        print(f'avg(RM) = {RM} ± {np.nanstd(np.array(RMs)[mask])}')
        RM_err = np.nanstd(np.array(RMs)[mask]) * u.rad/u.m**2
        λ = c/f
        ψ = RM*λ**2 - parallactic_angle # Apply correction here
        L = Q + U*1j
        Lrot = L*np.exp(-2j*ψ/u.rad)
        Qrot = np.real(Lrot)
        Urot = np.imag(Lrot)

        ax_RM.axhline(RM.value, c='r', zorder=-100)
        #ax_RM.axhspan((RM - RM_err).value, (RM + RM_err).value, color='r', zorder=-1000, alpha=0.2)

        I_lc = safe_nanmean(I, axis=-1)
        Q_lc = safe_nanmean(Qrot, axis=-1)
        U_lc = safe_nanmean(Urot, axis=-1)
        L_lc = np.sqrt(Q_lc**2 + U_lc**2)
        ψ_lc = 0.5*np.arctan2(U_lc, Q_lc)
        V_lc = safe_nanmean(V, axis=-1)

        # Remove baseline from Stokes I, if needed
        if 'baseline' in data[i].keys():
            phase = (t_base / ephem['period']).decompose()
            M, C = data[i]['baseline']
            I_lc -= M*phase + C

        ax_PA.scatter(t_base[mask], np.rad2deg(ψ_lc[mask]), c='k', s=2)
        ax_PA.scatter(t_base[mask], np.rad2deg(ψ_lc[mask]) + 180*u.deg, c='k', s=2)
        ax_lc.plot(t_base, I_lc*1e3, color='k', label="Total intensity")
        ax_lc.plot(t_base, L_lc*1e3, color='r', label="Lin. pol.")
        ax_lc.plot(t_base, V_lc*1e3, color='b', label="Circ. pol.")

        ax_RM.set_xticks([])
        ax_PA.set_xticks([])
        ax_al.set_xticks([])

        ax_RM.set_title(data[i]['title'])

        '''
        # Quick inspection of the de-rotated dynamic spectrum
        fig_ds, axs_ds = plt.subplots(1, 2)
        axs_ds[0].pcolormesh(Qrot.T)
        axs_ds[1].pcolormesh(Urot.T)
        fig_ds.savefig(f'{ds_file}.QU.png')
        plt.close(fig_ds)
        '''

        #plot_RM(ax, f, Qcol, Ucol, popt=popt, pcov=pcov, fscrunch_factor=args.fscrunch_factor)

        ax_RM.set_ylim([952, 970])
        ax_PA.set_ylim([-90, 270])
        ax_al.set_ylim([-5, 1])
        ax_lc.set_ylim([-12, 99])
        ax_lc.set_xlabel("Time since ToA (s)")

        ax_RM.set_xlim(data[i]['xlim'])
        ax_PA.set_xlim(data[i]['xlim'])
        ax_al.set_xlim(data[i]['xlim'])
        ax_lc.set_xlim(data[i]['xlim'])

    axs_RM[0].set_ylabel("RM (rad/m²)")
    axs_PA[0].set_ylabel("PA (deg)")
    axs_lc[0].set_ylabel("Flux density (mJy)")
    axs_al[0].set_ylabel("Lin. pol.\nspectral index")
    axs_lc[-1].legend()

    axs_PA[0].set_yticks([-90, 0, 90, 180, 270])

    for i in range(1, len(data)):
        for ax in [axs_RM[i], axs_PA[i], axs_lc[i], axs_al[i]]:
            ax.set_yticks([])

    if args.output_image is not None:
        plt.tight_layout()
        plt.savefig(args.output_image)
    else:
        plt.show()

if __name__ == '__main__':
    main()
