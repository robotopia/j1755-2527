import numpy as np
from numpy import cos, sin, exp, pi
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import yaml
import astropy.units as u
from astropy.constants import c

import argparse
import os
from importlib.resources import files
import dedisperse_dynspec as dd
from pa_uncertainties import pa_uncertainty

def main():
    parser = argparse.ArgumentParser(description="Plot polarisation profile")
    parser.add_argument('yaml_files', nargs=4, help="Iyaml Qyaml Uyaml Vyaml")
    parser.add_argument('--RM', type=float, help='Defaraday rotate first to this rotation measure (in rad/m²)')
    parser.add_argument('--DM', type=float, help='Dedisperse first using this dispersion measure (in pc/cm³)')
    parser.add_argument('--off_pulse_lims', help='Comma-separated list of inclusive time range to count as "off-pulse" for noise estimation purposes. Each time range is in the format X:Y. E.g. "0:40,125:140". Interpreted as absolute times (e.g. GPS seconds) unless --rel_time is used. If provided, then L and ψ (and their errors) are calculated as discussed in The Pulsar Handbook Section 7.4.3 Polarisation profiles. If not provided, then L=√(Q²+U²) and ψ=0.5*atan(U/Q) are used.')
    parser.add_argument('--outfile', help='Save plot to this file (default is to plt.show())')
    parser.add_argument('--outdynspec', help='Save de-faraday rotated spectra plots to this file (default is to *not* make this plot)')
    parser.add_argument('--rel_time', action='store_true', help='Display times relative to beginning of observation (default is to use absolute time as given in the YAML file)')
    parser.add_argument('--pdv', help='Output the polarisation profile to a pdv-style ASCII text file (where pdv is the name of the utility from PSRCHIVE)')

    args = parser.parse_args()

    with open(args.yaml_files[0], 'r') as yI:
        Iparams = dd.parse_yaml(yI)

    with open(args.yaml_files[1], 'r') as yQ:
        Qparams = dd.parse_yaml(yQ)

    with open(args.yaml_files[2], 'r') as yU:
        Uparams = dd.parse_yaml(yU)

    with open(args.yaml_files[3], 'r') as yV:
        Vparams = dd.parse_yaml(yV)

    I_flag_bins = Iparams['mask_time_bins'] if 'mask_time_bins' in Iparams.keys() else []
    Q_flag_bins = Qparams['mask_time_bins'] if 'mask_time_bins' in Qparams.keys() else []
    U_flag_bins = Uparams['mask_time_bins'] if 'mask_time_bins' in Uparams.keys() else []
    V_flag_bins = Vparams['mask_time_bins'] if 'mask_time_bins' in Vparams.keys() else []
    I_flag_chans = Iparams['mask_freq_bins'] if 'mask_freq_bins' in Iparams.keys() else []
    Q_flag_chans = Qparams['mask_freq_bins'] if 'mask_freq_bins' in Qparams.keys() else []
    U_flag_chans = Uparams['mask_freq_bins'] if 'mask_freq_bins' in Uparams.keys() else []
    V_flag_chans = Vparams['mask_freq_bins'] if 'mask_freq_bins' in Vparams.keys() else []

    flag_bins = list(set(I_flag_bins) | set(Q_flag_bins) | set(U_flag_bins) | set(V_flag_bins))
    flag_chans = list(set(I_flag_chans) | set(Q_flag_chans) | set(U_flag_chans) | set(V_flag_chans))

    # Force input (i.e. the dynspec data files) to be absolute path pointing to the same directory as yaml files
    Iparams['input'] = os.path.join(os.path.dirname(args.yaml_files[0]), Iparams['input'])
    Qparams['input'] = os.path.join(os.path.dirname(args.yaml_files[1]), Qparams['input'])
    Uparams['input'] = os.path.join(os.path.dirname(args.yaml_files[2]), Uparams['input'])
    Vparams['input'] = os.path.join(os.path.dirname(args.yaml_files[3]), Vparams['input'])

    Idynspec = dd.Dynspec(**Iparams)
    Qdynspec = dd.Dynspec(**Qparams)
    Udynspec = dd.Dynspec(**Uparams)
    Vdynspec = dd.Dynspec(**Vparams)

    if args.DM:
        Idynspec.set_freq_ref('centre')
        Qdynspec.set_freq_ref('centre')
        Udynspec.set_freq_ref('centre')
        Vdynspec.set_freq_ref('centre')

        Idynspec.dedisperse(args.DM)
        Qdynspec.dedisperse(args.DM)
        Udynspec.dedisperse(args.DM)
        Vdynspec.dedisperse(args.DM)

    if args.RM:
        RM = args.RM * u.m**2  # The radians are implicit
        f = Qdynspec.f * u.MHz
        L = Qdynspec.dynspec + Udynspec.dynspec*1j
        Lfr = L * np.exp(-2j*(args.RM*(c/f[:,np.newaxis])**2).decompose().value)
        Qdynspec.dynspec = np.real(Lfr)
        Udynspec.dynspec = np.imag(Lfr)

    if args.outdynspec:
        fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,5), sharex=True, sharey=True)
        Idynspec.plot(axs[0])
        Qdynspec.plot(axs[1])
        Udynspec.plot(axs[2])
        Vdynspec.plot(axs[3])
        plt.tight_layout()
        plt.savefig(args.outdynspec)
        plt.close(fig)

    # Frequency scrunch everything
    Idynspec.fscrunch()
    Qdynspec.fscrunch()
    Udynspec.fscrunch()
    Vdynspec.fscrunch()

    # Form plot
    t = Idynspec.t
    f = Idynspec.f

    # Use relative times, if requested
    if args.rel_time:
        t -= t[0]

    I = Idynspec.fscrunched
    Q = Qdynspec.fscrunched
    U = Udynspec.fscrunched
    V = Vdynspec.fscrunched

    L = np.sqrt(Q**2 + U**2)
    ψ = np.rad2deg(0.5*np.angle(Q + U*1j)) # (in deg)

    ΔQ = np.std(Qdynspec.dynspec, axis=Qdynspec.FREQAXIS)
    ΔU = np.std(Udynspec.dynspec, axis=Udynspec.FREQAXIS)

    off_pulse = []
    if args.off_pulse_lims:
        # Parse given time ranges
        for timerange_str in args.off_pulse_lims.split(','):
            starttime_str, endtime_str = timerange_str.split(':')
            starttime, endtime = float(starttime_str), float(endtime_str)
            mask = np.logical_and(t >= starttime, t <= endtime)
            off_pulse += list(I[mask])
        σI = np.std(off_pulse)
        Ltrue = np.zeros(L.shape)
        pa_mask = (L/σI >= 1.57)
        Ltrue[pa_mask] = σI*np.sqrt((L[pa_mask]/σI)**2 - 1)
        P0 = Ltrue / σI
        Δψ = np.rad2deg(pa_uncertainty(P0))
    else:
        Ltrue = L # Ltrue is a misnomer now, but can't do any better without off-pulse stats given
        pa_mask = (Ltrue >= 0) # i.e. keep everything

        # Calculate PA error
        ΔL = np.sqrt(ΔQ**2 + ΔU**2)/len(f) # First order estimate only
        Δψ = 360*(ΔL/L) # (in deg)

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    axs[1].plot(t, I, 'k', label="I")
    axs[1].plot(t, Ltrue, 'r', label="|L|")
    axs[1].plot(t, V, 'b', label="V")

    axs[0].errorbar(t[pa_mask], ψ[pa_mask], yerr=Δψ[pa_mask], fmt='k.')
    axs[0].errorbar(t[pa_mask], ψ[pa_mask]+180, yerr=Δψ[pa_mask], fmt='k.')
    axs[0].set_ylim([-90, 270])

    axs[0].set_ylabel("PA (deg)")
    axs[1].set_xlabel("Time (s)")
    axs[1].set_ylabel("Flux density (Jy)")
    axs[1].legend()
    axs[0].set_yticks([-90, 0, 90, 180, 270])

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

    if args.pdv:
        Ltrue[np.logical_not(pa_mask)] = np.NAN
        ψ[np.logical_not(pa_mask)] = np.NAN
        Δψ[np.logical_not(pa_mask)] = np.NAN
        pdv_data = np.array([
            np.zeros(I.shape), # "subint"
            np.zeros(I.shape), # "frequency" (but only f-scrunched output is supported here)
            np.arange(len(I)), # "phase bin number"
            I, Q, U, V,
            Ltrue, ψ, Δψ
        ])
        header = "# subint chan bin I Q U V L PA PA_err"
        np.savetxt(args.pdv, pdv_data.T, header=header, fmt="%.5e")

if __name__ == "__main__":
    main()
