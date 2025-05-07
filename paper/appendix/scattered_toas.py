import numpy as np
from scipy.optimize import root
from scipy.special import erf, erfc
from numpy import pi as π
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ctypes import *
import os
import argparse

def load_erfcx(so_file=os.path.join(os.path.dirname(__file__), 'erfcxinv.so')):

    erf_lib = CDLL(so_file)
    erf_lib.erfcx.argtypes = [c_double]
    erf_lib.erfcx.restype = c_double
    erf_lib.erfcxinv.argtypes = [c_double]
    erf_lib.erfcxinv.restype = c_double

    return erf_lib

erf_lib = load_erfcx()
erfcx, erfcxinv = erf_lib.erfcx, erf_lib.erfcxinv

def emg(t, h, μ, σ, τ):
    z = (t - μ)/σ
    args = 1/np.sqrt(2) * (σ/τ - z)
    if not isinstance(args, np.ndarray):
        args = [args]
    return h * np.exp(-0.5*z**2) * σ/τ * np.sqrt(π/2) * np.array([erfcx(arg) for arg in args])

def emg_mode(h, μ, σ, τ):
    arg = τ/σ * np.sqrt(2/π)
    t_m = μ - np.sqrt(2)*σ * erfcxinv(arg) + σ**2/τ
    return t_m, emg(t_m, h, μ, σ, τ)

# For the matched filter approach, we do it numerically
def matched_filter(ts, lag, h, μ, σ, τ):
    '''
    When comparing with the derivation in the notes, be aware that
        t -> t_0   and lag -> t
    '''
    z0 = (lag - ts - μ)/σ
    integrand = emg(ts, h, μ, σ, τ) * np.exp(-0.5*z0**2) * (-z0) / σ

    # Alternative, equivalent formulation
    #z = (ts - μ)/σ
    #integrand = (h/τ * np.exp(-0.5*z**2) - emg(ts, h, μ, σ, τ)/τ) * np.exp(-0.5*z0**2)

    # Here we just hope that things are centred enough so that the contribution
    # to the integral outside the given range of ts is negligible. We clearly are also
    # neglecting the factor of dt, so the result is dimensionally incorrect. However,
    # we're just using this to find a root, so none of this should matter.
    return np.sum(integrand)

def slope_for_logspace_data(x, y):
    '''
    Should only use this if x is logarithmically spaced.
    '''
    x0 = x[:-1]
    x1 = x[1:]
    x_gm = np.sqrt(x0*x1) # Geometric mean ("mid"-points of logspaced data)

    y0 = y[:-1]
    y1 = y[1:]
    y_gm = np.sqrt(y0*y1)
    return x_gm, y_gm/x_gm * np.log(y1/y0) / np.log(x1/x0)

def main():
    parser = argparse.ArgumentParser(description='Make plots for paper appendix')
    parser.add_argument('--output_image', help='File name of output image. Supports same image formats as "plt.savefig()". If not provided, will plt.show().')
    args = parser.parse_args()

    τs = np.logspace(-4, 4, 500)
    σ = 1
    h = 1
    μ = 0

    '''
    ts = np.linspace(-6, 6, 1000)
    #plt.plot(ts, [emg(t, h, μ, σ, 2) for t in ts])
    plt.plot(ts, [erfcxinv(t) for t in ts])
    #plt.xscale('log')
    plt.yscale('log')
    plt.show()
    '''

    LEHM_roots = np.array([root(lambda t: emg(t, h, μ, σ, τ) - 0.5*emg_mode(h, μ, σ, τ)[1], -0.1).x for τ in τs]).squeeze()
    #LEHM_roots_approx = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*erfcxinv(τ/σ*np.sqrt(2/π)) - σ/τ)**2), -0.1).x for τ in τs]).squeeze()
    #LEHM_roots_approx = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*erfcxinv(τ/σ*np.sqrt(2/π)))**2), -0.1).x for τ in τs]).squeeze()
    #LEHM_roots_approx = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*np.sqrt(np.log(τ/σ*np.sqrt(1/(2*π)))) - σ/τ)**2), -0.1).x for τ in τs]).squeeze()
    #LEHM_roots_approx = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*np.sqrt(np.log(τ/σ*np.sqrt(1/(2*π)))))**2), -0.1).x for τ in τs]).squeeze()
    IPLE_roots = np.array([root(lambda t: σ/τ + (t - μ)/σ - (σ/τ)**2*np.sqrt(π/2) * erfcx(1/np.sqrt(2)*(σ/τ - (t[0] - μ)/σ)), -0.1).x for τ in τs]).squeeze()
    PEAK_roots = np.array([emg_mode(h, μ, σ, τ)[0] for τ in τs]).squeeze()
    ts = np.linspace(-10, 10, 1000)
    TMPL_roots = np.array([root(lambda lag: matched_filter(ts, lag, h, μ, σ, τ), 1.0).x for τ in τs]).squeeze()

    # Get the rate of change per frequency (relative to the frequency at which τ=σ)
    νs = τs**-0.25

    _, LEHM_roots_slope = slope_for_logspace_data(νs, LEHM_roots)
    _, IPLE_roots_slope = slope_for_logspace_data(νs, IPLE_roots)
    _, PEAK_roots_slope = slope_for_logspace_data(νs, PEAK_roots)
    νs_gm, TMPL_roots_slope = slope_for_logspace_data(νs, TMPL_roots)

    # Equivalent DM from above slopes. Choose ν_s = 200 MHz and σ = 1 min. Will change later as needed
    import astropy.units as u
    σ_ref = 25 * u.s
    ν_s = (((70*u.ms)/σ_ref)**(0.25)*1*u.GHz).to('MHz')
    print(f'{ν_s = }')
    D = 4148.808 * u.MHz**2 / u.pc * u.cm**3 * u.s

    LEHM_DM = (LEHM_roots * σ_ref*νs**2*ν_s**2/D).to('pc cm-3')
    IPLE_DM = (IPLE_roots * σ_ref*νs**2*ν_s**2/D).to('pc cm-3')
    PEAK_DM = (PEAK_roots * σ_ref*νs**2*ν_s**2/D).to('pc cm-3')
    TMPL_DM = (TMPL_roots * σ_ref*νs**2*ν_s**2/D).to('pc cm-3')

    LEHM_DM_slope = (0.5 * νs_gm**3 * LEHM_roots_slope * σ_ref * ν_s**2 / D).to('pc cm-3')
    IPLE_DM_slope = (0.5 * νs_gm**3 * IPLE_roots_slope * σ_ref * ν_s**2 / D).to('pc cm-3')
    PEAK_DM_slope = (0.5 * νs_gm**3 * PEAK_roots_slope * σ_ref * ν_s**2 / D).to('pc cm-3')
    TMPL_DM_slope = (0.5 * νs_gm**3 * TMPL_roots_slope * σ_ref * ν_s**2 / D).to('pc cm-3')
    #                ^^^^^ Meaning of this is dimensionless ratios ν/ν_s, evaluated at "mid-way" points

    # Plots!
    fig = plt.figure(figsize=(10, 7))
    gs = GridSpec(2, 2, hspace=0)
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[1,0], sharex=ax0)
    ax2 = fig.add_subplot(gs[0,1], sharex=ax0)
    ax3 = fig.add_subplot(gs[1,1], sharex=ax0)

    # Define some colors, linestyles, etc
    colors = {'LEHM': 'tab:blue', 'IPLE': 'tab:blue', 'PEAK': 'tab:red', 'TMPL': 'tab:red'} # Using color to indicate whether ΔToA is positive or negative
    linestyles = {'LEHM': '--', 'IPLE': '-', 'PEAK': '-', 'TMPL': '--'} # Using linestyle to indicate whether ΔToA is better or worse, for the given sign
    plot_styles = {
        'LEHM': {'c': 'tab:blue', 'ls': '--', 'label': 'Leading edge at half max (LEHM)'},
        'IPLE': {'c': 'tab:blue', 'ls': '-' , 'label': 'Inflection point on leading edge (IPLE)'},
        'PEAK': {'c': 'tab:red',  'ls': '-' , 'label': 'Peak position (PEAK)'},
        'TMPL': {'c': 'tab:red',  'ls': '--', 'label': 'Matched filter (TMPL)'},
    }

    ax0.plot(νs, np.abs(LEHM_roots), **plot_styles['LEHM'])
    ax0.plot(νs, np.abs(IPLE_roots), **plot_styles['IPLE'])
    ax0.plot(νs, np.abs(PEAK_roots), **plot_styles['PEAK'])
    ax0.plot(νs, np.abs(TMPL_roots), **plot_styles['TMPL'])
    #ax0.plot(νs, τs, 'k')
    #ax0.plot(νs, np.sqrt(np.log(τs**2/2/π)) + 1/τs, 'g')
    ax0.plot(νs, np.exp(0.253314137315500*np.log(τs) + 0.110222420968151*np.log(τs)**2), 'g')
    #print(f'{erfcxinv(1) = }')
    #print(f'{erfcx(0) = }')
    #ax0.plot(τs, np.full(τs.shape, np.sqrt(2*np.log(2))), 'k--', alpha=0.2, label="expected max deviation")
    #ax0.plot(τs, -(np.arctan(np.log(τs)) - π/2)*np.sqrt(2*np.log(2))/π, label="sandbox function (arctan)")
    #ax0.plot(τs, 1/(τs + 1)*np.sqrt(2*np.log(2)), label="sandbox function (logistic)")
    #ax0.plot(τs, 4*τs**(-0.97), 'k--', alpha=0.2)

    #######
    # A quick test of asymptotic behaviour of LEHM
    #
    #Zs = -LEHM_roots #(np.logspace(-10, 0, 1000) - μ)/σ
    #arg = Zs/np.sqrt(2)
    #τ_σ = (Zs*erfc(-arg) + np.sqrt(2/π) * np.exp(-arg**2)) / erf(arg)
    #τ_σ = 1/np.array([(erfcx(arg[i]) - np.exp(arg[i]**2))/(Zs[i]*erfcx(arg[i]) - 2*Zs[i]*np.exp(arg[i]**2) - 2) for i in range(len(Zs))])
    #ax0.plot(τ_σ, Zs, 'k--', label="Asymptote (LEHM)")
    #ax0.plot(τ_σ, 1/τ_σ, 'g--')
    #
    #
    #######

    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylim([1e-5, 1e2])
    #ax0.set_xlim([1e-4, 1e4])
    #ax0.set_xlabel("$\\tau/\\sigma$")
    ax0.set_ylabel("$\\dfrac{|\\Delta{\\rm ToA}|}{\\sigma}$")
    #ax0.legend()

    ax1.plot(νs_gm, np.abs(LEHM_roots_slope), **plot_styles['LEHM'])
    ax1.plot(νs_gm, np.abs(IPLE_roots_slope), **plot_styles['IPLE'])
    ax1.plot(νs_gm, np.abs(PEAK_roots_slope), **plot_styles['PEAK'])
    ax1.plot(νs_gm, np.abs(TMPL_roots_slope), **plot_styles['TMPL'])

    #ax1.set_xlabel("$\\nu/\\nu_s$")
    ax1.set_ylabel("$\\left.d\\left|\\dfrac{\\Delta{\\rm ToA}}{\\sigma}\\right|\\middle/d\\left(\\dfrac{\\nu_s}{\\nu}\\right)\\right.$")
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    #ax1.set_xlim([0.1, 10])
    #ax1.set_ylim([-1, 1])
    ax1.legend()

    ax2.plot(νs, np.abs(LEHM_DM), **plot_styles['LEHM'])
    ax2.plot(νs, np.abs(IPLE_DM), **plot_styles['IPLE'])
    ax2.plot(νs, np.abs(PEAK_DM), **plot_styles['PEAK'])
    ax2.plot(νs, np.abs(TMPL_DM), **plot_styles['TMPL'])

    ax2.set_ylim([-20, 370])
    ax2.set_ylabel("$|{\\rm DM}_{\\rm hi-freq}|$ (${\\rm pc}\,{\\rm cm}^{-3}$)")

    ax3.plot(νs_gm, np.abs(LEHM_DM_slope), **plot_styles['LEHM'])
    ax3.plot(νs_gm, np.abs(IPLE_DM_slope), **plot_styles['IPLE'])
    ax3.plot(νs_gm, np.abs(PEAK_DM_slope), **plot_styles['PEAK'])
    ax3.plot(νs_gm, np.abs(TMPL_DM_slope), **plot_styles['TMPL'])

    ax1.set_xlabel("$\\nu/\\nu_s$")
    ax3.set_xlabel("$\\nu/\\nu_s$")
    ax3.set_ylabel("$|{\\rm DM}_{\\rm in-band}|$ (${\\rm pc}\,{\\rm cm}^{-3}$)")
    ax3.set_xscale('log')
    #ax3.set_yscale('log')
    ax3.set_xlim([0.1, 10])
    #ax3.set_ylim([-500, 500])

    ax0.text(0.97, 0.95, "(a)", transform=ax0.transAxes, va='top', ha='right')
    ax1.text(0.97, 0.95, "(c)", transform=ax1.transAxes, va='top', ha='right')
    ax2.text(0.97, 0.95, "(b)", transform=ax2.transAxes, va='top', ha='right')
    ax3.text(0.97, 0.95, "(d)", transform=ax3.transAxes, va='top', ha='right')

    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    def forward(x): return x**-4
    def inverse(x): return x**-0.25
    secax0 = ax0.secondary_xaxis('top', functions=(forward, inverse))
    secax0.set_xlabel("$\\tau/\\sigma$")
    secax2 = ax2.secondary_xaxis('top', functions=(forward, inverse))
    secax2.set_xlabel("$\\tau/\\sigma$")

    if args.output_image is not None:
        plt.tight_layout()
        plt.savefig(args.output_image)
    else:
        plt.show()

if __name__ == '__main__':
    main()
