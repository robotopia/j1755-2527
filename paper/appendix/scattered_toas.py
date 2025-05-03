import numpy as np
from scipy.optimize import root
from scipy.special import erf, erfc
from numpy import pi as π
import matplotlib.pyplot as plt
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

def main():
    parser = argparse.ArgumentParser(description='Make plots for paper appendix')
    parser.add_argument('--output_image', help='File name of output image. Supports same image formats as "plt.savefig()". If not provided, will plt.show().')
    args = parser.parse_args()

    τs = np.logspace(2, 10, 500)
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
    #LEHM_roots = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*erfcxinv(τ/σ*np.sqrt(2/π)) - σ/τ)**2), -0.1).x for τ in τs]).squeeze()
    #LEHM_roots = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*erfcxinv(τ/σ*np.sqrt(2/π)))**2), -0.1).x for τ in τs]).squeeze()
    #LEHM_roots = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*np.sqrt(np.log(τ/σ*np.sqrt(1/(2*π)))) - σ/τ)**2), -0.1).x for τ in τs]).squeeze()
    #LEHM_roots = np.array([root(lambda t: h*σ/τ*np.sqrt(π/2)*np.exp(-0.5*((t-μ)/σ)**2)*erfcx(1/np.sqrt(2)*(σ/τ - (t[0]-μ)/σ)) - 0.5*h*np.exp(-0.5*(np.sqrt(2)*np.sqrt(np.log(τ/σ*np.sqrt(1/(2*π)))))**2), -0.1).x for τ in τs]).squeeze()
    IPLE_roots = np.array([root(lambda t: σ/τ + (t - μ)/σ - (σ/τ)**2*np.sqrt(π/2) * erfcx(1/np.sqrt(2)*(σ/τ - (t[0] - μ)/σ)), -0.1).x for τ in τs]).squeeze()
    MODE_roots = np.array([emg_mode(h, μ, σ, τ)[0] for τ in τs]).squeeze()
    ts = np.linspace(-10, 10, 1000)
    MCHF_roots = np.array([root(lambda lag: matched_filter(ts, lag, h, μ, σ, τ), 1.0).x for τ in τs]).squeeze()

    # Get the rate of change per frequency (relative to the frequency at which τ=σ)
    νs = τs**4
    x = νs; x0 = x[:-1]; x1 = x[1:];
    νs_gm = np.sqrt(x0*x1)

    y = LEHM_roots; y0 = y[:-1]; y1 = y[1:];
    LEHM_roots_slope = np.sqrt(y0*y1/(x0*x1)) * np.log(y1/y0) / np.log(x1/x0)

    y = IPLE_roots; y0 = y[:-1]; y1 = y[1:];
    IPLE_roots_slope = np.sqrt(y0*y1/(x0*x1)) * np.log(y1/y0) / np.log(x1/x0)

    y = MODE_roots; y0 = y[:-1]; y1 = y[1:];
    MODE_roots_slope = np.sqrt(y0*y1/(x0*x1)) * np.log(y1/y0) / np.log(x1/x0)

    y = MCHF_roots; y0 = y[:-1]; y1 = y[1:];
    MCHF_roots_slope = np.sqrt(y0*y1/(x0*x1)) * np.log(y1/y0) / np.log(x1/x0)

    # Plots!

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6,10))

    axs[0].plot(τs, np.abs(LEHM_roots), label="Leading edge at half max (LEHM)")
    axs[0].plot(τs, np.abs(IPLE_roots), label="Inflection point on leading edge (IPLE)")
    axs[0].plot(τs, np.abs(MODE_roots), label="Mode")
    axs[0].plot(τs, np.abs(MCHF_roots), label="Matched filter")
    #axs[0].plot(τs, np.full(τs.shape, np.sqrt(2*np.log(2))), 'k--', alpha=0.2, label="expected max deviation")
    #axs[0].plot(τs, -(np.arctan(np.log(τs)) - π/2)*np.sqrt(2*np.log(2))/π, label="sandbox function (arctan)")
    #axs[0].plot(τs, 1/(τs + 1)*np.sqrt(2*np.log(2)), label="sandbox function (logistic)")
    #axs[0].plot(τs, 4*τs**(-0.97), 'k--', alpha=0.2)

    #######
    # A quick test of asymptotic behaviour of LEHM
    #
    Zs = LEHM_roots #(np.logspace(-10, 0, 1000) - μ)/σ
    arg = Zs/np.sqrt(2)
    #τ_σ = (Zs*erfc(-arg) + np.sqrt(2/π) * np.exp(-arg**2)) / erf(arg)
    τ_σ = 1/np.array([(erfcx(arg[i]) - np.exp(arg[i]**2))/(Zs[i]*erfcx(arg[i]) - 2*Zs[i]*np.exp(arg[i]**2) - 2) for i in range(len(Zs))])
    axs[0].plot(τ_σ, Zs, 'k--')
    #axs[0].plot(τ_σ, 1/τ_σ, 'g--')
    #
    #
    #######

    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    #axs[0].ylim([None, 4])
    axs[0].set_xlabel("$\\tau/\\sigma$")
    axs[0].set_ylabel("$|\\Delta$ToA$|/\\sigma$")
    axs[0].legend()

    axs[1].plot(νs_gm, LEHM_roots_slope, label="Leading edge at half max (LEHM)")
    axs[1].plot(νs_gm, IPLE_roots_slope, label="Left inflection point")
    axs[1].plot(νs_gm, MODE_roots_slope, label="Mode")
    axs[1].plot(νs_gm, MCHF_roots_slope, label="Matched filter")

    axs[1].set_xlabel("$\\nu/\\nu_s$")
    axs[1].set_ylabel("$\\Delta$ToA$/(\\nu\\sigma)$")
    axs[1].set_xscale('log')
    #axs[1].set_yscale('log')

    if args.output_image is not None:
        plt.tight_layout()
        plt.savefig("foo.png")
    else:
        plt.show()

if __name__ == '__main__':
    main()
