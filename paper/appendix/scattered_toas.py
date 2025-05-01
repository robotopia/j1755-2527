import numpy as np
from scipy.optimize import root
from numpy import pi as π
import matplotlib.pyplot as plt
from ctypes import *
import os

so_file = os.path.join(os.path.dirname(__file__), 'erfcxinv.so')

erf_lib = CDLL(so_file)
erf_lib.erfcx.argtypes = [c_double]
erf_lib.erfcx.restype = c_double
erf_lib.erfcxinv.argtypes = [c_double]
erf_lib.erfcxinv.restype = c_double

def emg(t, h, μ, σ, τ):
    z = (t - μ)/σ
    args = 1/np.sqrt(2) * (σ/τ - z)
    if not isinstance(args, np.ndarray):
        args = [args]
    return h * np.exp(-0.5*z**2) * σ/τ * np.sqrt(π/2) * np.array([erf_lib.erfcx(arg) for arg in args])

def emg_mode(h, μ, σ, τ):
    arg = τ/σ * np.sqrt(2/π)
    t_m = μ - np.sqrt(2)*σ * erf_lib.erfcxinv(arg) + σ**2/τ
    return t_m, emg(t_m, h, μ, σ, τ)

τs = np.logspace(-3, 3, 500)
σ = 1
h = 1
μ = 0

'''
ts = np.linspace(-6, 6, 1000)
#plt.plot(ts, [emg(t, h, μ, σ, 2) for t in ts])
plt.plot(ts, [erf_lib.erfcxinv(t) for t in ts])
#plt.xscale('log')
plt.yscale('log')
plt.show()
'''

# For the matched filter approach, we do it numerically
def matched_filter(ts, lag, h, μ, σ, τ):
    '''
    When comparing with the derivation in the notes, be aware that
        t -> t_0   and lag -> t
    '''
    z = (ts - μ)/σ
    z0 = (lag - ts - μ)/σ
    integrand = (h/τ * np.exp(-0.5*z**2) - emg(ts, h, μ, σ, τ)/τ) * np.exp(-0.5*z0**2)

    # Here we just hope that things are centred enough so that the contribution
    # to the integral outside the given range of ts is negligible. We clearly are also
    # neglecting the factor of dt, so the result is dimensionally incorrect. However,
    # we're just using this to find a root, so none of this should matter.
    return np.sum(integrand)

LAHM_roots = np.array([root(lambda t: emg(t, h, μ, σ, τ) - emg_mode(h, μ, σ, τ)[1]/2, -0.1).x for τ in τs]).squeeze()
INFL_roots = np.array([root(lambda t: σ/τ + (t - μ)/σ - (σ/τ)**2*np.sqrt(π/2) * erf_lib.erfcx(1/np.sqrt(2)*(σ/τ - (t[0] - μ)/σ)), -0.1).x for τ in τs]).squeeze()
MODE_roots = np.array([emg_mode(h, μ, σ, τ)[0] for τ in τs]).squeeze()
ts = np.linspace(-10, 10, 1000)
MCHF_roots = np.array([root(lambda lag: matched_filter(ts, lag, h, μ, σ, τ), 1.0).x for τ in τs]).squeeze()

plt.plot(τs, np.abs(LAHM_roots), label="Left at half max")
plt.plot(τs, np.abs(INFL_roots), label="Left inflection point")
plt.plot(τs, np.abs(MODE_roots), label="Mode")
plt.plot(τs, np.abs(MCHF_roots), label="Matched filter")
#plt.plot(τs, np.full(τs.shape, np.sqrt(2*np.log(2))), 'k--', alpha=0.2, label="expected max deviation")
#plt.plot(τs, -(np.arctan(np.log(τs)) - π/2)*np.sqrt(2*np.log(2))/π, label="sandbox function (arctan)")
#plt.plot(τs, 1/(τs + 1)*np.sqrt(2*np.log(2)), label="sandbox function (logistic)")
#plt.plot(τs, 4*τs**(-0.97), 'k--', alpha=0.2)

plt.xscale('log')
plt.yscale('log')
plt.ylim([None, 4])
plt.xlabel("$\\tau/\\sigma$")
plt.ylabel("$\\Delta$ToA$/\\sigma$")
plt.legend()
plt.tight_layout()
plt.savefig("foo.png")
#plt.show()
