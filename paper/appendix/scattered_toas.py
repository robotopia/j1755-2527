import numpy as np
from scipy.optimize import root
from numpy import pi as π
import matplotlib.pyplot as plt
from ctypes import *

so_file = 'erfcxinv.so' # See Makefile for generating this
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

τs = np.logspace(-5, 5, 500)
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

roots = np.array([root(lambda t: emg(t, h, μ, σ, τ) - emg_mode(h, μ, σ, τ)[1]/2, -0.1).x for τ in τs])

#plt.plot(τs, roots, label="numerical \"truth\"")
plt.plot(τs, np.abs(roots), label="Left edge matching")
plt.plot(τs, np.full(τs.shape, np.sqrt(2*np.log(2))), 'k--', alpha=0.2, label="expected max deviation")
#plt.plot(τs, -(np.arctan(np.log(τs)) - π/2)*np.sqrt(2*np.log(2))/π, label="sandbox function (arctan)")
#plt.plot(τs, 1/(τs + 1)*np.sqrt(2*np.log(2)), label="sandbox function (logistic)")
plt.plot(τs, 5*τs**(-1), 'k--', alpha=0.2)

plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-5, 4])
plt.legend()
plt.tight_layout()
plt.savefig("foo.png")
#plt.show()
