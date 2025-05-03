import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import sys
sys.path.insert(0, '..')
from scattered_toas import *

xs = np.logspace(-1, 8, 1000)
def expr(x):
    return np.exp(-0.5*(np.sqrt(2)*erfcxinv(x*np.sqrt(2/np.pi)) - 1/x)**2)

plt.plot(xs, [expr(x) for x in xs], label="$\\exp\\left(-\\frac{1}{2}\\left( \\sqrt{2} {\\rm erfcxinv}\\left(x\\sqrt{\\frac{2}{\\pi}}\\right) - \\frac{1}{x}\\right)^2\\right)$")
plt.plot(xs, np.sqrt(2*np.pi)/xs, label="$\\frac{\\sqrt{2\\pi}}{x}$")

def expr2(x, Z):
    return x*np.sqrt(2*np.pi)*erfcx(1/np.sqrt(2) * (x-Z))

def expr3(x, Z):
    e1 = np.sqrt(2*np.pi)*np.exp(Z**2/2)*(erfc(Z/np.sqrt(2)) - 2)
    e2 = x**2*(e1*Z - 2) - x*e1
    return e2

#for Z in [2]:
#    plt.plot(xs, [expr2(1/x, Z) for x in xs], label="Truth")
#    plt.plot(xs, [expr3(1/x, Z) for x in xs], label="Approx")

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$x$')
plt.ylabel('$f(x)$')
#plt.ylim([None, 2])
plt.show()

