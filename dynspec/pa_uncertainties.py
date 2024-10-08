from numpy import sqrt, exp, cos, linspace, arange, cumsum, nonzero
from numpy import empty, rad2deg, deg2rad, save, array, load, nan, logical_and
from numpy import pi as π
from scipy.special import erf
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from importlib.resources import files

def precalc_pa_unc():
    # From Everett and Weisberg (2001):
    # https://articles.adsabs.harvard.edu/pdf/1993A%26A...274..968N

    n = 100000  # This effectively sets the precision of the error
    dψ = (π/2) / n
    ψ = (arange(n) + 0.5)*dψ
    C = sqrt(0.5)*cos(2*ψ)

    P0s = linspace(0, 10, 1000)

    Δψ1 = empty(P0s.shape)
    Δψ2 = empty(P0s.shape)
    Δψ3 = empty(P0s.shape)

    for i in range(len(P0s)):
        P0 = P0s[i]
        η = C*P0

        G = exp(-0.5*P0**2)/sqrt(π) * (1/sqrt(π) + η*exp(η**2)*(1 + erf(η)))
        G *= π/n  # This effectively normalises the (discrete) prob. distr. so that it would sum to 1 over the range 0 to π/2

        intG = cumsum(G)

        σ1_idx = nonzero(intG > 0.6827)[0][0]
        σ2_idx = nonzero(intG > 0.9545)[0][0]
        σ3_idx = nonzero(intG > 0.9973)[0][0]

        Δψ1[i] = ψ[σ1_idx]
        Δψ2[i] = ψ[σ2_idx]
        Δψ3[i] = ψ[σ3_idx]

    # Save out a plot
    plt.plot(P0s, rad2deg(Δψ1), label="1σ")
    plt.plot(P0s, rad2deg(Δψ2), label="2σ")
    plt.plot(P0s, rad2deg(Δψ3), label="3σ")
    #plt.plot(P0s, Δψ3/Δψ1, label="3σ/1σ")
    #plt.plot(P0s, Δψ2/Δψ1, label="2σ/1σ")

    plt.xlabel("$P_0$")
    plt.ylabel("Δψ (deg)")

    plt.legend()

    plt.savefig("pa_uncertainties.png")

    # Save out numpy arrays
    save("1sigma.npy", array([P0s, Δψ1]))
    save("2sigma.npy", array([P0s, Δψ2]))
    save("3sigma.npy", array([P0s, Δψ3]))


def pa_uncertainty(P0, sigma=1):
    try:
        unc = load(f'{sigma}sigma.npy')
    except:
        raise(f"{sigma} is not a supported confidence interval. Must be 1, 2, or 3")

    Δψ = empty(P0.shape)
    Δψ[P0 > 10] = deg2rad(28.65)/P0[P0 > 10]
    Δψ[P0 < 0] = nan

    # Mask for the other P0s
    mask = logical_and(P0 >= 0, P0 <= 10)
    P0s, Δψs = unc
    func = interp1d(P0s, Δψs)
    Δψ[mask] = func(P0[mask])
    return Δψ

if __name__ == "__main__":
    precalc_pa_unc()
