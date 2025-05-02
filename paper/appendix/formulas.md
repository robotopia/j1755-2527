# Formulas relating to exponentially modified Gaussians (emg)

## Basic properties

Most of the below is summarised on [Wikipedia](https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution).

### Scattered pulse shape

$$
\begin{aligned}
{\rm emg}(t; h, \mu, \sigma, \tau) &=
    \frac{h\sigma}{\tau} \sqrt{\frac{\pi}{2}} \exp \left(
        \frac{1}{2} \left(\frac{\sigma}{\tau}\right)^2 - \frac{t - \mu}{\sigma}
    \right)
    {\rm erfc}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{t - \mu}{\sigma}\right)
    \right) \\
    &= \frac{h\sigma}{\tau} \sqrt{\frac{\pi}{2}} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right)
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{t - \mu}{\sigma}\right)
    \right)
\end{aligned}
$$

### Mode

Assuming $\tau \ge 0$,

$$
t_m = \mu - \sqrt{2}\sigma \, {\rm erfcxinv}\left(\frac{\tau}{\sigma} \sqrt{\frac{2}{\pi}} \right) + \frac{\sigma^2}{\tau}
$$

$$
{\rm emg}(t_m) = h \exp \left(-\frac{1}{2} \left(\frac{t_m - \mu}{\sigma}\right)^2\right)
$$

## Derivatives (and other useful properties)

### erfc

$$
\frac{d}{dt} {\rm erfc}(t) = -\frac{2}{\sqrt{\pi}} \exp(-t^2)
$$

### erfcx

$$
\begin{aligned}
{\rm erfcx}(t) &= \exp(t^2) \cdot {\rm erfc}(t) \\
\frac{d}{dt}{\rm erfcx(t)} &= 2t \exp(t^2) \cdot {\rm erfc}(t) +
\exp(t^2) \cdot \left(-\frac{2}{\sqrt{\pi}} \exp(-t^2) \right) \\
    &= 2t\,{\rm erfcx}(t) - \frac{2}{\sqrt{\pi}}
\end{aligned}
$$

Some useful values and asymptotes:

| | $t \rightarrow -\infty$ | $t = 0$ | $t = 1$ | $t \rightarrow +\infty$ |
| :--: | :--: | :--: | :--: | :--: |
| ${\rm erfcx}(t)$ | $\sim 2 \exp(t^2)$ | $1 - \frac{2t}{\sqrt{\pi}} + \mathcal{O}(t^2)$ |  | $\sim \frac{1}{t\sqrt{\pi}}$ |
| ${\rm erfcxinv}(t)$ | - | $(t \rightarrow 0^+)$: $\sim \frac{1}{t\sqrt{\pi}}$ | $\frac{\sqrt{\pi}(1 - t)}{2} + \mathcal{O}(t^2)$ | $\sim -\sqrt{\ln\left(\frac{t}{2}\right)}$ |

### emg

$$
\begin{aligned}
\frac{d}{dt} {\rm emg}(t) &=
    \frac{h\sigma}{\tau} \sqrt{\frac{\pi}{2}} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) \left(-\frac{t - \mu}{\sigma}\right) \frac{1}{\sigma} \cdot
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{t - \mu}{\sigma}\right)
    \right)  +{} \\ &\quad
    \frac{h\sigma}{\tau} \sqrt{\frac{\pi}{2}} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) \cdot \left[
    2 \left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{t - \mu}{\sigma}\right)
    \right) \,
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{t - \mu}{\sigma}\right)
    \right) - \frac{2}{\sqrt{\pi}}
    \right] \left(-\frac{1}{\sqrt{2}\,\sigma}\right) \\
    &= {\rm emg}(t) \left(-\frac{t - \mu}{\sigma^2}\right) +
        {\rm emg}(t) \left(-\frac{1}{\tau} + \frac{t - \mu}{\sigma^2}\right) +
        \frac{h}{\tau} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
        \right) \\
    &= \frac{h}{\tau} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
        \right) - \frac{1}{\tau} \, {\rm emg}(t)
\end{aligned}
$$

The last line is verified in Sage:
```
var('h mu sigma tau t')

erfcx(x) = exp(x^2)*erfc(x)
z = (t - mu)/sigma
Z = sigma/tau - z
emg(t) = h*sigma/tau * sqrt(pi/2) * exp(-z^2/2) * erfcx(Z/sqrt(2))

bool(diff(emg(t), t) == h/tau * exp(-1/2*z^2) - emg(t)/tau)
```

----------------------

$$
\begin{aligned}
\frac{d^2}{dt^2}{\rm emg}(t)
&= \frac{d}{dt} \left[ \frac{h}{\tau} \exp \left(
    -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) - \frac{1}{\tau} \, {\rm emg}(t) \right] \\
&= \frac{h}{\tau} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) \cdot \left(-\frac{t - \mu}{\sigma}\right)\frac{1}{\sigma} -
    \frac{1}{\tau} \left[ \frac{h}{\tau} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) - \frac{1}{\tau} \, {\rm emg}(t) \right] \\
&= -\frac{h}{\tau\sigma} \exp \left(
    -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) \left( \frac{\sigma}{\tau} + \frac{t - \mu}{\sigma} \right) + \frac{1}{\tau^2} \, {\rm emg}(t) \\
&= -\frac{h}{\tau\sigma} \exp \left(
    -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) \left( \frac{\sigma}{\tau} + \frac{t - \mu}{\sigma} \right) +{} \\
&\qquad\qquad \frac{h\sigma}{\tau^3} \sqrt{\frac{\pi}{2}} \exp \left(
        -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right)
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{t - \mu}{\sigma}\right)
    \right) \\
&= -\frac{h}{\tau\sigma} \exp \left(
    -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) \left[ \left( \frac{\sigma}{\tau} + \frac{t - \mu}{\sigma} \right) - \frac{\sigma^2}{\tau^2} \sqrt{\frac{\pi}{2}}
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{t - \mu}{\sigma}\right)
    \right)
    \right]
\end{aligned}
$$

```
...

bool(diff(emg(t), t, 2) == -h/tau/sigma*exp(-1/2*z^2)*(sigma/tau + z - (sigma/tau)^2*sqrt(pi/2)*erfcx(Z/sqrt(2))))
```

## Defining ToAs on scattered pulses

### Leading edge at half maximum (LEHM)

This method identifies ToAs with the location where the leading edge reaches half of the (scattered) pulse's maximum value. The maximum value is found by the mode, so this ToA is thus defined at the location which satisfies

$$
{\rm emg}({\rm ToA_{LEHM}})
    = \frac{1}{2}{\rm emg}(t_m)
$$

This can solved numerically.

#### Asymptotic behaviour

For small scattering timescales ($\tau \ll \sigma$), ${\rm ToA_{LEHM}}$ approaches $\frac12$FWHM of the unscattered pulse, which occurs at

$$
\lim_{\tau/\sigma \rightarrow 0} {\rm ToA_{LEHM}} = \mu - \sigma\sqrt{2\ln 2}.
$$

##### Proof of asymptotic behaviour

The equation which ${\rm ToA_{LEHM}}$ satisfies can be expanded to

$$
\begin{aligned}
\frac{h\sigma}{\tau} \sqrt{\frac{\pi}{2}} & {} \exp \left(
        -\frac{1}{2} \left(\frac{{\rm ToA_{LEHM}} - \mu}{\sigma}\right)^2
    \right)
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - \frac{{\rm ToA_{LEHM}} - \mu}{\sigma}\right)
    \right) \\
&= \frac{1}{2} h \exp \left(-\frac{1}{2} \left(\sqrt{2}\, {\rm erfcxinv}\left(\frac{\tau}{\sigma} \sqrt{\frac{2}{\pi}} \right) - \frac{\sigma}{\tau}\right)^2\right)
\end{aligned}
$$

Cancelling a few factors, and adopting the shorthand
$$
Z \equiv \frac{{\rm ToA_{LEHM}} - \mu}{\sigma},
$$
the above becomes
$$
\frac{\sigma}{\tau} \sqrt{2\pi} \exp \left(
        -\frac{Z^2}{2} 
    \right)
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - Z\right)
    \right)
= \exp \left(-\frac{1}{2} \left(\sqrt{2}\, {\rm erfcxinv}\left(\frac{\tau}{\sigma} \sqrt{\frac{2}{\pi}} \right) - \frac{\sigma}{\tau}\right)^2\right) 
$$

When $\tau \ll \sigma$, the above asymptotically approaches

$$
\begin{aligned}
2 \exp \left(
        -\frac{Z^2}{2} 
    \right) &= 1 \\
Z^2 &= 2 \ln 2 \\
{\rm ToA_{LEHM}} &= \mu - \sigma\sqrt{2\ln 2},
\end{aligned}
$$
where the sign of the square root was chosen appropriate for the half-maximum point on the leading edge. Since this is a regime in which this definition of ToAs will never be chosen, the first order (constant) behaviour of the asymptote is sufficient.

When $\tau \gg \sigma$, the left hand side of the earlier expression can be expanded about $Z$ to

$$
\begin{aligned}
\frac{\sigma}{\tau} \sqrt{2\pi} & {} \exp \left( -\frac{Z^2}{2} \right)
    {\rm erfcx}\left(
        \frac{1}{\sqrt{2}} \left(\frac{\sigma}{\tau} - Z\right)
    \right) \\
    &\approx \frac{\sigma}{\tau} \sqrt{2\pi} \exp \left( -\frac{Z^2}{2} \right)
    \left[
        {\rm erfcx}\left( -\frac{Z}{\sqrt{2}} \right)
        + \frac{1}{\sqrt{2}}\frac{\sigma}{\tau}\left(-2\frac{Z}{\sqrt{2}} \, {\rm erfcx}\left( -\frac{Z}{\sqrt{2}} \right) - \frac{2}{\sqrt{\pi}} \right)
    \right] \\
    &= \frac{\sigma}{\tau} \sqrt{2\pi} \exp \left( -\frac{Z^2}{2} \right)
    \left[
        {\rm erfcx}\left( -\frac{Z}{\sqrt{2}} \right) \left(
        1 - \frac{\sigma}{\tau}Z \right) - \frac{\sigma}{\tau}\sqrt{\frac{2}{\pi}}
    \right] \\
    &= \frac{\sigma}{\tau} \sqrt{2\pi} 
    \left[
        {\rm erfc}\left( -\frac{Z}{\sqrt{2}} \right) \left(
        1 - \frac{\sigma}{\tau}Z \right) - \frac{\sigma}{\tau}\sqrt{\frac{2}{\pi}} \exp \left( -\frac{Z^2}{2} \right)
    \right]
\end{aligned}
$$

and the right hand side,

$$
\begin{aligned}
\exp \left(-\frac{1}{2} \left(\sqrt{2}\, {\rm erfcxinv}\left(\frac{\tau}{\sigma} \sqrt{\frac{2}{\pi}} \right) - \frac{\sigma}{\tau} \right)^2\right)
    &\approx \exp \left(-\frac{1}{2} \left(\sqrt{2}\, \sqrt{\ln\left(\frac{\tau}{\sigma} \sqrt{\frac{1}{2\pi}} \right)}\right)^2\right) \\
    &= \exp \left(\ln\left(\frac{\sigma}{\tau} \sqrt{2\pi} \right)\right) \\
    &= \frac{\sigma}{\tau} \sqrt{2\pi}
\end{aligned}
$$

Equating both sides, we find

$$
1 = {\rm erfc}\left( -\frac{Z}{\sqrt{2}} \right) \left(
        1 - \frac{\sigma}{\tau}Z \right) - \frac{\sigma}{\tau}\sqrt{\frac{2}{\pi}} \exp \left( -\frac{Z^2}{2} \right)
$$

I'm finding this one a bit hard to solve for $Z$ (and hence, ${\rm ToA_{LEHM}}$), but the inverse (solving for $\tau/\sigma$) is
$$
\begin{aligned}
\frac{\tau}{\sigma}
    &= \frac{Z \, {\rm erfc}\left( -\frac{Z}{\sqrt{2}} \right) + \sqrt{\frac{2}{\pi}} \exp \left( -\frac{Z^2}{2} \right)}{{\rm erf}\left( \frac{Z}{\sqrt{2}} \right)}
\end{aligned}
$$

> [!warning]
> When I plot this up, it doesn't seem to match LEHM. Instead, weirdly, it seems to match IPLE, which in turn seems to just show the asymptotic behaviour ${\rm ToA}/\sigma \approx (\tau/\sigma)^{-1}$.
> 
> Sage confirms my expansion ($x = \sigma/\tau$):
> ```
> bool(taylor(exp(-Z^2/2)*erfcx(1/sqrt(2)*(x-Z)), x, 0, 1) == erfc(-Z/sqrt(2))*(1 - x*Z) - x*sqrt(2/pi)*exp(-Z^2/2))
> ```
> which returns `True`.

### Inflection point on leading edge (IPLE)

\[TODO\]

### Matched filter

Towards a matched filter approach, we compute the convolution of the emg with another Gaussian matched to the unscattered pulse. 

$$
{\rm emg}(t) \ast \exp \left[ -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2 \right]
    = \int_{-\infty}^{\infty} {\rm emg}(t_0) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \, dt_0
$$

where I'm using $t_0$ for the lag parameter. We want the place where this is maximised, so we differentiate with respect to $t$ and set to zero.

$$
\begin{aligned}
0   &= \frac{d}{dt} \int_{-\infty}^{\infty} {\rm emg}(t_0) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \, dt_0 \\
    &= \int_{-\infty}^{\infty} \frac{\partial}{\partial t} \left( {\rm emg}(t_0) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \right) \, dt_0 \\
    &= \int_{-\infty}^{\infty} {\rm emg}(t_0) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \cdot \left( -\frac{(t - t_0) - \mu}{\sigma} \right) \frac{1}{\sigma} \, dt_0
\end{aligned}
$$

The following is a continuation of the above, arriving at an alternative, but equivalent expression:

$$
\begin{aligned}
    &= -\int_{-\infty}^{\infty} {\rm emg}(t_0) \frac{\partial}{\partial t_0} \left( \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \right) \, dt_0 \\
    &= \left[ {\rm emg}(t_0) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \right]_{t_0 = -\infty}^{t_0 = \infty} +{} \\
    &\qquad\qquad \int_{-\infty}^{\infty} \frac{d}{dt_0}\left({\rm emg}(t_0)\right) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \, dt_0 \\
    &= \int_{-\infty}^{\infty} \left( \frac{h}{\tau} \exp \left(
        -\frac{1}{2} \left(\frac{t_0 - \mu}{\sigma}\right)^2
        \right) - \frac{1}{\tau} \, {\rm emg}(t_0) \right)
         \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \, dt_0
\end{aligned}
$$



## Numerical derivatives

In constructing plots using the above, I have used logarithmically spaced NumPy arrays (`np.logspace()`) for the abscissae. Now that I am also interested in plotting the derivatives of these functions, it is easier to just compute the derivatives numerically from their already-derived values rather than 