# Formulas relating to exponentially modified Gaussians (emg)

## PDF

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

## Derivatives

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

--------------------

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
    &= \int_{-\infty}^{\infty} {\rm emg}(t_0) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \cdot \left( -\frac{(t - t_0) - \mu}{\sigma} \right) \frac{1}{\sigma} \, dt_0 \\
    &= -\int_{-\infty}^{\infty} {\rm emg}(t_0) \frac{\partial}{\partial t_0} \left( \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \right) \, dt_0 \\
    &= \left[ {\rm emg}(t_0) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \right]_{t_0 = -\infty}^{t_0 = \infty} +{} \\
    &\qquad\qquad \int_{-\infty}^{\infty} \frac{d}{dt_0}\left({\rm emg}(t_0)\right) \exp \left[ -\frac{1}{2} \left(\frac{(t - t_0) - \mu}{\sigma}\right)^2 \right] \, dt_0
\end{aligned}
$$