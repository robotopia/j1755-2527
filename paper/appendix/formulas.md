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

$$
\begin{aligned}
\frac{d^2}{dt^2}{\rm emg}(t)
&= \frac{d}{dt} \left[ \frac{h}{\tau} \exp \left(
    -\frac{1}{2} \left(\frac{t - \mu}{\sigma}\right)^2
    \right) - \frac{1}{\tau} \, {\rm emg}(t) \right] \\
&= 
\end{aligned}
$$
