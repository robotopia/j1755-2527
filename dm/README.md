# Dispersion measure
Contains code/data relating to measuring the dispersion measure

## In-band measurements

The in-band measurement of the original ASKAP pulse (as reported in [Dobie et al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.535..909D/abstract)), is

$$
{\rm DM} = 710 ^{+200}_{-180} \, {\rm pc}/{\rm cm}^{3}.
$$

Interestingly, the in-band measurement of one of the bright MWA pulses is 1221 ± 257 pc/cm³.

Scattering, which becomes significant at such high DMs (towards this source) at these low frequencies, could explain this discrepancy. In particular, suppose we modelled the intrinsic pulses as Gaussians. Then the scattered pulses could be modelled as exponentially modified Gaussians (EMGs) whose means get shifted (relative to their unscattered counterparts) approximately proportionally to the scattering timescale, $\tau_{\rm sc}$, which itself goes as $\nu^{-4}$ (or $\nu^{-4.4}$). To keep it general, let the spectral index of the scattering timescale be $\beta$, so that the scattering timescale is

$$
\tau_{\rm sc} = \tau_{\rm sc,1 GHz} \bigg( \frac{\nu}{1\,{\rm GHz}} \bigg)^\beta
$$

The total observed delay is then made up of both a contribution from dispersion as well as a contribution from scattering:

$$
\begin{aligned}
    \tau &= \tau_{\rm DM} + \tau_{\rm sc} \\
        &= 4148.8\,{\rm s} \times {\rm DM}\,\bigg(\frac{\nu}{\rm 1 MHz}\bigg)^{-2} + \tau_{\rm sc,1 GHz} \bigg( \frac{\nu}{1\,{\rm GHz}} \bigg)^\beta
\end{aligned}
$$

Electron density models provide a way of relating $\tau_{\rm sc,1 GHz}$ to the DM. Thus, we can find out (to first order), for each model, what scattering timescale and associated spectral index would create the observed slope at 185 MHz that I previously "mistook" for a pure DM delay.

$$
\frac{d\tau}{d\nu} = -8297.6 \, {\rm s}\,{\rm MHz}^{-1} \times{\rm DM} \, \bigg(\frac{\nu}{\rm 1 MHz}\bigg)^{-3} + \beta \, \frac{\tau_{\rm sc,1 GHz}}{\rm 1\,GHz} \bigg( \frac{\nu}{1\,{\rm GHz}} \bigg)^{\beta - 1}.
$$

I've done some calculations in [scattering_calc.py](dm/scattering_calc.py), just to get a ballpark idea of what kind of DM could possible account for the large delay slope measured for the MWA pulse discussed above. I find that for a DM of ~600 MHz, the slope from the DM is comparable to the slope from scattering (assuming an index of -4).
