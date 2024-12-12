Here is Dougal's plot of which periods were originally ruled out by non-detections in VAST

![[vast/askap-nondetections.png]]

Here is the lightcurve plot from https://dev.pipeline.vast-survey.org/sources/34260695/#measurementsHeader:

![[vast/vast-lightcurve.png]]

Now what I want is a plot showing when each of these observations occurred in J1755-2527's pulse phase. I've downloaded the data containing the start time of each of these observations to [[vast/VAST-Pipeline-Source-ID-34260695.csv]], also obtainable from the above link.

## Re-testing whether any VAST observation should contain a pulse

The point here is to use both barycentering and DM delay calcs to retest whether any VAST observations should have contained a pulse.

Ensure that you have exported an environment variable called `PLANETARY_EPHEMERIS` which points to the location of an ephemeris such as `de430.bsp` or `de440.bsp`. Then, run

```
python fold.py
```

This will produce `fold.png`, a copy of which here:

![[vast/fold_fixed.png]]
Pulse #0 is the bright one. Otherwise, there are only two other observations where the pulse is squarely inside the observation:

1. 2023-04-06
2. 2023-10-18

Looking at the VAST lightcurve above, the 2023-10-18 one (second from the right) looks absolutely empty. But the 2023-04-06 point is tantalisingly higher S/N than the surrounding points. I would like to take a closer look at that one.

---------------------------------

Dougal reports that obs was already checked in detail, and showed up nothing.

In the meantime, I've been hackily trying to get the ASKAP lightcurves myself. For the original ASKAP detection, I managed to extract the lightcurve from the paper directly. It turns out that that particular plot is encoded in SVG. So (following ChatGPT's advice), I opened the PDF in Inkscape and saved out just the plot as a separate SVG file, saved here as [[vast/stae2376_g12.svg]]. Opening this SVG in a browser, and using the developer tools, I identified the path element that corresponds to the Stokes I lightcurve, and cut and paste it to a separate file.

```
echo "0,0 2.813,-4.507 2.813,1.952 2.813,2.432 2.813,-2.398 2.813,3.209 2.813,-1.381 2.813,0.549 2.813,-3.938 2.813,5.535 2.813,-1.58 2.813,-3.052 2.813,2.056 2.813,-2.411 2.813,6.348 2.813,-3.726 2.813,1.595 2.813,-2.093 2.813,-0.582 2.813,-2.078 2.813,2.829 2.813,-0.583 2.813,3 2.813,4.675 2.813,-4.072 2.813,2.716 2.813,-0.219 2.813,0.416 2.813,2.079 2.813,-5.105 2.813,2.019 2.813,-0.444 2.813,2.364 2.813,-0.314 2.813,-2.609 2.813,4.792 2.813,-5.565 2.813,4.135 2.812,-2.02 2.813,1.13 2.813,-0.727 2.813,-0.865 2.813,-0.525 2.813,4.111 2.813,1.67 2.813,-3.856 2.813,6.008 2.813,-0.745 2.813,3.197 2.813,0.577 2.813,3.693 2.813,6.273 2.813,14.166 2.813,28.152 2.813,31.194 2.813,31.868 2.813,25.983 2.813,35.615 2.813,5.785 2.813,-15.115 2.813,-22.056 2.813,-52.573 2.813,-28.369 2.813,-45.126 2.813,-12.88 2.813,-5.803 2.813,-6.94 2.813,-0.428 2.813,-0.152" | tr ' ' '\n' | sed 's/^[^,]*,//' > askap1-mpath.txt
```

The path element is expressed as a relative path, but a bit of NumPy sorts that out:

**parse_askap1.py**
```
import numpy as np
import matplotlib.pyplot as plt

pulse_diff = np.loadtxt('askap1-mpath.txt')
pulse = np.cumsum(pulse_diff) / 1e3  # Because the ULP site expects Jy, not mJy

np.savetxt('askap1-I.txt', pulse)

plt.plot(np.arange(len(pulse))*10, pulse)
plt.xlabel("Time from start of observations (s)")
plt.ylabel("Flux density (Jy)")
plt.savefig('askap1-I.png')
```

![[vast/askap1-I.png]]

The only thing to remember is that each time step is 10 seconds. But happily, the y-axis looks correct!

The metadata needed for uploading to the ULP website (https://ulp.duckdns.org/data/timing/lightcurve_add/10) is summarised here:

> [!ASKAP #1 metadata]
> **Polarisations**: I
> **Telescope**: ASKAP
> **Frequency (MHz)**: 887.5
> **Bandwidth (MHz)**: 288
> **MJD of first sample**: 59965.03588541667 (= 2023-01-21T00:51:40.5 UTC)
> **Duration of samples (s)**: 10
> **DM used (pc/cm^3)**: 0
> **DM referemce frequency (MHz)**: 887.5  (I *think*, since "DM used" is 0, that this doesn't matter...)
> **Values (Jy)**: `askap1-I.txt`

------------------------------------

The second ASKAP pulse is in `j1755_burst2_unsub_lc.csv`, but it's the lightcurve without model subtraction. I'll use this until the subtracted one comes along.

Parse into a format that my website can take:

```
from astropy.time import Time
with open("j1755_burst2_unsub_lc.csv", "r") as f:
    lines = f.readlines()
    prev_t = Time(5227513916.09279/86400.0, scale='utc', format='mjd')
    for line in lines[1:]: # Skip header row
        tokens = line.split(',')
        t = Time(float(tokens[0])/86400.0, scale='utc', format='mjd')
        dt = t - prev_t
        prev_t = t
        I = float(tokens[1])
        print(t.mjd, dt.to('s').value, I)
```

Output saved in [[vast/j1755_burst2_unsub_lc_summary.txt]].

> [!ASKAP #2 Metadata]
> - **Polarisations**: `_ _ I`
> - **Telescope**: ASKAP
> - **Frequency (MHz)**: 887.5
> - **MJD of first sample**: 60503.633288111094
> - **Duration of samples (s)**: 9.95328
> - **DM used (pc/cm^3)**: 0
> - **DM reference frequency (MHz)**: 
