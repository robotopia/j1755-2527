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