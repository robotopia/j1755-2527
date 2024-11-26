# Running scripts

Ensure that you have an environment variable called `PLANETARY_EPHEMERIS` which points to the location of an ephemeris such as `de430.bsp` or `de440.bsp`. Then, run

```
python fold.py
```

This will produce `fold.png`, a copy of which here:

![[vast/fold_fixed.png]]
# Discussion

Here is Dougal's plot of which periods were originally ruled out by non-detections in VAST

![[vast/askap-nondetections.png]]

Here is the lightcurve plot from https://dev.pipeline.vast-survey.org/sources/34260695/#measurementsHeader:

![[vast/vast-lightcurve.png]]

Now what I want is a plot showing when each of these observations occurred in J1755-2527's pulse phase. I've downloaded the data containing the start time of each of these observations to [[vast/VAST-Pipeline-Source-ID-34260695.csv]], also obtainable from the above link.