# Data processing notes

### Pull the data from Acacia

In the parent directory to this (`/dynspec`):

```
rclone copy mwasci:askapj1755/J1755-25_1729341386_scan18_images_models.tar.gz .
tar -xzf J1755-25_1729341386_scan18_images_models.tar.gz \
  J1755-25_1729341386_scan18-I.csv \
  J1755-25_1729341386_scan18-Q.csv \
  J1755-25_1729341386_scan18-U.csv \
  J1755-25_1729341386_scan18-V.csv \
  J1755-25_1729341386_scan18-I.yaml \
  J1755-25_1729341386_scan18-Q.yaml \
  J1755-25_1729341386_scan18-U.yaml \
  J1755-25_1729341386_scan18-V.yaml
```

### Rename the files to be consistent with the naming convention

Also in the parent directory to this (`/dynspec`):

```
mv J1755-25_1729341386_scan18-I.yaml 1413381294-I.yaml
mv J1755-25_1729341386_scan18-Q.yaml 1413381294-Q.yaml
mv J1755-25_1729341386_scan18-U.yaml 1413381294-U.yaml
mv J1755-25_1729341386_scan18-V.yaml 1413381294-V.yaml
mv J1755-25_1729341386_scan18-I.csv 1413381294-I.csv
mv J1755-25_1729341386_scan18-Q.csv 1413381294-Q.csv
mv J1755-25_1729341386_scan18-U.csv 1413381294-U.csv
mv J1755-25_1729341386_scan18-V.csv 1413381294-V.csv
```

I had to manually edit the YAML file to fix some typos, update it with the new file names, correct the "transpose" option for this data set, and also flag one of the time bins at the end. But once you make the changes, say, to the Stokes I YAML file, the others can be immediately changed with

```
for s in Q U V
do
  sed "s/-I/-$s/" < 1413381294-I.yaml > 1413381294-${s}.yaml
done
```

## Checking the zero-DM lightcurve

Uses https://github.com/GLEAM-X/dynspec_tools
```
dedisperse_dynspec --yaml 1413381294-I.yaml --dms 0 \
  --lightcurve 1413381294/1413381294-I_DM0_lightcurve.png \
  --dynspec_image 1413381294/1413381294-I_DM0.png
```

![[dynspec/1413381294/1413381294-I_DM0_fixed.png]]

Looking ok, except that the baseline is off. (Also, I believe the flux density units should be in Jy.) Interestingly, there *is* low-level, but apparently real emission in the original detected pulse in the few hundred sections prior to the 2-minute-long bright burst, so perhaps the true baseline is where it levels off on the right hand side of (i.e. after) the bright pulse above. 

Let's test that by looking at the polarisation. First, the dynamic spectra.

```
for S in Q U V
do
  dedisperse_dynspec --yaml 1413381294-${S}.yaml --dms 0 \
    --lightcurve 1413381294/1413381294-${S}_DM0_lightcurve.png \
    --dynspec_image 1413381294/1413381294-${S}_DM0.png
done
```

And then to save the results to this repo:

```
for S in Q U V
do
  cp 1413381294/1413381294-${S}_DM0.png 1413381294/1413381294-${S}_DM0_fixed.png
done
```

![[dynspec/1413381294/1413381294-Q_DM0_fixed.png]]

![[dynspec/1413381294/1413381294-U_DM0_fixed.png]]

![[dynspec/1413381294/1413381294-V_DM0_fixed.png]]

So they look ok in the sense that their baselines look sensible, as expected. Doing baseline subtraction for Stokes I and computing the polarisation profile:

> [!warning]
> Requires `>=v0.3.29` of [dynspec_tools](https://github.com/GLEAM-X/dynspec_tools).

```
polprof --RM -961 \
  --off_pulse_lims 1413381750:1413381890 \
  --outfile 1413381294/1413381294-DM0-polprof.png \
  --subtract_I_baseline \
  --PA_sigma_limit 1.5 \
  1413381294-[IQUV].yaml
```

> [!warning]
> I used **-961** (i.e. *negative*) instead of the positive value reported by Dougal. Might be my sign convention is wrong in my code.

![[dynspec/1413381294/1413381294-DM0-polprof_fixed.png]]

So, I can't see much polarisation in the first couple hundred seconds, but hey, look at that PA curve! There's an open question in my mind whether the RM could be changing across the pulse(s), and whether that might be shifting the respective PAs up and down. Let me look at the spectra.

```
for RM in $(seq 940 980)
do
  polprof --RM -${RM} \
    --off_pulse_lims 1413381750:1413381890 \
    --outfile 1413381294/1413381294-DM0-RM${RM}-polprof.png \
    --subtract_I_baseline \
    --PA_sigma_limit 1.5 \
    --outdynspec_binning 8 64 \
    --outdynspec 1413381294/1413381294-DM0-RM${RM}-dynspec.png \
    1413381294-[IQUV].yaml
done
magick 1413381294/1413381294-DM0-RM*-dynspec.png \
  -delay 5 -loop 0 1413381294/1413381294-DM0-RMs-dynspec.gif
magick 1413381294/1413381294-DM0-RM*-polprof.png \
  -delay 5 -loop 0 1413381294/1413381294-DM0-RMs-polprof.gif
```

![[1413381294-DM0-RMs-dynspec.gif]]
![[1413381294-DM0-RMs-polprof.gif]]

Actually, now that I look at this animation more closely, I'm not so sure that the PA doesn't change shape. At the beginning of the animation, the right hand component is more negatively sloped, while at the end of the animation, it's more positively sloped. That means the trailing side of the component changes its PA more quickly as a function of RM correction than the right hand side.

What would do that? If I imagine the PA as a function of $\lambda^2$, it's hard to imagine how a deviation from a straight line would have this effect. But I can think of one idea: the final reported PA is the frequency-scrunched version, so if the lower frequencies were weighted differently compared to the higher frequencies across the pulse, then the PAs will be also weighted differently. So I reckon what we're seeing is evidence of a different spectral index (at least, in the linear polarisation) across the pulse.

Let's now test this by measuring the spectral index in different bins. But I'll do it for both stokes I and L.