It has already been reported in [Dobie et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024MNRAS.535..909D/abstract) that no detections were made in GPM 2022. That was before we had an ephemeris. Now we can check whether there were any dim pulses that occurred at the expected times. We will also need to redo GPM2024 because I need the dynamic spectra.

The script [[gpm/get_obsids.py]] does the cross-matching that produces a list of ObsIDs in GPM2024 in which J1755-2527 was (approximately) within 15 degrees of the pointing centre, and in which coincide with the pulse ToAs (actually, within ~1 min either side of the ToAs, because of the pulse width being about ~2 mins). The resulting list is stored in [[gpm/gpm2024_J1755-2527_observations.txt]].

When using this list with the GPM pipeline, make sure to strip the header, either manually, or with something like
```
sed -n '/^[^#]/ p' < gpm2024_J1755-2527_observations.txt > some_other_file.txt
```

Calibration solutions will have to be downloaded for these observations as well. You can use the pipeline's `gpm_track.py` utility to get the list of calibrators needed for these specific observations. It's a bit clunky (as of this writing), but it can done like this:
```
while read obsid
do
  singularity exec $GPMCONTAINER $GPMBASE/gpm_track.py obs_calibrator --obs_id $obsid
done < gpm2024_J1755-2527_observations.txt | sort -u | tee gpm2024_J1755-2527_calibrations.txt
```
It takes a minute or two to run, so I really need to make the `obs_calibrator` directive work with the `--obs_file` option. If I find I have to do it again lots, I'll make sure to implement this. I'll also save the output of this operation here: [[gpm/gpm2024_J1755-2527_calibrations.txt]].