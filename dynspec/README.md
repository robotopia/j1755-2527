# Various processing notes

## Stacking spectra together

```
python stack_ds.py --bin_size 4 --output_average_ds_image GPM_stacked.png $(cat GPM_observations.txt | sed 's/$/.pkl/') --dedisperse 0 --fscrunch_factor 24

python stack_ds.py --bin_size 4 --output_average_ds_image D0042_185MHz_stacked.png $(cat D0042_observations_185MHz.txt | sed 's/$/.pkl/') --dedisperse 0 --fscrunch_factor 24

python stack_ds.py --bin_size 2 --output_average_ds_image D0042_200MHz_stacked.png $(cat D0042_observations_200MHz.txt | sed 's/$/.pkl/') --dedisperse 0 --fscrunch_factor 96
```
