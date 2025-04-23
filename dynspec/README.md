# Various processing notes

## Stacking spectra together

```
python stack_ds.py $(cat GPM_observations.txt | sed 's/$/.pkl/') --dedisperse 0 --output_average_ds GPM_observations_stacked.pkl --output_average_ds_image GPM_stacked.png

python stack_ds.py $(cat D0042_observations_185MHz.txt | sed 's/$/.pkl/') --dedisperse 0 --output_average_ds D0042_observations_185MHz_stacked.pkl --output_average_ds_image D0042_185MHz_stacked.png

python stack_ds.py $(cat D0042_observations_200MHz.txt | sed 's/$/.pkl/') --dedisperse 0 --output_average_ds D0042_observations_200MHz_stacked.pkl --output_average_ds_image D0042_200MHz_stacked.png
```
