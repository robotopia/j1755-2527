# Various processing notes

## Stacking spectra together

```
python stack_ds.py $(cat GPM_observations.txt) --dedisperse 0 --output_average_ds GPM_observations_stacked.pkl --output_average_ds_image GPM_stacked.png

python stack_ds.py $(cat D0042_observations_185MHz.txt) --dedisperse 0 --output_average_ds D0042_observations_185MHz_stacked.pkl --output_average_ds_image D0042_185MHz_stacked.png

python stack_ds.py $(cat D0042_observations_200MHz.txt) --dedisperse 0 --output_average_ds D0042_observations_200MHz_stacked.pkl --output_average_ds_image D0042_200MHz_stacked.png

python stack_ds.py $(cat 2023_meerkat_observations.txt) --dedisperse 0 --output_average_ds 2023_meerkat_stacked.pkl --output_average_ds_image 2023_meerkat_stacked.png

python stack_ds.py 1413381294_meerkat.pkl --dedisperse 0 --output_average_ds 2024_meerkat_stacked.pkl --output_average_ds_image 2024_meerkat_stacked.png

python plot_stacked_spectra.py
```

## Making a pulsestack

```
python pulsestack_all.py --ncols 2 --output_pulsestack_image pulsestack.pdf --xlim -300 300 1*pkl
```
