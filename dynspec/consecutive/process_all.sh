#!/bin/bash

# Put consecutive observations together
python stitch_consecutive_obs_together.py

# The above script generated new YAML files, but the "Input file" field has to be updated
ls *_and_*.pkl | cut -d. -f1 | while read obsids
do
  sed -i "s/[0-9]*.pkl/${obsids}.pkl/" ${obsids}-I.yaml
done
