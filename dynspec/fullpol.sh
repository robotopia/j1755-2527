#!/bin/bash

for stokes in I Q U V
do
  # One of the pulses straddles the boundary between observations 1410777264 and 1410777560, so we join them together like so:
  cat 1410777264-${stokes}.csv 1410777560-${stokes}.csv > 1410777264-${stokes}_and_1410777560-${stokes}.csv

  # Same for 1410773120 and 1410773416
  cat 1410773120-${stokes}.csv 1410773416-${stokes}.csv > 1410773120-${stokes}_and_1410773416-${stokes}.csv

done

# The corresponding YAML files have been made and are git-tracked:
# - 1410773120-I_and_1410773416-I.yaml
# - 1410777264-I_and_1410777560-I.yaml
# But we also need to broadcast these to the other polarisations
for stokes in Q U V
do
  sed "s/-I/-${stokes}/g" < 1410773120-I_and_1410773416-I.yaml > 1410773120-${stokes}_and_1410773416-${stokes}.yaml
  sed "s/-I/-${stokes}/g" < 1410777264-I_and_1410777560-I.yaml > 1410777264-${stokes}_and_1410777560-${stokes}.yaml
done
