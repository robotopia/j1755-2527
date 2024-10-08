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

# Now just run the polprof.py script on all the obs with pulses
python polprof.py --rel_time --off_pulse_lims 0:50,200:300 --DM 1100 --RM 961 --outdynspec 1410785848-dynspec.png --outfile 1410785848-fullpol.png 1410785848-?.yaml

python polprof.py --rel_time --DM 1100 --RM 961 --outdynspec 1410773120_and_1410773416-dynspec.png --outfile 1410773120_and_1410773416-fullpol.png 1410773120-?_and_1410773416-?.yaml
