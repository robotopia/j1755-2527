#!/bin/bash

# One of the pulses straddles the boundary between observations 1410777264 and 1410777560, so we join them together like so:
cat 1410777264-I.csv 1410777560-I.csv > 1410777264-I_and_1410777560-I.csv

# Same for 1410773120 and 1410773416
cat 1410773120-I.csv 1410773416-I.csv > 1410773120-I_and_1410773416-I.csv

# The corresponding YAML files have been made and are git-tracked:
# - 1410773120-I_and_1410773416-I.yaml
# - 1410777264-I_and_1410777560-I.yaml

