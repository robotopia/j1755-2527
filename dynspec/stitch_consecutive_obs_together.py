#!/usr/bin/python

import numpy as np
from glob import glob
from copy import deepcopy
import shutil
import pickle

pklfiles = glob('??????????.pkl')
pklfiles.sort()

for i in range(len(pklfiles)-1):
    pklfile = pklfiles[i]
    pklfile_next = pklfiles[i+1]

    data = np.load(pklfile, allow_pickle=True)
    data_next = np.load(pklfile_next, allow_pickle=True)

    dt = data['TIMES'][1] - data['TIMES'][0]
    dt_next = data_next['TIMES'][1] - data_next['TIMES'][0]

    # Don't bother if the sample times are too different (differ by more than 0.01%)
    tol = 1e-4
    if np.abs(dt - dt_next) > dt*tol:
        #print(f"{np.abs(dt - dt_next) = }")
        continue

    # Don't bother if the time difference between the last sample of the first pickle
    # and the first sample of the second pickle is not close enough to the sample time
    if np.abs(data_next['TIMES'][0] - data['TIMES'][-1] - dt) > dt*tol:
        #print(f"{np.abs(data_next['TIMES'][0] - data['TIMES'][-1]) = }")
        continue

    # Don't bother if the first frequency isn't the same
    if np.abs(data_next['FREQS'][0] - data['FREQS'][0]) > data['FREQS'][0]*tol:
        continue

    # If we got this far, then *this* obs (data) and the *next* obs (data_next) are
    # consecutive and should be packed together.
    new_pklfile = f"{pklfile[:-4]}_and_{pklfile_next}"
    print(f"Packing {pklfile} and {pklfile_next} into {new_pklfile}")
    new_data = deepcopy(data)

    new_data['TIMES'] = np.hstack((data['TIMES'], data_next['TIMES']))
    new_data['DS'] = np.vstack((data['DS'], data_next['DS']))
    #print(f"{new_data['DS'].shape = }")

    with open(new_pklfile, 'wb') as f:
        pickle.dump(new_data, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Generate a new YAML to go along with it
    for stokes in ['I']: # Can expand to other Stokes later, if needed
        old_yaml = f"{pklfile[:-4]}-{stokes}.yaml"
        new_yaml = f"{new_pklfile[:-4]}-{stokes}.yaml"
        shutil.copyfile(old_yaml, new_yaml)

        # The contents of the YAML will have to be changed slightly, but this is
        # far easier to do with sed than here in Python.

