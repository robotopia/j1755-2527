import numpy as np
import pickle
from astropy.time import Time

##########
# First light curve (from 2023)

I = np.loadtxt('askap1-I.txt')

# Save to pickle file
ds = np.empty((len(I), 1, 4))
ds[:,0,0] = I # Stokes I
ds[:,:,1:] = np.nan # Stokes Q, U, V
dat = {
    'DS': ds,
    'DS_MED': ds,
    'DS_STD': ds,
    'POLS': ['I', 'Q', 'U', 'V'],
    'TIMES': np.arange(len(I))*10 + 59965.03588541667*86400,
    'FREQS': np.array([887500000]),
    'TELESCOPE': 'ASKAP',
}
pkl = f"../dynspec/{Time(59965.03588541667, scale='utc', format='mjd').gps:.0f}_askap.pkl"
with open(pkl, 'wb') as f:
    pickle.dump(dat, f)

##########
# Second light curve (from 2024)

ts = []
Is = []
Qs = []
Us = []
Vs = []
with open('j1755_burst2_unsub_lc.csv', 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:
        t, I, L, Q, U, V, pa = line.split(',')
        ts.append(float(t))
        Is.append(float(I))
        Qs.append(float(Q))
        Us.append(float(U))
        Vs.append(float(V))

# Save to pickle file
ds = np.empty((len(ts), 1, 4))
ds[:,0,0] = Is
ds[:,0,1] = Qs
ds[:,0,2] = Us
ds[:,0,3] = Vs
dat = {
    'DS': ds,
    'DS_MED': ds,
    'DS_STD': ds,
    'POLS': ['I', 'Q', 'U', 'V'],
    'TIMES': np.array(ts),
    'FREQS': np.array([887500000]),
    'TELESCOPE': 'ASKAP',
}
pkl = f"../dynspec/{Time(ts[0]/86400, scale='utc', format='mjd').gps:.0f}_askap.pkl"
with open(pkl, 'wb') as f:
    pickle.dump(dat, f)

