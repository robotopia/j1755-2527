import numpy as np
import pickle
from astropy.time import Time

pkl = '1685306788_sdp_l0_1024ch_J1755-2527_min_bl150.pkl'

dat = np.load(pkl, allow_pickle=True)

time_ranges = [
    [5192026003, 5192026601],
    [5192030691, 5192030890],
    [5192034613, 5192035212],
    [5192038693, 5192039292],
]

output_dicts = [{} for time_range in time_ranges]

# Pull out the bits to be extracted, and save them temporarily in memory
for i in range(len(time_ranges)):
    time_range = time_ranges[i]
    mask = (dat['TIMES'] >= time_range[0]) & (dat['TIMES'] <= time_range[1])

    output_dicts[i]['TIMES']  = dat['TIMES'][mask]
    output_dicts[i]['DS']     = dat['DS'][mask,:,:]
    output_dicts[i]['DS_STD'] = dat['DS_STD'][mask,:,:]
    output_dicts[i]['DS_MED'] = dat['DS_MED'][mask,:,:]

# Use the existing dat struct as a template for writing out the (four) new, smaller pickle files
for i in range(len(time_ranges)):
    output_dict = output_dicts[i]

    dat['TIMES']  = output_dict['TIMES']
    dat['DS']     = output_dict['DS']
    dat['DS_STD'] = output_dict['DS_STD']
    dat['DS_MED'] = output_dict['DS_MED']

    start_gps = Time(dat['TIMES'][0]/86400, scale='utc', format='mjd').gps
    output_pkl = f'{int(np.round(start_gps))}_meerkat.pkl'
    with open(output_pkl, 'wb') as f:
        pickle.dump(dat, f)
