import numpy as np
import pickle
import argparse
from astropy.io import fits
from timing import *

def main():

    parser = argparse.ArgumentParser(description="Calculate the primary beam correction and save the correction to the same pickle file (key = 'PB_CORR'). Currently only supports MWA observations. This will overwrite the 'PB_CORR' field in the pickle file.")
    parser.add_argument('ds_file', help="The name of the pickle file to read and write to.")
    parser.add_argument('metafits', help="The (MWA-style) metafits file associated with this observation. This is used to get the observations''GRIDNUM'.")
    parser.add_argument('--output_file', help="Write to this other file instead of overwriting the input file")

    args = parser.parse_args()

    dat = np.load(args.ds_file, allow_pickle=True)

    if dat['TELESCOPE'] != 'MWA':
        raise NotImplementedError(f"Primary beam correction calculation not implemented for {dat['TELESCOPE']}")

    # TODO: check that pols = ['XX', ...]

    with fits.open(args.metafits) as hdul:
        gridnum = hdul[0].header["GRIDNUM"]

    coord = get_J1755_ephemeris()['coord']
    MWA = EarthLocation.from_geodetic(lat=-26.703319 * u.deg, lon=116.67081 * u.deg, height=377 * u.m)
    freqs = dat['FREQS']
    t = Time(dat['TIMES'][0]/86400, scale='utc', format='mjd')

    rX, rY = np.array([beam_lookup_1d(coord.fk5.ra.deg, coord.fk5.dec.deg, gridnum, t, freq) for freq in freqs]).T

    dat['PB_CORR'] = np.array([rX, np.full(rX.shape, np.nan), np.full(rX.shape, np.nan), rY]).T

    with open(args.output_file or args.ds_file, 'wb') as pkl:
        pickle.dump(dat, pkl)


if __name__ == '__main__':
    main()
