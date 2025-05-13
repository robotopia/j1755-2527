import numpy as np
import pickle
import argparse
import os
from mwa_pb_lookup.lookup_beam import beam_lookup_1d

from astropy.io import fits
import astropy.units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord
from timing import *

MeerKAT_D = 13.5 * u.m

def gaussian_1d(x, mu, sigma):
    """
    Calculates the Gaussian function for each element in x.

    Parameters:
      x (array_like): The input array.
      mu (float): The mean of the Gaussian distribution.
      sigma (float): The standard deviation of the Gaussian distribution.

    Returns:
      ndarray: An array of the same shape as x, containing the Gaussian values.
    """
    #return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma)**2)
    return np.exp(-0.5 * ((x - mu) / sigma)**2)

def main():

    parser = argparse.ArgumentParser(description="Calculate the primary beam correction and save the correction to the same pickle file (key = 'PB_CORR'). Currently only supports MWA observations. This will overwrite the 'PB_CORR' field in the pickle file.")
    parser.add_argument('ds_file', help="The name of the pickle file to read and write to.")
    parser.add_argument('--metafits', help="The (MWA-style) metafits file associated with this observation. This is used to get the observations''GRIDNUM'. Required for MWA calculations.")
    parser.add_argument('--output_file', help="Write to this other file instead of overwriting the input file.")
    parser.add_argument('--overwrite', action='store_true', help="Write to this other file instead of overwriting the input file")
    parser.add_argument('--pointing', help="\"HH:MM:SS ±DD:MM:SS\" of telescope pointing direction. Needed for MeerKAT PB calculation")

    args = parser.parse_args()

    output_file = args.output_file or args.ds_file

    # Only proceed if the output file doesn't already have a PB_CORR column in it, OR if overwrite flag is set
    if os.path.exists(output_file):
        try:
            dat = np.load(output_file, allow_pickle=True)
        except:
            response = input(f"Unable to open {output_file} as a pickled dynamic spectrum. Continuing will overwrite this file completely. Do you want to continue? [y/N]: ").strip().lower()
            if response != 'y':
                print("Aborting.")
                exit()
        if 'PB_CORR' in dat.keys():
            if args.overwrite == False:
                print(f"{output_file} already contains a PB_CORR column. Skipping.")
                exit()

    # Read the input file. If the output file *is* the input file, then we've already read it in, above
    if output_file != args.ds_file:
        dat = np.load(args.ds_file, allow_pickle=True)

    coord = get_J1755_ephemeris()['coord']

    if dat['TELESCOPE'] == 'MWA':

        # TODO: check that pols = ['XX', ...]

        if not args.metafits:
            raise ValueError("--metafits option required for MWA PB calculation")

        with fits.open(args.metafits) as hdul:
            gridnum = hdul[0].header["GRIDNUM"]

        MWA = EarthLocation.from_geodetic(lat=-26.703319 * u.deg, lon=116.67081 * u.deg, height=377 * u.m)
        freqs = dat['FREQS']
        t = Time(dat['TIMES'][0]/86400, scale='utc', format='mjd')

        rX, rY = np.array([beam_lookup_1d(coord.fk5.ra.deg, coord.fk5.dec.deg, gridnum, t, freq) for freq in freqs]).T

        dat['PB_CORR'] = np.array([rX, np.full(rX.shape, np.nan), np.full(rX.shape, np.nan), rY]).T

    elif dat['TELESCOPE'] == 'MeerKAT':
        if not args.pointing:
            raise ValueError("--pointing option required for MeerKAT PB calculation")

        try:
            pointing = SkyCoord(args.pointing, unit=(u.hour, u.deg), frame='fk5')
        except:
            raise ValueError(f"Could not parse '{args.pointing}' as HH:MM:SS ±DD:MM:SS")

        freqs = dat['FREQS'] * u.Hz

        # Calibrated to the plots of Villiers 2023
        fwhms = 1.13*(((c/freqs) / MeerKAT_D)*u.rad).to('deg')
        sigmas = fwhms/2.355
        sep = pointing.separation(coord)

        dat['PB_CORR'] = np.broadcast_to([gaussian_1d(sep.deg, 0, sigma.value) for sigma in sigmas], (4, len(freqs))).T

    else:
        raise NotImplementedError(f"Primary beam correction calculation not implemented for {dat['TELESCOPE']}")

    with open(output_file, 'wb') as pkl:
        pickle.dump(dat, pkl)


if __name__ == '__main__':
    main()
