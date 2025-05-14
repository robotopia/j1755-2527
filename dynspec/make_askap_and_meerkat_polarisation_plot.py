import numpy as np
import matplotlib.pyplot as plt

import argparse

def main():
    parser = argparse.ArgumentParser(description="Make a plot showing the different polarisations of the original ASKAP pulse vs the more recent ASKAP and MeerKAT pulses")

    parser.add_argument('--output_plot', help="Save the output plot to the named file. If not supplied, default action is to plt.show()")

    args = parser.parse_args()

    # Hardcode which dynamic spectra we're opening for plotting
    pkls = ['1358297519_askap.pkl', '1404832334_askap.pkl', '1413381294_meerkat.pkl']

    for pkl in pkls:
        dat = np.load(dat, allow_pickle=True)

        # Get Stokes I and dedisperse
        t = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=EarthLocation.of_site(dat['TELESCOPE']))
        dt = t[1] - t[0]
        f = dat['FREQS'] * u.Hz
        try:
            df = f[1] - f[0]
        except:
            df = dat['BW'] * u.Hz
        nchans = len(f)

        f += 0.5*df # A hack because frequencies are labelled as lower edge of channels instead of middle of channels where they should be
        Ss = {}
        for pol in "IQUV":
            S, fmask, tmask = Stokes_ds(dat, pol=pol)

            # Dedispersing requires no nans
            S[np.isnan(I)] = 0.0

            # Dedisperse
            f_ref = np.nanmean(f)
            Sdd = dedisperse_ds(S, src_dm, f, f_ref, dt)

            # Reapply channel flags 
            if len(fmask) > 0:
                Sdd[:,fmask] = np.nan

            Ss[pol] = S





if __name__ == '__main__':
    main()
