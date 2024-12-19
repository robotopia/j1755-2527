import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys

src_coord = SkyCoord("17:55:34.87", "-25:27:49.1", unit=(u.hourangle, u.deg), frame='icrs')
src_period = 4186.330360402573 * u.s
src_pepoch = Time(59965.038453386056, scale='utc', format='mjd')
src_dm = 922.4853409645726 * u.pc / u.cm**3

# Binning parameters
bin_size = 4 * u.s
bin_size_phase = (bin_size / src_period).decompose()

def fold(times, period, pepoch):

    # Dedisperse to infinite frequency
    #D = 4.148808e3 * u.MHz**2 / u.pc * u.cm**3 * u.s
    #times -= D * self.dm / freqs_MHz**2

    pulses_phases = ((times - pepoch) / period).decompose()
    pulses, phases = np.divmod(pulses_phases + 0.5, 1)
    phases -= 0.5

    return pulses, phases

def barycentre(times):
    bc_correction = times.light_travel_time(src_coord, ephemeris='jpl')
    return times + bc_correction

def StokesI_ds(dat):
    if dat['POLS'][0] == 'XX' and dat['POLS'][3] == 'YY':
        return 0.5*np.real(dat['DS'][:,:,0] + dat['DS'][:,:,3])
    else:
        raise NotImplementedError(f"I haven't been taught how to deal with POLS = {ds['POLS']} yet")

def main():
    mins = []
    maxs = []

    # Work out the dimensions of the final output array
    # This has all been prepared in advance, painstakingly.
    # the minimum phase bin (using 4-s bins) is index = -78,
    # and the max is index = 97
    min_phase_bin = -78
    max_phase_bin = 97
    nphase_bins = max_phase_bin - min_phase_bin + 1
    nchans = 768

    output_ds = np.zeros((nphase_bins, nchans))
    counts = np.zeros((nphase_bins, nchans)) # For tracking how many items contribute to the same bin

    pkls = sys.argv[1:]
    for pkl in pkls:
        print(f"Opening {pkl}...")
        dat = np.load(pkl, allow_pickle=True)

        if int(dat['FREQS'][0]) != 184960000 or int(dat['FREQS'][1]) != 185000000:
            print(f"{pkl} ({int(dat['FREQS'][0]), int(dat['FREQS'][1])}...) is not in the desired frequency range (184960000...). Skipping...")
            continue

        I = StokesI_ds(dat)
        location = EarthLocation.of_site(dat['TELESCOPE'])
        times = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=location)
        barytimes = barycentre(times)

        _, phases = fold(barytimes, src_period, src_pepoch)
        phase_bins = np.round(phases / bin_size_phase).astype(int)

        mins.append(min(phase_bins))
        maxs.append(max(phase_bins))

        # All of the dynamic spectra I'm considering in this run have a sample time of 4 seconds, which is also what the bin size has been set to. This means that each output bin will only ever have one pixel from each dynamic spectrum added to it. If this assumption breaks, then the following method of adding each dynamic spectrum into the output dynamic spectrum won't work.
        output_ds[(phase_bins[0] - min_phase_bin):(phase_bins[0] - min_phase_bin + I.shape[0])] = np.nansum([output_ds[(phase_bins[0] - min_phase_bin):(phase_bins[0] - min_phase_bin + I.shape[0])], I], axis=0)
        counts[(phase_bins[0] - min_phase_bin):(phase_bins[0] - min_phase_bin + nphase_bins)] += ~np.isnan(counts[(phase_bins[0] - min_phase_bin):(phase_bins[0] - min_phase_bin + nphase_bins)])

        #plt.pcolormesh(I.T)
        #plt.show()
        #print(phase_bins)

    mean_ds = output_ds / counts

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
    axs[0].plot(np.nanmean(mean_ds, axis=1))
    axs[1].pcolormesh(mean_ds.T, vmin=-0.1, vmax=0.6)
    axs[2].pcolormesh(counts.T)
    plt.show()

if __name__ == '__main__':
    main()
