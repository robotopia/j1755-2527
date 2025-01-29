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

# Output filenames
pulsestack_png = "pulsestack.png"
average_ds_png = "average_ds.png"

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
    min_phase_bin = np.inf
    max_phase_bin = -np.inf

    pkls = sys.argv[1:]
    dats = []
    phase_bins_list = []

    # Load data, and convert the times to phases and phase bins
    for pkl in pkls:
        print(f"Opening {pkl}...")
        dat = np.load(pkl, allow_pickle=True)
        dat['OBSID'] = pkl[:-4].replace('_', ' ')

        location = EarthLocation.of_site(dat['TELESCOPE'])
        times = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=location)
        barytimes = barycentre(times)

        _, phases = fold(barytimes, src_period, src_pepoch)
        phase_bins = np.round(phases / bin_size_phase).astype(int)
        phase_bins_list.append(phase_bins)

        if min(phase_bins) < min_phase_bin:
            min_phase_bin = int(np.round(min(phase_bins)))
        if max(phase_bins) > max_phase_bin:
            max_phase_bin = int(np.round(max(phase_bins)))

        dat['PHASE_BINS'] = phase_bins
        dats.append(dat)

    # Work out the dimensions of the final output array
    nphase_bins = max_phase_bin - min_phase_bin + 1
    nchans = len(dats[0]['FREQS'])

    output_ds = np.zeros((nphase_bins, nchans))
    counts = np.zeros((nphase_bins, nchans)) # For tracking how many items contribute to the same bin

    # Go through the dynamic spectra and add them to the final output array. Also make a pulsestack
    fig, ax = plt.subplots(figsize=(8,12))
    i = 0
    yticks = []
    ylabels = []
    for dat in dats:

        I = StokesI_ds(dat)
        phase_bins = dat['PHASE_BINS']

        ax.plot(phase_bins*bin_size, np.nanmean(I, axis=1)+i, 'k', alpha=0.5)
        yticks.append(i)
        ylabels.append(dat['OBSID'])
        i += 1

        # All of the dynamic spectra I'm considering in this run have the same sample time, which is also what the bin size has been set to. This means that each output bin will only ever have one pixel from each dynamic spectrum added to it. If this assumption breaks, then the following method of adding each dynamic spectrum into the output dynamic spectrum won't work.
        start = phase_bins[0] - min_phase_bin
        end = start + I.shape[0]
        output_ds[start:end] = np.nansum([output_ds[start:end], I], axis=0)
        counts[start:end] += ~np.isnan(counts[start:end])

        #plt.pcolormesh(I.T)
        #plt.show()
        #print(phase_bins)

    mean_ds = output_ds / counts

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ObsID")
    plt.tight_layout()
    plt.savefig(pulsestack_png)
    plt.close(fig)

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8,12))
    axs[0].plot(np.nanmean(mean_ds, axis=1))
    axs[1].pcolormesh(mean_ds.T, vmin=-0.1, vmax=0.6)
    axs[2].pcolormesh(counts.T)
    plt.tight_layout()
    plt.savefig(average_ds_png)

if __name__ == '__main__':
    main()
