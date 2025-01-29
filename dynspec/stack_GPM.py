import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys
import argparse

def fold(times, period, pepoch):

    # Dedisperse to infinite frequency
    #D = 4.148808e3 * u.MHz**2 / u.pc * u.cm**3 * u.s
    #times -= D * self.dm / freqs_MHz**2

    pulses_phases = ((times - pepoch) / period).decompose()
    pulses, phases = np.divmod(pulses_phases + 0.5, 1)
    phases -= 0.5

    return pulses, phases

def barycentre(times, src_coord):
    bc_correction = times.light_travel_time(src_coord, ephemeris='jpl')
    return times + bc_correction

def StokesI_ds(dat):
    if dat['POLS'][0] == 'XX' and dat['POLS'][3] == 'YY':
        return 0.5*np.real(dat['DS'][:,:,0] + dat['DS'][:,:,3])
    else:
        raise NotImplementedError(f"I haven't been taught how to deal with POLS = {ds['POLS']} yet")

def main():
    parser = argparse.ArgumentParser(description="Stack dynamic spectra together according to a folding ephemeris")

    parser.add_argument('coord', help="RA/Dec of source, in the format \"HH:MM:SS.S ±DD:MM:SS.S\"")
    parser.add_argument('period', type=float, help="The period of the source, in seconds")
    parser.add_argument('PEPOCH', type=float, help="A reference epoch (MJD) to mark zero rotation phase")
    parser.add_argument('DM', type=float, help="The dispersion measure, in pc/cm³")

    parser.add_argument('bin_size', type=float, help="The size of a time bin, in seconds")
    parser.add_argument('--fscrunch_factor', type=int, help="How many frequency channels to average together in the final output plot")

    parser.add_argument('output_pulsestack_image', help="The filename to use for the pulsestack image")
    parser.add_argument('output_average_ds_image', help="The filename to use for the average DS image")

    parser.add_argument('ds_files', nargs='*', help="The names of the input .pkl files, without the .pkl extension")

    args = parser.parse_args()

    src_coord = SkyCoord(args.coord, unit=(u.hourangle, u.deg), frame='icrs')
    src_period = args.period * u.s
    src_pepoch = Time(args.PEPOCH, scale='utc', format='mjd')
    src_dm = args.DM * u.pc / u.cm**3

    # Binning parameters
    bin_size = args.bin_size * u.s
    bin_size_phase = (bin_size / src_period).decompose()

    fscrunch_factor = args.fscrunch_factor

    # Output filenames
    pulsestack_png = args.output_pulsestack_image
    average_ds_png = args.output_average_ds_image

    min_phase_bin = np.inf
    max_phase_bin = -np.inf

    pkls = args.ds_files
    dats = []
    phase_bins_list = []

    # Load data, and convert the times to phases and phase bins
    for pkl in pkls:
        print(f"Opening {pkl}...")
        dat = np.load(pkl, allow_pickle=True)
        dat['OBSID'] = pkl[:-4].replace('_', ' ')

        location = EarthLocation.of_site(dat['TELESCOPE'])
        times = Time(dat['TIMES']/86400.0, scale='utc', format='mjd', location=location)
        barytimes = barycentre(times, src_coord)

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
    profile = np.nanmean(mean_ds, axis=1)

    # Draw a curve representing the DM
    f = dats[0]['FREQS']/1e6 # MHz
    dmdelay = 4.148808e3 * src_dm.to('pc cm-3').value * (1/f**2 - 1/np.mean(f)**2)

    # Do a bit of extra averaging in frequency
    mean_ds = np.nanmean(np.reshape(mean_ds, (mean_ds.shape[0], -1, fscrunch_factor)), axis=-1)
    counts = np.nanmean(np.reshape(counts, (counts.shape[0], -1, fscrunch_factor)), axis=-1)
    fscr = np.nanmean(np.reshape(f, (-1, fscrunch_factor)), axis=-1)

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ObsID")
    plt.tight_layout()
    plt.savefig(pulsestack_png)
    plt.close(fig)

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(5,8))
    t = bin_size.to('s').value * (np.arange(nphase_bins) + min_phase_bin)

    axs[0].plot(t, profile)
    axs[1].pcolormesh(t, fscr, mean_ds.T, vmin=-0.1, vmax=0.6)

    axs[1].plot(dmdelay - 10, f, 'r--', lw=2, label=f"DM = {src_dm.to('pc cm-3').value} pc/cm³")
    axs[1].plot(dmdelay + 80, f, 'r--', lw=2)

    cax = axs[2].pcolormesh(t, fscr, counts.T)
    cbar = fig.colorbar(cax, ax=axs[2], orientation='horizontal')
    cbar.set_label("Number of dynamic spectra\nthat contributed to each point")
    axs[-1].set_xlabel("Time (s)")
    axs[0].set_ylabel("Flux density (Jy/beam)")
    axs[1].set_ylabel("Frequency (MHz)")
    axs[1].legend()
    axs[2].set_ylabel("Frequency (MHz)")
    plt.tight_layout()
    plt.savefig(average_ds_png)

if __name__ == '__main__':
    main()
