import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import sys
import argparse

def fold(times, period, pepoch):

    pulses_phases = ((times - pepoch) / period).decompose()
    pulses, phases = np.divmod(pulses_phases + 0.5, 1)
    phases -= 0.5

    return pulses, phases

def barycentre(times, src_coord):
    bc_correction = times.light_travel_time(src_coord, ephemeris='jpl')
    return times + bc_correction

def calc_dmdelay(dm, f, f_ref):
    return 4.148808e3 * dm * (1/f**2 - 1/f_ref**2)

def StokesI_ds(dat):
    if dat['POLS'][0] == 'XX' and dat['POLS'][3] == 'YY':
        return 0.5*np.real(dat['DS'][:,:,0] + dat['DS'][:,:,3])
    else:
        raise NotImplementedError(f"I haven't been taught how to deal with POLS = {ds['POLS']} yet")

def dedisperse_ds(ds, dm, f, f_ref, bin_size_seconds):
    dmdelays = calc_dmdelay(dm, f, f_ref)

    # It doesn't matter if the edges are wrapped, so use simple FFT to do
    # the DM delay adjustment (essentially, a "roll" with fractional bins)
    FFT = np.fft.rfft(ds, axis=0)
    FFT_freqs = np.fft.rfftfreq(ds.shape[0], d=bin_size_seconds)

    # Make a mesh grid of delays and frequencies with the same size as the FFT matrix
    DMDELAYS, FFT_FREQS = np.meshgrid(dmdelays, FFT_freqs)

    # Apply the phase ramp
    FFT *= np.exp(2j*np.pi*FFT_FREQS*DMDELAYS)

    # Invert the FFT
    return np.fft.irfft(FFT, n=ds.shape[0], axis=0)

def main():
    parser = argparse.ArgumentParser(description="Stack dynamic spectra together according to a folding ephemeris")

    parser.add_argument('coord', help="RA/Dec of source, in the format \"HH:MM:SS.S ±DD:MM:SS.S\"")
    parser.add_argument('period', type=float, help="The period of the source, in seconds")
    parser.add_argument('PEPOCH', type=float, help="A reference epoch (MJD) to mark zero rotation phase")
    parser.add_argument('--DM', type=float, help="The dispersion measure, in pc/cm³")

    parser.add_argument('bin_size', type=float, help="The size of a time bin, in seconds")
    parser.add_argument('--fscrunch_factor', type=int, help="How many frequency channels to average together in the final output plot")
    parser.add_argument('--draw_dm_sweep', type=float, nargs='*', help="Time(s) (at centre frequency) at which to draw DM sweep curves in the synamic spectra panel")

    parser.add_argument('output_pulsestack_image', help="The filename to use for the pulsestack image")
    parser.add_argument('output_average_ds_image', help="The filename to use for the average DS image")

    parser.add_argument('ds_files', nargs='*', help="The names of the input .pkl files, without the .pkl extension")

    args = parser.parse_args()

    src_coord = SkyCoord(args.coord, unit=(u.hourangle, u.deg), frame='icrs')
    src_period = args.period * u.s
    src_pepoch = Time(args.PEPOCH, scale='utc', format='mjd')
    if args.DM is not None:
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

    f = dats[0]['FREQS']/1e6 # MHz
    
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
    t = bin_size.to('s').value * (np.arange(nphase_bins) + min_phase_bin)

    # If a DM has been provided, dedisperse the composite dynamic spectrum
    if args.DM is not None:
        f_ref = np.nanmean(f)
        mean_ds = dedisperse_ds(mean_ds, args.DM, f, f_ref, bin_size.to('s').value)
        counts = dedisperse_ds(counts, args.DM, f, f_ref, bin_size.to('s').value)
        t -= calc_dmdelay(args.DM, f_ref, np.inf)

    # Form the profile
    profile = np.nanmean(mean_ds, axis=1)

    # Do a bit of extra averaging in frequency
    if fscrunch_factor is not None:
        mean_ds = np.nanmean(np.reshape(mean_ds, (mean_ds.shape[0], -1, fscrunch_factor)), axis=-1)
        counts = np.nanmean(np.reshape(counts, (counts.shape[0], -1, fscrunch_factor)), axis=-1)
        fscr = np.nanmean(np.reshape(f, (-1, fscrunch_factor)), axis=-1)
    else:
        fscr = f

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ObsID")
    plt.tight_layout()
    plt.savefig(pulsestack_png)
    plt.close(fig)

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(5,8))

    axs[0].plot(t, profile)
    axs[1].pcolormesh(t, fscr, mean_ds.T, vmin=-0.1, vmax=0.6)

    # Draw a curve representing the DM
    if args.DM is not None and args.draw_dm_sweep is not None:
        first = True
        dmdelay = -calc_dmdelay(args.DM, f, np.nanmean(f))
        # The minus sign ^^^ is because the dynamic spectrum will already be dedispersed
        for time in args.draw_dm_sweep:
            label = f"Zero DM curve" if first else None
            axs[1].plot(dmdelay + time, f, 'w', lw=2, label=label)
            first = False

    cax = axs[2].pcolormesh(t, fscr, counts.T)
    cbar = fig.colorbar(cax, ax=axs[2], orientation='horizontal')
    cbar.set_label("Number of dynamic spectra\nthat contribute to each point")
    axs[-1].set_xlabel("Time (s)")
    axs[0].set_ylabel("Flux density (Jy/beam)")
    axs[1].set_ylabel("Frequency (MHz)")
    if args.DM is not None and args.draw_dm_sweep is not None:
        axs[1].legend()
    if args.DM is not None:
        axs[0].set_title(f"Dedispersed to DM = {args.DM} pc/cm³")
    axs[2].set_ylabel("Frequency (MHz)")

    plt.tight_layout()

    # Change the shape of the "counts" panel to make the x-label ("Time (s)")
    # not be obscured by the color bar
    pos = axs[2].get_position() # [x0, y0, width, height]
    delta = 0.02 # inches?
    pos.y0 += delta
    axs[2].set_position(pos)

    plt.savefig(average_ds_png)

if __name__ == '__main__':
    main()
