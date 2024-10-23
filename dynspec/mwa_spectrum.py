import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, correlate2d, correlation_lags
from scipy.optimize import curve_fit

# Load the raw dynamic spectrum of interest
ds = np.loadtxt("1410773120-I_and_1410773416-I.csv").T
f0 = 169.6 # Centre of lowest channel (MHz)
df = 0.04 # Channel width (MHz)
bw = 30.72
flo = f0 - df/2
fhi = flo + bw
fctr = (flo + fhi)/2
dt = 4.0 # Sample time (s)

print(f"original dynamic spectrum shape: {ds.shape}")

# Get the normalised lightcurve
lc = np.nanmean(ds, axis=0)

# Smooth it a bit
lc_smoothed = savgol_filter(lc, 15, 3) # window size 51, polynomial order 3
lc_smoothed /= np.nanmean(lc_smoothed)

# Avg the data into "coarse channels" (cc)
nchan, nbins = ds.shape
ds_cc = np.nanmean(ds.reshape(24, -1, nbins), axis=1)
cc_nchans = ds_cc.shape[0]
cc = np.arange(cc_nchans)
print(f"averaging coarse channels. New shape: {ds_cc.shape}")

# Cross correlate with smoothed lightcurve
correlated = correlate2d(ds_cc, [lc_smoothed], mode='same') / len(lc_smoothed)
lags = correlation_lags(ds_cc.shape[-1], len(lc_smoothed), mode='same')

# Find the peaks nearest to lag=0
zero_lag_idx = np.where(lags==0)[0][0]
local_peaks = np.logical_and(correlated[:,1:-1] > correlated[:,2:],
                             correlated[:,1:-1] > correlated[:,:-2])
local_peak_idxs = [np.where(local_peaks[c,:])[0] for c in cc]
local_peak_distances = [np.abs(local_peak_idxs[c] - zero_lag_idx).squeeze() for c in cc]
nearest_local_peak_idxs_idx = np.array([np.argmin(local_peak_distances[c]) for c in cc])
#print(local_peak_distances)
#print(nearest_local_peak_idxs_idx)

# Convert "nearest_local_peak_idxs_idx" (which are idxs into the shortlists of idxs)
# into the idxs that can refer to the original ds_cc

#print(local_peak_idxs)
nearest_local_peak_idxs = np.array([local_peak_idxs[c][nearest_local_peak_idxs_idx[c]] for c in cc])
#print(nearest_local_peak_idxs)

# Fit a curve to the positions to get a nominal DM
cc_chan_width = bw / cc_nchans # Width of coarse channels
cc_f0 = flo + cc_chan_width/2 # Centre of lowest coarse channel
cc_freqs = np.arange(cc_nchans)*cc_chan_width + cc_f0 # List of frequencies
times = nearest_local_peak_idxs*dt

def dmdelay(freq, t0, dm):
    return t0 + 4.148808e3*dm/freq**2

#print(local_peak_idxs)
#print("------------")
#print(local_peak_distances)

# time at infinite frequency, from start of obs
#       |    DM
#       |    |
p0 = (50*dt, 150)
popt, pcov = curve_fit(dmdelay, cc_freqs, times, p0=p0)
#print(f"{popt = }")
#print(f"{np.sqrt(pcov) = }")

#plt.plot(lc)
#plt.plot(lc_smoothed)
plt.pcolormesh((lags + zero_lag_idx)*dt, cc_freqs, ds_cc)
plt.plot(nearest_local_peak_idxs*dt, cc_freqs, 'rx', label="Correlation peaks")
plt.plot(dmdelay(cc_freqs, *popt), cc_freqs, 'w--', label=f"DM = {popt[1]:.0f} ± {np.sqrt(pcov[1,1]):.0f} pc/cm^3")
plt.xlabel("Time from start of obs (s)")
plt.ylabel("Frequency (MHz)")
plt.legend()
plt.savefig("1410773120-I_and_1410773416-I_dm_fit.png")
plt.clf()


# Get the predicted pulse centres according to the DM fit
dm_times = dmdelay(cc_freqs, *popt)

# Convert these to time bin indexes
dm_idxs = np.round(dm_times/dt).astype(int)

# Get the correlation value at those points as a measurement of brightness of pulse
# in each channel. The correlation kernel (i.e. the smoothed lightcurve) was normalised
# so that the correlation values are still in Jy
cc_fluxes = np.array([correlated[c,dm_idxs[c]] for c in cc])

# Throw away a few outliers
cc_fluxes[9:12] = np.nan
not_nan_mask = ~np.isnan(cc_fluxes)
cc_goodfluxes = cc_fluxes[not_nan_mask]
cc_goodfreqs = cc_freqs[not_nan_mask]

# Fit a spectral index, coz that's what we like to do
def powerlaw(nu_GHz, S1GHz, alpha):
    return S1GHz*nu_GHz**alpha

p0 = (0.05, -1)
popt, pcov = curve_fit(powerlaw, cc_goodfreqs/1e3, cc_goodfluxes, p0=p0)
print(popt, pcov)

freqs = np.logspace(2, 3, 10) # Used for plotting

plt.plot(cc_freqs, cc_fluxes, 'k.', label="MWA")
plt.plot([887.5], [0.2], 'bx', label="ASKAP") # Add Dougal's point, and draw their slope as well
#plt.plot(freqs, powerlaw(freqs/1e3, *popt), 'k--', alpha=0.5, label=f"$\\alpha = {popt[1]:.1f} ± {np.sqrt(pcov[1,1]):.1f}$")


plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("1410773120-I_and_1410773416-I_spectrum.png")
