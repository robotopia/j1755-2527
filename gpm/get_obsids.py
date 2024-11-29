import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import sys

# Define source location
src_coord = SkyCoord("17h55m34.9s -25d27m49.1s", frame='icrs')
print(f"Source coordinates: {src_coord}")

# Get MWA obsids in GPM within 15 degrees of source (very coarse, will be refined later)
obs = []
page = 1
while True:
    print(f"Attempting to retrieve page {page} from https://ws.mwatelescope.org/")
    response = requests.get(f"https://ws.mwatelescope.org/metadata/find/?pagesize=500&notdeleted=on&projectid=G0080&minra=234&maxra=304&mindec=-40&maxdec=-10&page={page}")
    results = response.json()
    if len(results) == 0:
        print("Empty page. Stopping loop")
        break
    page += 1
    obs += results

print(f"Total of {len(obs)} observations retrieved")
obsids, _, _, _, RAs, Decs = list(map(list, zip(*obs)))

# Convert obsids to Time object
obsids = Time(obsids, scale='utc', format='gps')

pointings = SkyCoord(RAs, Decs, unit=(u.deg, u.deg), frame='icrs')
separations = src_coord.separation(pointings)

# Now read in the predicted ToAs from the ULP website
predicted_ToAs = Time(np.loadtxt('predicted_toas_gpm2024.txt'), scale='utc', format='gps')

# The pulses are about 2 minutes wide, so we want to keep any MWA obs which includes
# the timestamp ~1 minute before, and ~1 minute after, each predicted time. Because
# the duration of each pulse is less than the duration of an observation, this strategy
# is sufficient; it'll never be the case were I miss part of a pulse, even though I'm
# only checking if the start and end times fall within an observation
print("Cross-matching with predicted ToAs...")
obs_duration = 296 * u.s
pulse_duration = 2 * u.min

start_times = predicted_ToAs - pulse_duration/2
end_times = predicted_ToAs + pulse_duration/2

start_caught = np.array([np.any(np.logical_and(obsid <= start_times, start_times <= obsid + obs_duration)) for obsid in obsids], dtype=bool)
end_caught = np.array([np.any(np.logical_and(obsid <= end_times, end_times <= obsid + obs_duration)) for obsid in obsids], dtype=bool)

start_or_end_caught = np.logical_or(start_caught, end_caught)

#print(obsids.mjd)
#print(obsids[start_or_end_caught].gps.astype(int))

outfile = 'gpm2024_J1755-2527_observations.txt'
print(f"Saving results to {outfile}")
header = f'Created with:\n  {' '.join(sys.argv)}\n\nObsIDs:'
np.savetxt(outfile, obsids[start_or_end_caught].gps, fmt='%d', header=header)
