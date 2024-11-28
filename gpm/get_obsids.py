import numpy as np
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

# Define source location
src_coord = SkyCoord("17h55m34.9s -25d27m49.1s", frame='icrs')
print(f"Source coordinates: {src_coord}")

# Get obsids in GPM within 15 degrees of source (very coarse, will be refined later)
obs = []
for page in [3]:
    print(f"Retrieving page {page} from https://ws.mwatelescope.org/")
    response = requests.get(f"https://ws.mwatelescope.org/metadata/find/?pagesize=500&notdeleted=on&projectid=G0080&minra=234&maxra=304&mindec=-40&maxdec=-10&page={page}")
    obs += response.json()

obsids, _, _, _, RAs, Decs = list(map(list, zip(*obs)))

pointings = SkyCoord(RAs, Decs, unit=(u.deg, u.deg), frame='icrs')
separations = src_coord.separation(pointings)
print(separations)
