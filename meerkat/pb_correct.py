#/usr/bin/env python

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys

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

c = 3.e8*u.m/u.s
D = 13.5*u.m

coords_1 = SkyCoord(sys.argv[1], unit=(u.hour, u.deg), frame='fk5') # hmsdms
coords_2 = SkyCoord(sys.argv[2], unit=(u.hour, u.deg), frame='fk5') # hmsdms
nu = sys.argv[3]*u.MHz #MHz

# Calibrated to the plots of Villiers 2023
fwhm = 1.13*np.degrees(((c/nu) / D)*u.rad)
print(f"Beam FWHM = {fwhm}")
sigma = fwhm/2.355
print(f"Beam sigma = {sigma}")

sep = coords_1.separation(coords_2)
print(f"Distance from pointing centre = {sep.deg} deg")

print(gaussian_1d(sep.deg, 0, sigma.value))

from matplotlib import pyplot as plt

x = np.linspace(0, 1.5*fwhm.value)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(x, gaussian_1d(x, 0, sigma.value), label='primary beam', color='black')
ax.scatter(sep.deg, gaussian_1d(sep.deg, 0, sigma.value), label='source position', color='red')
ax.set_xlabel("Distance from boresight / deg")
ax.set_ylabel("Primary beam value (normalised to 1.0)")
ax.set_title(f"Primary beam at {nu}")
ax.legend()
fig.savefig("Primary beam correction.png", bbox_inches="tight")
