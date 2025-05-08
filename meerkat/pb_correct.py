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
  return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma)**2)

c = 3.e8*u.m/u.s
D = 13.5*u.m

coords_1 = SkyCoord(sys.argv[1], unit=(u.hour, u.deg), frame='fk5') # hmsdms
coords_2 = SkyCoord(sys.argv[2], unit=(u.hour, u.deg), frame='fk5') # hmsdms
nu = sys.argv[3]*u.MHz #MHz

# Calibrated to the plots of Villiers 2023
fwhm = 1.13*np.degrees(((c/nu) / D)*u.rad)
theta = fwhm/2.355

sep = coords_1.separation(coords_2)

print(gaussian_1d(sep.deg, 0, theta.value))
