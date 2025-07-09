#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "05/04/2025"

import os
import sys
import shutil

import matplotlib
matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['serif']})

from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.nddata import Cutout2D

import numpy as np

# Didn't work
#from astropy.visualization import make_lupton_rgb
#rgb_default = make_lupton_rgb(red, green, blue, filename="temp.jpeg")

from astropy.visualization import PercentileInterval
from astropy.visualization import AsinhStretch

def normalize(arr, vmin, vmax):
    nor = (arr - vmin) / (vmax - vmin)
    nor[np.where(nor<0.0)] = 0.0
    nor[np.where(nor>1.0)] = 1.0
    return nor

#coords_old = SkyCoord("17:55:34.87 -25:27:49.1", frame="fk5", unit=(u.hourangle, u.deg))
# Offset is about 3.5" North
# Maybe 3" East
#coords = coords_old.spherical_offsets_by(3*u.arcsec, 3.5*u.arcsec)
# Highest S/N measurement -- but worst PSF
#coords_m1= SkyCoord("07:15:07.37 -11:41:53.72", frame="fk5", unit=(u.hourangle, u.deg))
# Lower S/N measurements -- with better PSF
#coords_m2= SkyCoord("07:15:07.38 -11:41:53.59", frame="fk5", unit=(u.hourangle, u.deg))
#coords_m3= SkyCoord("07:15:07.35 -11:41:53.81", frame="fk5", unit=(u.hourangle, u.deg))
# Average of those two
#coord_askap_orig= SkyCoord("07:15:07.365 -11:41:53.745", frame="fk5", unit=(u.hourangle, u.deg))
# Unweighted average from three measurements
#coord_askap_orig= SkyCoord("07:15:07.3704 -11:41:53.7072", frame="fk5", unit=(u.hourangle, u.deg))
# Weighted average of 16 measurements at UHF

coord_askap_orig= SkyCoord("17:55:34.87 -25:27:49.1", frame="fk5", unit=(u.hourangle, u.deg))
err_askap_a = Angle('00h00m00.1s').to('deg')
err_askap_b = Angle('00d00m00.1s').to('deg')

coord_mkt= SkyCoord("17:55:34.87 -25:27:49.94", frame="fk5", unit=(u.hourangle, u.deg))
err_mkt_a = Angle('00h00m00.02s').to('deg')
err_mkt_b = Angle('00d00m00.4s').to('deg')

# Weighted average of 3 measurements at L-band
#coords_mktl= SkyCoord("07h15m07.37365723s -11d41m53.28048706s", frame="fk5", unit=(u.hourangle, u.deg))
#err_mktl = 0.33*u.arcsec
# Single highest-S/N measurement
coords_mktl1= SkyCoord("17:55:34.87 -25:27:49.1", frame="fk5", unit=(u.hourangle, u.deg))
#print(coords.to_string('hmsdms', sep=":"))
#print(coords.to_string('hmsdms'))
#err_ra_mktl1 = 3.7291164e-5*u.deg
#err_dec_mktl1 = 3.3163615E-5*u.deg


fits_red = "k.fits"
fits_green = "h.fits"
fits_blue = "j.fits"

red = fits.open(fits_red)
green = fits.open(fits_green)
blue = fits.open(fits_blue)

# Percentile interval for image normalisation
pct = 99.999
interval = PercentileInterval(pct)

stretch = AsinhStretch(a=0.1)

framesize = 0.3*u.arcmin

w_red = wcs.WCS(red[1].header)
w_green = wcs.WCS(green[1].header)
w_blue = wcs.WCS(blue[1].header)

red_data = Cutout2D(red[1].data, coord_mkt, framesize, wcs = w_red)
green_data = Cutout2D(green[1].data, coord_mkt, framesize, wcs = w_green)
blue_data = Cutout2D(blue[1].data, coord_mkt, framesize, wcs = w_blue)

i = interval.get_limits(red_data.data)
r = stretch(normalize(red_data.data, *i))
i = interval.get_limits(green_data.data)
g = stretch(normalize(green_data.data, *i))
i = interval.get_limits(blue_data.data)
b = stretch(normalize(blue_data.data, *i))

rgb = np.dstack([r,g,b])

# Should be the same for all images thanks to Cutout2D

fig = plt.figure(figsize=(5,5))

xstart = 0.2 # fraction of frame
ystart = 0.15 # fraction of frame
size = 0.7 # fraction of frame
size = 0.7 # fraction of frame
ax = fig.add_axes([xstart,ystart,size,size], projection = red_data.wcs)
img = ax.imshow(rgb,origin="lower")


# Grid overlay
#overlay = ax.get_coords_overlay("fk5")
#overlay.grid(axes = ax, color='white', ls = 'dotted')
#overlay["ra"].set_ticklabel_visible(False)
#overlay["dec"].set_ticklabel_visible(False)

# Axis labels
lon = ax.coords['ra']
lon.set_axislabel("Right Ascension (J2000)")
#lon.set_major_formatter('hh:mm')
lat = ax.coords['dec']
lat.set_axislabel("Declination (J2000)")
#lat.set_major_formatter('dd:mm')

# Transient source
# MeerKAT
width = size * err_mkt_a / framesize # in "fraction of the frame"
height = size * err_mkt_b / framesize # in "fraction of the frame"

ell = Ellipse((coord_mkt.ra.value, coord_mkt.dec.value), 
              width=2*err_mkt_a.value, height=2*err_mkt_b.value, angle=180-155.03435,
              edgecolor='#0f0',
              facecolor='none',
              transform=ax.get_transform("world")
              )
ax.add_patch(ell)
ax.invert_yaxis()

#rm = SphericalCircle((coord_askap_orig.ra, coord_askap_orig.dec), err_mkt, edgecolor='magenta', facecolor='none', transform=ax.get_transform("world"), alpha=0.5)
#ax.add_patch(rm)
#rm = SphericalCircle((coords_mktl.ra, coords_mktl.dec), err_mktl, edgecolor='black', facecolor='none', transform=ax.get_transform("world"), alpha=0.8)
#ax.add_patch(rm)

#ax.errorbar(coords_mktl1.ra, coords_mktl1.dec, xerr=err_ra_mktl1, yerr=err_dec_mktl1, color='red', transform=ax.get_transform("world"))
#print(err_ra_mktl1.to(u.arcsec))
#print(err_dec_mktl1.to(u.arcsec))

#r = SphericalCircle((coords_old.ra, coords_old.dec), 4*u.arcsec, edgecolor='red', facecolor='none', transform=ax.get_transform("world"))
#ax.add_patch(r)
# M81: D25
#m81_coords = SkyCoord("09:55:33.17306 +69:03:55.0610", frame="fk5", unit=(u.hourangle, u.deg))
# From Gemma:
# Log D25 major axis = 2.43 (in units of 0.1 arcmins)
# Log R25 (ratio between major and minor D25 axis) = 0.28
# The PA on this ellipse (I think!) = 157 deg (north through east)
#m81_d25_a = (10**2.43)/20. # arcmin
#m81_d25_b = m81_d25_a/(10**0.28) # arcmin
#m81_d25_b = 0.28*m81_d25_a # arcmin
#m81_d25_pa = 157. # deg
#framesize = 20. # arcmin

# But plotting an ellipse in matplotlib is actually the one thing that doesn't transform well -- you get infinite error messages if the a and b are in degrees. They have to be in fraction of the frame.

#width = size * m81_d25_a / framesize # in "fraction of the frame"
#height = size * m81_d25_b / framesize # in "fraction of the frame"

#ax.add_patch(Ellipse((m81_coords.ra.value, m81_coords.dec.value), 
#                     width=width, height=height, angle=m81_d25_pa,
#                     width=m81_d25_a, height=m81_d25_b, angle=m81_d25_pa, 
#                     edgecolor='green',
#                     facecolor='none',
#                     transform=ax.get_transform("world")
#                    ))

# Scatter makes the plot zoom out for some reason
#ax.set_xlim(0,xmax)
#ax.set_ylim(0,ymax)

fig.savefig("J1755_nir.pdf", dpi=1000, bbox_inches="tight")
#fig.savefig("M81_SDSS.eps", dpi=1000)
#fig.savefig("M81_SDSS.pdf")
