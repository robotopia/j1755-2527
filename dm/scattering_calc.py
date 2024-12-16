import numpy as np
import pygedm
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord

src = SkyCoord("17:55:34.9", "-25:27:49.1", unit=(u.hourangle, u.deg), frame='icrs')
gl = src.galactic.l
gb = src.galactic.b

DMunit = u.pc / u.cm**3

D = 4148.8 * u.s
nu1 = 185 * u.MHz
DM_at_887MHz = 600 * DMunit

apparent_DM_at_185MHz = 1221 * DMunit
apparent_slope_at_185MHz = -2*D*apparent_DM_at_185MHz/DMunit *(nu1/u.MHz)**(-3) / u.MHz

DM_slope_at_887MHz = -2*D*DM_at_887MHz/DMunit * (887*u.MHz/u.MHz)**(-3) / u.MHz
extrapolated_DM_slope = -2*D*DM_at_887MHz/DMunit * (nu1/u.MHz)**(-3) / u.MHz
inferred_scattering_slope = apparent_slope_at_185MHz - extrapolated_DM_slope

_, tau_sc_1GHz_ne2001 = pygedm.dm_to_dist(gl.deg, gb.deg, DM_at_887MHz, method='ne2001')
_, tau_sc_1GHz_ymw16 = pygedm.dm_to_dist(gl.deg, gb.deg, DM_at_887MHz, method='ymw16')

beta = -4.0
scattering_slope_ne2001_at_887MHz = beta * tau_sc_1GHz_ne2001 * (887*u.MHz/u.GHz)**(beta - 1) / u.GHz
scattering_slope_ne2001 = beta * tau_sc_1GHz_ne2001 * (nu1/u.GHz)**(beta - 1) / u.GHz
scattering_slope_ymw16 = beta * tau_sc_1GHz_ymw16 * (nu1/u.GHz)**(beta - 1) / u.GHz

print(f"{DM_slope_at_887MHz.to('s MHz-1') = }")
print(f"{apparent_slope_at_185MHz.to('s MHz-1') = }")
print(f"{extrapolated_DM_slope.to('s MHz-1') = }")
print(f"{inferred_scattering_slope.to('s MHz-1') = }")
print(f"{scattering_slope_ne2001_at_887MHz.to('s MHz-1') = }")
print(f"{scattering_slope_ne2001.to('s MHz-1') = }")
print(f"{scattering_slope_ymw16.to('s MHz-1') = }")
