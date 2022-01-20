import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import z_at_value
from astropy.constants import G as GNewton
from astropy.constants import c as speedOfLight
from params import Om0, zmax
from scipy.interpolate import interp1d

# Background quantities are computed with Astropy
cosmo   = FlatLambdaCDM(Om0=Om0, H0 = 100.0)                              # Mpc/h units
c2OverG = (speedOfLight**2/GNewton).to(u.Msun/u.Mpc).value/1e10           # c^2/G in 10^10Msun/h/(Mpc/h)
rho_m   = cosmo.Om0 * cosmo.critical_density0.to('Msun/Mpc^3').value/1e10 # <rho_m> in 10^10 Msun/h/(Mpc/h)^3

# Getting an interpolator for the distances and redshifts
_z = np.linspace(0, zmax+1, 200)
_d = cosmo.comoving_distance(_z).value

getZ = interp1d(_d, _z, fill_value="extrapolate")
