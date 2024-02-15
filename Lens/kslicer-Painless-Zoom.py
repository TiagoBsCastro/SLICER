import astropy.units as u
import numpy as np
import sys
import utils
import smr
from time import time
from scipy.integrate import quad
from astropy.io import fits
from astropy.cosmology import Flatw0waCDM, z_at_value 
from astropy.constants import G as GNewton
from astropy.constants import c as speedOfLight
from colossus.cosmology import cosmology
from glob import glob

############################## Input #################################

# Cosmology 
Om0 = 0.2400
Ob0 = 0.0370
H0  = 100.00
w0  = -1.000
wa  = 0.0000
s8  = 0.8000
ns  = 0.9600

# Output file
output_file = 'outputs.txt'

# Sim directory
base_dir = './'

# Sim details
lens_thickness = 100.0   # Mpc/h
lens_size      = 1_000.0 # Mpc/h

# Redshift source
zs = 1100.0

# pyFFTW threads
threads = 4

############################ Cosmology ###############################

# Background quantities are computed with Astropy
cosmo   = Flatw0waCDM(Om0=Om0, H0 = 100.0, w0=w0, wa=wa) # Mpc/h units
c2OverG = (speedOfLight**2/GNewton).to(u.Msun/u.Mpc).value/1e10            # c^2/G in 10^10Msun/h/(Mpc/h)
rho_m   = cosmo.Om0 * cosmo.critical_density0.to('Msun/Mpc^3').value/1e10  # <rho_m> in 10^10 Msun/h/(Mpc/h)^3
# Perturbation quantities are computed with Colossus
params = {'flat': True, 'H0': H0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': s8, 'ns': ns, 'relspecies':False, 'de_model':'w0wa', 'w0':w0, 'wa':wa}
pert = cosmology.setCosmology('MyCosmo', params)

########################### Mass Maps ################################

mapsNames = sorted(glob( base_dir + 'DensityMaps/snap*npy' ))
mapsNames = np.array(mapsNames)[::-1]
mapsNames = mapsNames.reshape(mapsNames.size//3, 3)
mapsNames = np.array(mapsNames)[1::]
zsnap     = 1.0/np.loadtxt( output_file, comments=',' )[:-1] - 1
zsnap     = zsnap[::-1]

################ Reconstruction of the Lens Planes ###################

dllow = np.empty((len(mapsNames), 2), dtype=np.float32)
dlup  = np.empty_like(dllow)
dllen = np.empty_like(dllow)

for i, fname in enumerate(mapsNames):

    print(i, fname)
    dllow[i,1] = cosmo.comoving_distance(zsnap[i]).to('Mpc').value
    dlup[i,1]  = dllow[i,1] + lens_thickness

    dllow[i,0] = z_at_value(cosmo.comoving_distance, dllow[i,1]*u.Mpc) if dllow[i,1] > 0.0 else 0.0
    dlup[i,0]  = z_at_value(cosmo.comoving_distance,  dlup[i,1]*u.Mpc)

# Defining the source redshifts for integration
ztab = dlup[:,0]

# Computing the Lens Planes (Distance Average)
for i, (d1, d2) in enumerate(zip(dllow, dlup)):

    dllen[i, 0]  = quad( lambda z: z * cosmo.comoving_distance(z).value, dllow[i, 0], dlup[i, 0] )[0]
    dllen[i, 0] /= quad( lambda z:     cosmo.comoving_distance(z).value, dllow[i, 0], dlup[i, 0] )[0]
    dllen[i, 1]  = cosmo.comoving_distance(dllen[i, 0]).value

print( "Plane list:     Dl_           Lens             Dl^"  )
for d1, dl, d2 in zip(dllow[:,1], dllen[:,1], dlup[:,1]):

    print("           {0: >10.4f}     {1: >10.4f}      {2: >10.4f} ".format(d1, dl, d2 ))

#################### Starting to Sum the Maps #######################

print( "\nIntegrating the convergence map: zs={:5.4f}\n".format(zs) )
planeIndexes = dlup[:,0]<= zs + 1e-4
maps = mapsNames[planeIndexes]
kappa = np.zeros_like( np.load( maps[0, 0] ) )

for i, fname in enumerate(maps):


    for ii in range(3):

        print( "\t{}".format(fname[ii]) )
        lensKernel = utils.lensKernel( dllen[i,0], zs, cosmo )
        mass       = (np.load(fname[ii]) + 1) * rho_m * lens_size**2 * lens_thickness
        mass      /= mass.shape[0] * mass.shape[1]
        mass      /= (lens_size**2/mass.shape[0]/mass.shape[1])
        mean       = rho_m * (dlup[i,1] - dllow[i,1])

        growthCorr = pert.growthFactor(dllen[i,0])/pert.growthFactor(zsnap[i])
        print( "\t<rho>_theory/<rho>_map={}".format(mean/mass.mean()) )
        print( "\tGrowth Factor Correction={}\n".format(growthCorr) )
        kappa = 4.0 * np.pi/c2OverG * lensKernel * ( mass - mass.mean() ) * growthCorr * (1.0+dllen[i,0])**2

        hdu = fits.PrimaryHDU()
        hdu.data = kappa
        ang_size = -1 # -1 to flag that this is not produced by SLICER
        utils.genericHeader(hdu.header, kappa.shape[0] , zs, ang_size)
        print("\nSaving the Map")
        hdu.writeto( base_dir + "LensMaps/kappa_{0:5.4f}_{1:d}.fits".format(dllen[i, 0], ii) , overwrite=True)
        print("Computing shear maps")
        smr.smr( base_dir + "LensMaps/kappa_{0:5.4f}_{1:d}.fits".format(dllen[i, 0], ii), threads=threads)
        print("Done\n")
