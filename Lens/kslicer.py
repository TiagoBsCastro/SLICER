from astropy.io import fits
import numpy as np
from scipy.integrate import quad
from astropy.cosmology import Flatw0waCDM
from astropy.cosmology import z_at_value
from astropy.constants import G as GNewton
from astropy.constants import c as speedOfLight
from colossus.cosmology import cosmology
import astropy.units as u
from glob import glob
import sys
import input
import utils
import smr

############################ Cosmology ###############################

# Background quantities are computed with Astropy
cosmo   = Flatw0waCDM(Om0=input.Om0, H0 = 100.0, w0=input.w0, wa=input.wa) # Mpc/h units
c2OverG = (speedOfLight**2/GNewton).to(u.Msun/u.Mpc).value/1e10            # c^2/G in 10^10Msun/h/(Mpc/h)
rho_m   = cosmo.Om0 * cosmo.critical_density0.to('Msun/Mpc^3').value/1e10  # <rho_m> in 10^10 Msun/h/(Mpc/h)^3
# Perturbation quantities are computed with Colossus
pert = cosmology.setCosmology('MyCosmo', input.pert_params)

########################### Mass Maps ################################

mapsNames = sorted(glob( input.massMapsNames.substitute(mapIndex=input.imap) ))
mapsNames = np.array(mapsNames)
zsnap     = np.loadtxt( input.planesFileNames.substitute(mapIndex=input.imap), usecols=[6] )

################ Reconstruction of the Lens Planes ###################

dllow = np.empty((len(mapsNames), 2), dtype=np.float)
dlup  = np.empty_like(dllow)
dllen = np.empty_like(dllow)

for i, fname in enumerate(mapsNames):

    dllow[i,1] = fits.getheader(fname)['DLLOW'] * input.H0/100
    dlup[i,1]  = fits.getheader(fname)['DLUP'] * input.H0/100

    dllow[i,0] = z_at_value(cosmo.comoving_distance, dllow[i,1]*u.Mpc) if dllow[i,1] > 0.0 else 0.0
    dlup[i,0]  = z_at_value(cosmo.comoving_distance,  dlup[i,1]*u.Mpc)

print(dlup)
# Defining the source redshifts for integration
if input.source_redshifts == 'All':

    ztab = dlup[:,0]

else:

    ztab = dlup[input.source_redshifts,0]

# Testing the Lens Planes
utils.planesSanityCheck(dllow, dlup)

# Computing the Lens Planes (Distance Average)
for i, (d1, d2) in enumerate(zip(dllow, dlup)):

    dllen[i, 0]  = quad( lambda z: z * cosmo.comoving_distance(z).value, dllow[i, 0], dlup[i, 0] )[0]
    dllen[i, 0] /= quad( lambda z:     cosmo.comoving_distance(z).value, dllow[i, 0], dlup[i, 0] )[0]
    dllen[i, 1]  = cosmo.comoving_distance(dllen[i, 0]).value

print( "Plane list:     Dl_           Lens             Dl^"  )
for d1, dl, d2 in zip(dllow[:,1], dllen[:,1], dlup[:,1]):

    print("           {0: >10.4f}     {1: >10.4f}      {2: >10.4f} ".format(d1, dl, d2 ))

#################### Starting to Sum the Maps #######################

for zs in ztab:

    print( "\nIntegrating the convergence map: zs={:5.4f}\n".format(zs) )
    planeIndexes = dlup[:,0]<= zs + 1e-4
    maps = mapsNames[planeIndexes]
    kappa = np.zeros_like( fits.getdata(maps[0]) )
    header = fits.getheader(maps[0])

    for i, fname in enumerate(maps):

        print( "\t{}".format(fname) )
        lensKernel = utils.lensKernel( dllen[i,0], zs, cosmo )
        A          = utils.A( dllen[i,1], header['PHYSICALSIZE'] )
        A         /= header['NAXIS1'] * header['NAXIS2']
        mass       = fits.getdata(fname)/A
        mean       = rho_m * (dlup[i,1] - dllow[i,1])

        growthCorr = pert.growthFactor(dllen[i,0])/pert.growthFactor(zsnap[i])
        print( "\t<rho>_theory/<rho>_map={}".format(mean/mass.mean()) )
        print( "\tGroth Factor Correction={}\n".format(growthCorr) )
        kappa += 4.0 * np.pi/c2OverG * lensKernel * ( mass - mass.mean() ) * growthCorr * (1.0+dllen[i,0])**2

    hdu = fits.PrimaryHDU()
    hdu.data = kappa
    utils.genericHeader(hdu.header, header['NAXIS1'] , zs, header['PHYSICALSIZE'])
    print("\nSaving the Map")
    hdu.writeto( input.kappaMapsNames.substitute(mapIndex=input.imap,zs=np.round(zs,4)) , overwrite=True)
    print("Computing shear maps")
    smr.smr( input.kappaMapsNames.substitute(mapIndex=input.imap,zs=np.round(zs,4)) )
    print("Done\n")
