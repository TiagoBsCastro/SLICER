import numpy as np
import params
from print import print
from cosmology import cosmo, rho_m, getZ

############### EXTRA INPUT ###################
_NBINS = 20  # Number of radial Bins
_NPMIN = 100 # Minimum number of particles
_Npx   = 1   # Multiplicative factor for the minimum number of mass (see SAP)
###############################################

print('Reading PLC: {}'.format(params.name.format('plc')))
plc = np.load(params.name.format('plc'))

print('Computing the HMF in {} equal-volume radial bins'.format(_NBINS))
rmin = 0.0
rmax = cosmo.comoving_distance(params.zmax).value
print('Rmin.: {0:.3f}'.format(rmin))
print('Rmax.: {0:.3f}'.format(rmax))
if rmin == 0:
    rbins = np.linspace(0.0, rmax**3, _NBINS)**(1/3)
else:
    rbins = np.linspace(rmin**3, rmax**3, _NBINS)**(1/3)
vbins = cosmo.comoving_volume( getZ(rbins) ).value * params.skyfrac
vbins = vbins[1:] - vbins[:-1]

print("Reading the features of the individual halos")
mvir  = plc['f38'][plc['f5'] == -1]
r     = plc['f89'][plc['f5'] == -1]
z     = plc['f82'][plc['f5'] == -1] 
mpart = rho_m * params.boxsize**3/params.npart * 1e10
npart = np.round(mvir/mpart).astype(np.int32)

print("Binning in mass as in SAP")
print(" (see: https://github.com/TiagoBsCastro/SAP/blob/ML/scripts/HMF.py)")
NBINS   = 10*int(np.ceil(np.log10( (1e16//(mpart*_Npx))/_NPMIN ) )) + 1
mbins   = np.geomspace(_NPMIN, (1e16//(mpart*_Npx)  + 1), NBINS).astype(int)
mbins   = np.unique(mbins)
mbins  *= _Npx
dlogm   = np.log(mbins[1:]/mbins[:-1])
print('Mmin.: {0:.3f}'.format(mbins.min()))
print('Mmax.: {0:.3f}'.format(mbins.max()))

print('Looping over the shells')
for v, r1, r2 in zip(vbins, rbins[:-1], rbins[1:]):

    print("# Filtering halos {0:.3f}<=r<{1:.3f}".format(r1, r2))
    cut = (r >= r1) & (r < r2)
    print("# Building the histogram")
    counts, bins    = np.histogram(npart[cut], bins=mbins)
    mass, bins      = np.histogram(npart[cut], bins=mbins, weights=mvir[cut])
    mass[counts>0]  = (mass[counts>0]/counts[counts>0])
    mass[counts==0] = (mbins[:-1] * mbins[1:] / (mbins[:-1] + mbins[1:]) * 2 * mpart)[counts==0]

    redshift = z[cut].mean()

    print("# Saving the HMF at mean z={0:.3f}".format(redshift))
    np.save(params.name.format('hmf-a={0:.3f}'.format(1/(redshift + 1))), np.transpose([mass, counts/v/dlogm, mbins[:-1]*mpart, mbins[1:]*mpart]))
