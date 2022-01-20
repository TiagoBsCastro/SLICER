import numpy as np
import params
from print import print
from cosmology import cosmo, rho_m, getZ
import healpy as hp
import matplotlib.pyplot as plt

############### EXTRA INPUT ###################
_NSIDE  = 512  # Resolution of the healpix map
_NBINS  = 20   # Number of radial Bins
###############################################
NPIX = hp.nside2npix(_NSIDE)

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

print("Creating mask")
tm, pm = hp.pixelfunc.pix2ang(_NSIDE, np.arange(NPIX))
mask = (tm<=params.fovinradians)

print("Starting loop over the layers")
for d1, d2 in zip(rbins[:-1], rbins[1:]):

    print("## Filtering halos in layer {0:.3f}<=r<{1:.3f}".format(d1, d2))
    cut = (plc['f89'] >= d1) & (plc['f89'] < d2)
    z = plc['f82'][cut].mean()

    print("## Pixelizing the halos")
    indices = hp.ang2pix(_NSIDE, plc['f90'][cut], plc['f91'][cut])
    idx, counts_hp = np.unique(indices, return_counts=True)

    print("## Creating Healpix map")
    hpx_map = np.zeros(NPIX,dtype=int)
    hpx_map[idx] = counts_hp
    print("## Pixel occupation, [min, max, mean, total]: [",round(hpx_map.min(),2), round(hpx_map.max(),2), round(hpx_map.mean(),2), hpx_map.sum(), "]")
    mean_pix=np.mean(hpx_map[mask])
    delta=hpx_map/mean_pix-1
    delta[~mask] = hp.UNSEEN

    print("## Saving map at z={:.3f}".format(z))
    hp.mollview(delta, unit=r"$\delta$", min=np.quantile(delta[mask], 0.1), max=np.quantile(delta[mask], 0.9), title="z={:.3f}".format(z))
    plt.savefig("mollweide-z={:.3f}.pdf".format(z), bbox_inches='tight')
    plt.close()
