import numpy as np
import pymaster as nmt
from astropy.io import fits
import matplotlib.pyplot as plt
import smr

# This script describes the functionality of the flat-sky version of pymaster

# Dimensions:
# First, a flat-sky field is defined by four quantities:
#  - Lx and Ly: the size of the patch in the x and y dimensions (in radians)
Lx = 10.0 * np.pi/180
Ly = 10.0 * np.pi/180


mpt  = fits.getdata("kappa_fft_rec.fits")
mpt2 = fits.getdata("gamma_fft.fits")
#  - Nx and Ny: the number of pixels in the x and y dimensions
Nx = mpt.shape[0]
Ny = mpt.shape[1]

print("Is close: {}".format(np.allclose(mpt, mpt2)))

mask = np.ones((Nx, Ny))

# Bins:
# For flat-sky fields, bandpowers are simply defined as intervals in ell, and
# pymaster doesn't currently support any weighting scheme within each interval.
l0_bins = np.arange(Nx/4) * 4 * np.pi/Lx
lf_bins = (np.arange(Nx/4)+1) * 4 * np.pi/Lx
b = nmt.NmtBinFlat(l0_bins, lf_bins)
# The effective sampling rate for these bandpowers can be obtained calling:
ells_uncoupled = b.get_effective_ells()

# Fields:
# Once you have maps it's time to create pymaster fields.
# Note that, as in the full-sky case, you can also pass
# contaminant templates and flags for E and B purification
# (see the documentation for more details)
res = []
for f0 in [nmt.field.NmtFieldFlat(Lx, Ly, mask, [mpt]), nmt.field.NmtFieldFlat(Lx, Ly, mask, [mpt2])]:

    # Workspaces:
    # As in the full-sky case, the computation of the coupling matrix and of
    # the pseudo-CL estimator is mediated by a WorkspaceFlat case, initialized
    # by calling its compute_coupling_matrix method:
    w00 = nmt.NmtWorkspaceFlat()
    w00.compute_coupling_matrix(f0, f0, b)

    # Computing power spectra:
    # As in the full-sky case, you compute the pseudo-CL estimator by
    # computing the coupled power spectra and then decoupling them by
    # inverting the mode-coupling matrix. This is done in two steps below,
    # but pymaster provides convenience routines to do this
    # through a single function call
    cl00_coupled = nmt.compute_coupled_cell_flat(f0, f0, b)

    # Let's look at the results!
    plt.plot(ells_uncoupled, cl00_coupled[0], label='Coupled')
    res += [cl00_coupled[0].tolist()]

plt.loglog()
plt.legend()
plt.show()

res = np.array(res)

plt.plot(ells_uncoupled, res[1]/res[0]-1)
plt.xscale('log')
plt.show()
