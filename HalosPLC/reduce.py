import numpy as np

############### input ################

fname = "03/piccolo-C0-plc.npy"

######################################

# Blocks
# id, rvir, rs, mvir, m200b, m200c, m500c, m2500c, b2a, c2a, redshift, x, y, z, vx, vy, vz, r, theta, phi
infos  = [1, 11, 12, 38, 39, 40, 41, 42, 51, 52, 82, 86, 87, 88, 20, 21, 22, 89, 90, 91]
blocks = ['id', 'rvir', 'rs', 'mvir', 'm200b', 'm200c', 'm500c', 'm2500c', 'b2a', 'c2a', 'redshift', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'r', 'theta', 'phi']
# PID block to filter only structures
_pid   = 5

# Loading all catalog
cat = np.load(fname)
pid = cat[f'f{_pid}']
nhalos = (pid == -1).sum()

# Getting the dtype of the reduced array
dtype = np.dtype([(blk, cat.dtype[i]) for i, blk in zip(infos, blocks)])

# Getting only halos
cati = np.empty(nhalos, dtype=dtype)

for i, blk in zip(infos, blocks):

    cati[blk] = cat[f'f{i}'][pid==-1]

np.save(fname.replace('.npy', '-reduced.npy'), cati)
