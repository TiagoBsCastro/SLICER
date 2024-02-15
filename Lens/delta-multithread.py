import numpy as np
import MAS_library as MASL
import readsnap as rs
from os.path import exists
from glob import glob
from multiprocessing import Pool

######################################### INPUT ##################################################
snapshots = [f'/dataA02/asaro/Clusters1x/D1/snap_{i:03d}' for i in range(92) ] #snapshot name
grid      = 4096                                                               #grid size
ptypes    = [0, 1, 2, 3, 4, 5]                                                 #CDM + Gas + Low-Res + Stars + BH
MAS       = 'CIC'                                                              #3th order interpolation kernel
verbose   = False                                                              #whether print information on the progress
size      = 1_000_000                                                          #plane size
thickness = 100_000                                                            #lens thickness
num_procs = 48
##################################################################################################

def density_map (snapshot):

    print("Working on {}".format(snapshot))
    # reading the particle positions
    npart = rs.snapshot_header(snapshot).npart
    BoxSize = rs.snapshot_header(snapshot).boxsize
    pos  = np.empty((npart.sum(), 3), dtype=np.float32)
    mass = np.empty(npart.sum(), dtype=np.float32)

    # Getting the main group position
    dir = '/'.join(snapshot.split("/")[:-1])
    snap = snapshot.split("_")[-1]
    GPOS = rs.read_block(snapshot.split("snap")[0]+"/groups_090/sub_090", "GPOS")[0]

    nparti = 0  
    for p in ptypes:

        pos[nparti:nparti+npart[p]]  = rs.read_block(snapshot, "POS ", p)
        mass[nparti:nparti+npart[p]] = rs.read_block(snapshot, "MASS", p)
        nparti += npart[p]

    for axis, name in zip([[0,1], [0,2], [1,2]], ['xy', 'xz', 'yz']):

        # Filtering around the center
        perp = np.setxor1d([0, 1, 2], axis)[0]
        sel  = (pos[:, axis[0]] - GPOS[axis[0]] <= size/2) & (pos[:, axis[0]] - GPOS[axis[0]]>= -size/2) &\
		   (pos[:, axis[1]] - GPOS[axis[1]] <= size/2) & (pos[:, axis[1]] - GPOS[axis[1]]>= -size/2) &\
		   (pos[:, perp] - GPOS[perp] <= thickness/2) & (pos[:, perp] - GPOS[perp]>= -thickness/2)

        pos2d = pos[sel][:, axis] - (GPOS[axis].astype(np.float32) - size/2)
        W = mass[sel]

        # create the delta array
        delta = np.zeros((grid, grid), dtype=np.float32)

        # construct 2D density field
        MASL.MA(pos2d, delta, size, MAS, verbose=verbose, W=W)

        delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0
        fout = "DensityMaps/"+snapshot.split("/")[-1]+"_"+name+".npy"
        np.save(fout, delta)

with Pool(processes=num_procs) as pool:
        # Map the list of snapshots to the processing function
        pool.map(density_map, snapshots)
