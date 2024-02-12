import numpy as np
import MAS_library as MASL
import g3read
from os.path import exists
from glob import glob
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

snapshots = glob('lensing/w-1p1/snapdir_*/*.0') #snapshot name
snapshots += glob('lensing/w-09/snapdir_*/*.0')
grid      = 16384                               #grid size
BoxSize   = 250.0                               #Mpc/h
ptypes    = [1]                                 #CDM + neutrinos
MAS       = 'TSC'                               #3th order interpolation jernel
verbose   = False                               #whether print information on the progress

for snapshot in snapshots:

    if not rank:

        print("Working on {}".format(snapshot[:-2]))
        if g3read.GadgetFile(snapshot).header.num_files != size:

            print("Number of files different than the number of ranks!")
            comm.Abort()

    comm.Barrier()
    # reading the particle positions
    pos = g3read.read_new(snapshot[:-1]+str(rank), "POS ", 1)

    for axis, name in zip([[0,1], [0,2], [1,2]], ['xy', 'xz', 'yz']):

        pos2d = pos[:, axis].astype(np.float32)

        # create the delta array
        delta = np.zeros((grid, grid), dtype=np.float32)
        
        # construct 2D density field
        MASL.MA(pos2d, delta, BoxSize, MAS, verbose=verbose)

        # the 'totals' array will hold the sum of each 'data' array
        if not rank:
            # only processor 0 will actually get the data
            totals = np.zeros_like(delta)
        else:
            totals = None

        # use MPI to get the totals 
        comm.Reduce(
            [delta, MPI.FLOAT],
            [totals, MPI.FLOAT],
            op = MPI.SUM,
            root = 0
        )

        if not rank:
            delta = totals
            # compute density contrast: delta = rho/<rho> - 1
            delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0
            np.save(snapshot[:-2]+"_"+name+".npy", delta)

        comm.Barrier()
