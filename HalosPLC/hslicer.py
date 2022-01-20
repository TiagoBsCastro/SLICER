import numpy as np
import numexpr as ne
import params
from print import print
######################## Cosmology ########################
print("Setting cosmology")
from cosmology import cosmo, getZ
######################## Geometry ########################
print("Setting replication geometry")
from geometry import geometry
print("Saving Geometry: {}".format(params.name.format('geometry')))
np.save(params.name.format('geometry'), geometry)
##########################################################
####                         Main                     ####
##########################################################

print("Building the layers")
layerdistances = cosmo.comoving_distance(params.zsnapshots[1:]/2 + params.zsnapshots[:-1]/2)
dmax           = cosmo.comoving_distance(params.zmax)
layerdistances = layerdistances[layerdistances<dmax]
layerdistances = np.append(layerdistances, dmax)
layerdistances = np.append(0*dmax, layerdistances)

print("Looping over the layers")
for i, (d1, d2, z) in enumerate(zip(layerdistances[:-1], layerdistances[1:], params.zsnapshots)):

    print("# Reading Snapshot {0:03d}/{1:03d} at z={2:4.3f}".format(i, layerdistances.size-2, z))
    # Reading Snapshot
    snapi = np.genfromtxt(params.snapshots[i], dtype=params.dtype)

    if i == 0:

	# Copying the input dtype and adding some extra features
        outtype  = [('f{}'.format(x), snapi.dtype[x]) for x in range(len(snapi.dtype))] 
        outtype += [('f82', float)]   # Redshifts
        outtype += [('f83', int)]     # Repetition on x
        outtype += [('f84', int)]     # Repetition on y
        outtype += [('f85', int)]     # Repetition on z
        outtype += [('f86', float)]   # x
        outtype += [('f87', float)]   # y
        outtype += [('f88', float)]   # z
        outtype += [('f89', float)]   # r
        outtype += [('f90', float)]   # theta
        outtype += [('f91', float)]   # phi

        plc = np.empty(0, outtype)

    repinside = geometry[ ((d1.value>=geometry['nearestpoint']) & (d1.value<geometry['farthestpoint'])) | ((d2.value>=geometry['nearestpoint']) & (d2.value<geometry['farthestpoint'])) ]
    print("# Total number of replications inside the layer {}".format(repinside.shape[0]))
    print("# Looping over the Replications")
    print("##                          x   y   z      nearest        farthest")
    for repi in repinside:

        print("## Working on replication {}".format(repi))

        x = snapi['f17'] + repi['x'] * params.boxsize - params.boxsize/2
        y = snapi['f18'] + repi['y'] * params.boxsize - params.boxsize/2
        z = snapi['f19'] + repi['z'] * params.boxsize - params.boxsize/2

        print("## Computing angular coordinates")
        xy2 = ne.evaluate('x**2 + y**2')
        r = ne.evaluate('sqrt(xy2 + z**2)') # radius
        t = ne.evaluate('arccos(z/r)') # theta [0, pi]
        p = ne.evaluate('arctan2(y, x)+pi', local_dict={'pi':np.pi}) # phi [0, 2pi]

        print("## Creating angular and layer filter")
        cut = (r>=d1.value) & (r<d2.value) & (t<=params.fovinradians)
        nhalos = cut.sum()
        print("## Total number of halos inside the PLC layer: {}".format(nhalos))

        # Inflating the plc array
        plci = np.empty(nhalos, dtype=outtype)
        for i in range(82):

            plci['f{}'.format(i)] = snapi[cut]['f{}'.format(i)]

        plci['f82'] = getZ(r[cut])
        plci['f83'] = repi['x']
        plci['f84'] = repi['y']
        plci['f85'] = repi['z']
        plci['f86'] = x[cut]
        plci['f87'] = y[cut]
        plci['f88'] = z[cut]
        plci['f89'] = r[cut]
        plci['f90'] = t[cut]
        plci['f91'] = p[cut]

        print("## Appending to the PLC")
        plc = np.append(plc, plci)

    print("## All done for this snap\n")

print("Saving PLC: {}".format(params.name.format('plc')))
np.save(params.name.format('plc'), plc)
print("All done!")
