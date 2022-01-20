import numpy as np
import params
import cosmology

_N = 100

# List of Vertices of a unit cube centered at 0.0, 0.0, 0.0
_vertices  = np.array([[0.,0.,0.],[0.,0.,1.],
                       [0.,1.,0.],[0.,1.,1.],
                       [1.,0.,0.],[1.,0.,1.],
                       [1.,1.,0.],[1.,1.,1.]]) - 0.5

# List of the position of the faces centers of a unit cube centered at 0.0, 0.0, 0.0
_faces     =   np.array([[0.0,0.0,-0.5],[0.0,0.0,0.5],
                         [0.0,-0.5,0.0],[0.0,0.5,0.0],
                         [-0.5,0.0,0.0],[0.5,0.0,0.0]])

# List of the position of the "arestas" centers of a unit cube centered at 0.0, 0.0, 0.0
_arestas   =   np.array([[ 0.0,-0.5,-0.5],[ 0.0, 0.5,-0.5],[-0.5, 0.0,-0.5],[ 0.5, 0.0,-0.5],
                         [ 0.0,-0.5, 0.5],[ 0.0, 0.5, 0.5],[-0.5, 0.0, 0.5],[ 0.5, 0.0, 0.5],
                         [-0.5,-0.5, 0.0],[ 0.5,-0.5, 0.0],
                         [-0.5, 0.5, 0.0],[ 0.5, 0.5, 0.0]])

def getSurface():
    '''
    Returns N^2 points over the surface of the main repetition
    '''
    return np.concatenate((
    np.array( [[   x,   y,-0.5] for x in np.linspace(-0.5, 0.5,_N) for y in np.linspace(-0.5,0.5,_N) ]),
    np.array( [[   x,   y, 0.5] for x in np.linspace(-0.5, 0.5,_N) for y in np.linspace(-0.5,0.5,_N) ]),
    np.array( [[   x,-0.5,   y] for x in np.linspace(-0.5, 0.5,_N) for y in np.linspace(-0.5,0.5,_N) ]),
    np.array( [[   x, 0.5,   y] for x in np.linspace(-0.5, 0.5,_N) for y in np.linspace(-0.5,0.5,_N) ]),
    np.array( [[-0.5,   x,   y] for x in np.linspace(-0.5, 0.5,_N) for y in np.linspace(-0.5,0.5,_N) ]),
    np.array( [[ 0.5,   x,   y] for x in np.linspace(-0.5, 0.5,_N) for y in np.linspace(-0.5,0.5,_N) ])))

_surface = getSurface()

def getBorder():
    '''
    Returns N^2 points over the border of the main repetition
    '''

    return np.concatenate(
       (np.array( [[   x,-0.5,-0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[   x, 0.5,-0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[-0.5,   x,-0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[ 0.5,   x,-0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[   x,-0.5, 0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[   x, 0.5, 0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[-0.5,   x, 0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[ 0.5,   x, 0.5] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[-0.5,-0.5,   x] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[ 0.5,-0.5,   x] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[-0.5, 0.5,   x] for x in np.linspace(-0.5, 0.5,_N) ]),
        np.array( [[ 0.5, 0.5,   x] for x in np.linspace(-0.5, 0.5,_N) ]) ))

_border = getBorder()

def getVertices(repx,repy,repz):
    '''
    Returns the vertice's positions of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    xc = repx
    yc = repy
    zc = repz

    return np.array( [vertice + [xc, yc, zc] for vertice in _vertices])

def getArestas(repx, repy, repz):
    '''
    Returns the aresta's positions of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''

    xc = repx
    yc = repy
    zc = repz

    return np.array( [aresta + [xc, yc, zc] for aresta in _arestas])

def getFaces(repx, repy, repz):
    '''
    Returns the face's positions of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    xc = repx
    yc = repy
    zc = repz

    return np.array( [face + [xc, yc, zc] for face in _faces])

def getClosestAresta(repx, repy,repz):
    '''
    Returns the Closest aresta's position of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    arestas   = getArestas(repx, repy, repz)
    distances = np.sqrt( (arestas**2).sum(axis=1) )

    return arestas[distances.argmin()]


def getFarthestAresta(repx, repy,repz):
    '''
    Returns the Farthest aresta's position of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    arestas   = getArestas(repx, repy, repz)
    distances = np.sqrt( (arestas**2).sum(axis=1) )

    return arestas[distances.argmax()]


def getClosestFace(repx, repy,repz):
    '''
    Returns the Closest face's center position of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    faces     = getFaces(repx, repy, repz)
    distances = np.sqrt( (faces**2).sum(axis=1) )

    return faces[distances.argmin()]


def getFarthestFace(repx, repy,repz):
    '''
    Returns the Farthest Face's center position of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    faces     = getFaces(repx, repy, repz)
    distances = np.sqrt( (faces**2).sum(axis=1) )

    return faces[distances.argmax()]


def getClosestVertice(repx, repy,repz):
    '''
    Returns the Closest vertice's positions of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    vertices  = getVertices(repx, repy, repz)
    distances = np.sqrt( (vertices**2).sum(axis=1) )

    return vertices[distances.argmin()]


def getFarthestVertice(repx, repy,repz):
    '''
    Returns the Farthest vertice's positions of the repetition
    with cordinates repx, repy, repz
    repx, repy, repz should be in Box Size units
    '''
    vertices  = getVertices(repx, repy, repz)
    distances = np.sqrt( (vertices**2).sum(axis=1) )

    return vertices[distances.argmax()]

def getClosestIntersection(repx,repy,repz):
    '''
    Returns the Closest Intersection of the repetition
    with cordinates repx, repy, repz with the Past-Light Cone
    repx, repy, repz should be in Box Size units
    '''

    if (repx == 0) & (repy == 0) & (repz == 0):

        return 0.

    rep = np.array([repx, repy, repz])
    if np.sum(rep == 0.0) == 1:

        return np.sqrt( (getClosestAresta(repx,repy,repz)**2).sum() )

    elif np.sum( rep == 0.0 ) == 2:

        return np.sqrt( (getClosestFace(repx,repy,repz)**2).sum() )

    else:

        return np.sqrt( (getClosestVertice(repx,repy,repz)**2).sum() )


def getFarthestIntersection(repx,repy,repz):
    '''
    Returns the Farthest Intersection of the repetition
    with cordinates repx, repy, repz with the Past-Light Cone
    repx, repy, repz should be in Box Size units
    '''
    return  np.sqrt( (getFarthestVertice(repx,repy,repz)**2).sum() )


repmax = int( np.ceil( (cosmology.cosmo.comoving_distance(params.zmax).value - params.boxsize/2 )/params.boxsize) )

geometry = np.array(
     [ [repx, repy, repz, getClosestIntersection(repx,repy,repz), getFarthestIntersection(repx, repy, repz) ]
              for repx in range(-repmax, repmax+1) for repy in range(-repmax, repmax+1) for repz in range(-repmax, repmax+1)
     ])

if params.fovindeg < 180.0:

    repinfov = np.zeros(len(geometry), dtype=bool)

    for i, rep in enumerate(geometry[:,0:3]):

        for point in _surface + rep:

            theta = np.arccos(point[2]/np.linalg.norm(point))
            if theta < 0:
                theta += np.pi/2.0

            if theta <= params.fovinradians:

                repinfov[i] = True
                break

            else:

                repinfov[i] = False

    geometry = geometry[repinfov]

dt = np.dtype([('x', float), ('y', float),('z', float),('nearestpoint', float),('farthestpoint', np.float)])
geometry.dtype = dt
geometry['nearestpoint']  *= params.boxsize
geometry['farthestpoint'] *= params.boxsize

del _vertices, _arestas, _faces, _border, _surface
