import numpy as np

#4th order accurate gradient function based on 2nd order version from http://projects.scipy.org/scipy/numpy/browser/trunk/numpy/lib/function_base.py
def gradientO4(f, *varargs):

    """Calculate the fourth-order-accurate gradient of an N-dimensional scalar function.

    Uses central differences on the interior and first differences on boundaries

    to give the same shape.

    Inputs:

      f -- An N-dimensional array giving samples of a scalar function

      varargs -- 0, 1, or N scalars giving the sample distances in each direction

    Outputs:

      N arrays of the same shape as f giving the derivative of f with respect

       to each dimension.

    """

    N = len(f.shape)  # number of dimensions

    n = len(varargs)

    if n == 0:

        dx = [1.0]*N

    elif n == 1:

        dx = [varargs[0]]*N

    elif n == N:

        dx = list(varargs)

    else:

        raise SyntaxError( "invalid number of arguments" )

    # use central differences on interior and first differences on endpoints

    #print dx

    outvals = []

    # create slice objects --- initially all are [:, :, ..., :]

    slice0 = [slice(None)]*N

    slice1 = [slice(None)]*N

    slice2 = [slice(None)]*N

    slice3 = [slice(None)]*N

    slice4 = [slice(None)]*N

    otype = f.dtype.char

    if otype not in ['f', 'd', 'F', 'D']:

        otype = 'd'

    for axis in range(N):       

        # select out appropriate parts for this dimension

        out = np.zeros(f.shape, f.dtype.char)

        slice0[axis] = slice(2, -2)

        slice1[axis] = slice(None, -4)

        slice2[axis] = slice(1, -3)

        slice3[axis] = slice(3, -1)

        slice4[axis] = slice(4, None)

        # 1D equivalent -- out[2:-2] = (f[:4] - 8*f[1:-3] + 8*f[3:-1] - f[4:])/12.0

        out[slice0] = (f[slice1] - 8.0*f[slice2] + 8.0*f[slice3] - f[slice4])/12.0

        slice0[axis] = slice(None, 2)

        slice1[axis] = slice(1, 3)

        slice2[axis] = slice(None, 2)

        # 1D equivalent -- out[0:2] = (f[1:3] - f[0:2])

        out[slice0] = (f[slice1] - f[slice2])

        slice0[axis] = slice(-2, None)

        slice1[axis] = slice(-2, None)

        slice2[axis] = slice(-3, -1)

        ## 1D equivalent -- out[-2:] = (f[-2:] - f[-3:-1])

        out[slice0] = (f[slice1] - f[slice2])

        # divide by step size

        outvals.append(out / dx[axis])

        # reset the slice object in this dimension to ":"

        slice0[axis] = slice(None)

        slice1[axis] = slice(None)

        slice2[axis] = slice(None)

        slice3[axis] = slice(None)

        slice4[axis] = slice(None)

    if N == 1:

        return outvals[0]

    else:

        return outvals       

def laplacian_O3(f, *varargs):

    """ Third order accurate, 5-point formula. Not really laplacian, but second derivative along each axis. Sum the outputs to get laplacian."""

    N = len(f.shape)  # number of dimensions

    n = len(varargs)

    if n == 0:

        dx = [1.0]*N

    elif n == 1:

        dx = [varargs[0]]*N

    elif n == N:

        dx = list(varargs)

    else:

        raise SyntaxError( "invalid number of arguments" )

    # use central differences on interior and first differences on endpoints

    outvals = []

    # create slice objects --- initially all are [:, :, ..., :]

    slice0 = [slice(None)]*N

    slice1 = [slice(None)]*N

    slice2 = [slice(None)]*N

    slice3 = [slice(None)]*N

    slice4 = [slice(None)]*N    

    otype = f.dtype.char

    if otype not in ['f', 'd', 'F', 'D']:

        otype = 'd'

    for axis in range(N):       

        # select out appropriate parts for this dimension

        out = np.zeros(f.shape, f.dtype.char)

        # http://www.sitmo.com/eq/262 (3rd order accurate)

        slice0[axis] = slice(2, -2)

        slice1[axis] = slice(None, -4)

        slice2[axis] = slice(1, -3)

        slice3[axis] = slice(3, -1)

        slice4[axis] = slice(4, None)

        # 1D equivalent -- out[2:-2] = (-f[:-4] + 16*f[1:-3] + -30*f[2:-2] + 16*f[3:-1] - f[4:])/12.0

        out[slice0] = (-f[slice1] + 16.0*f[slice2] - 30.0*f[slice0] + 16.0*f[slice3] - f[slice4])/12.0

        # http://www.sitmo.com/eq/260 (2nd order accurate; there's also a 3rd order accurate that requires 5 points) 

        slice0[axis] = slice(None, 2)

        slice1[axis] = slice(3, 5)

        slice2[axis] = slice(2, 4)

        slice3[axis] = slice(1, 3)

        # 1D equivalent -- out[0:2] = 2*f[0:2] - 5*f[1:3] + 4*f[2:4] - f[3:5]

        out[slice0] = (2.0*f[slice0] - 5.0*f[slice3] + 4.0*f[slice2] - f[slice1])

        slice0[axis] = slice(-2, None)

        slice1[axis] = slice(-5, -3)

        slice2[axis] = slice(-4, -2)

        slice3[axis] = slice(-3, -1)

        # 1D equivalent -- out[0:2] = 2*f[0:2] - 5*f[1:3] + 4*f[2:4] - f[3:5]

        out[slice0] = (2.0*f[slice0] - 5.0*f[slice3] + 4.0*f[slice2] - f[slice1])

        # divide by step size

        axis_dx = dx[axis]

        outvals.append(out / (axis_dx*axis_dx))

        # reset the slice object in this dimension to ":"

        slice0[axis] = slice(None)

        slice1[axis] = slice(None)

        slice2[axis] = slice(None)

        slice3[axis] = slice(None)

        slice4[axis] = slice(None)

    if N == 1:

        return outvals[0]

    else:

        return outvals

def test_laplacian_03():

    """ Simple sanity check on whether the centers and edges have signs correct with the right value"""

    increasing = np.arange(10)-5

    decreasing = -increasing

    ones = np.ones((10,))

    x_grid = ones[:9]*increasing[:,None]

    y_grid = ones[:, None]*increasing[None,:9]

    grid = x_grid + y_grid

    d2dx2, d2dy2 = laplacian_O3(grid**3.0, 1, 1)    

    assert (d2dx2 == 6.0*grid).all()

    assert (d2dy2 == 6.0*grid).all()

    d2dx2, d2dy2 = laplacian_O3(grid**2.0, 3.0, 5.0)    

    assert (d2dx2 == 2.0/3.0/3.0).all()

    assert (d2dy2 == 2.0/5.0/5.0).all()

    print("Laplacian checks ok")
