import numpy as np

class OddLensPlanes(Exception):
   """
   Base class for other exceptions
   """
   pass

def planesSanityCheck (dllow, dlup):
    """
    Checks the sanity of the lens planes:

    dllow: Array
        Numpy array with the comoving distance lower boundary of the plane
    dlup: Array
        Numpy array with the upper comoving distance boundary of the plane
    """
    size = dllow[0].size
   
    if size != dlup[0].size:

        raise OddLensPlanes("!DL_ and DL^ are not the same length!")

    for i in range(size-1):

        if all(dlup[i] != dllow[i+1]):

            raise OddLensPlanes("!There are holes in the Lens Planes!")

def genericHeader (header, naxis, source, angle):
    """
    Generic fits file header
    """
    header['NAXIS'] = 2
    header['NAXIS1'] = naxis
    header['NAXIS2'] = naxis
    header['ZSOURCE'] = source
    header['ANGLE'] = angle

def A (z, theta):
    """
    Returns the area of the squared-based cone

    z: float
        Height of the cone
    theta:
        Angular aperture
    """
    return (2 * z * np.tan( theta * np.pi / 180.0 / 2.0) )**2

def V (z, theta):
    """
    Returns the volume of the squared-based cone

    z: float
        Height of the cone
    theta:
        Angular aperture
    """
    
    return z/3.0 * A(z, theta)

def lensKernel (z, zs, cosmo):
    """
    Returns the so-called lensing efficiency.

    z: Float
        Redshift of the plane
    zs: Float
        Redshift of the source
    cosmo: Astropy cosmology instance
        Cosmology object to compute distances
    """

    return (cosmo.angular_diameter_distance(z) * cosmo.angular_diameter_distance_z1z2(z,zs) / cosmo.angular_diameter_distance(zs)).value

def lensKernelIntegral (z1, z2, zs, cosmo, f=None):
    """
    Returns the integral of the lensKernernel weighted by f.

    z1: Float
        Lower limit
    z2: Float
        Upper limit
    zs: Astropy cosmology instance
        Cosmology object to compute distances
    f: function (default is None)
        Wheight function
    """

    from scipy.integrate import quad

    xs = cosmo.angular_diameter_distance(zs) 

    if f is None:

        integrand = lambda z: (cosmo.angular_diameter_distance(z) * cosmo.angular_diameter_distance_z1z2(z,zs) \
                        / xs * cosmo.hubble_distance / cosmo.efunc(z)).value

    elif callable(f):

        integrand = lambda z: (cosmo.angular_diameter_distance(z) * cosmo.angular_diameter_distance_z1z2(z,zs) \
                        / xs * cosmo.hubble_distance.value / cosmo.efunc(z) * f(z)).value

    else:

        raise TypeError("The wheight function f is not calable")

    return quad(integrand, z1, z2)[0]
