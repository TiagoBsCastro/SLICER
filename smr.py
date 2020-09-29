import numpy as np
from numpy.fft import ifft2, fft2, fft, ifft, fftfreq, fftshift, ifftshift
from astropy.io import fits
from derivatives import laplacian_O3, gradientO4

##################################################  Utils  ########################################################

# Laplacian in 2D
def convergence_fft(potential, KX, KY):

    return - (ifft(KX**2*fft(potential/2, axis = 0), axis = 0) + ifft(KY**2*fft(potential/2, axis = 1), axis = 1) )

# Shear direct terms
def shear1_fft(potential, KX, KY):

    return - (ifft(KX**2*fft(potential/2, axis = 0), axis = 0) - ifft(KY**2*fft(potential/2, axis = 1), axis = 1) )

# Shear cross terms
def shear2_fft(potential, KX, KY):

    return - ifft2( KX*KY*fft2(potential) )

def lensing_potential (kappa, KX, KY):

    return ifft2( -2*fft2(kappa)/(KX**2 + KY**2) )

def unpad (array, n):

    return array[array.shape[0]//(2*n):-array.shape[0]//(2*n), array.shape[1]//(2*n):-array.shape[1]//(2*n)]

####################################################################################################################

def smr(fname, fout=None, derivative="FFT"):
    '''
    Computes the shear maps reconstructing the lensing potential
    from the convergence map.

    fname: str
        Convergence map file name.
    
    fout: str
        File name for the output map. 
        It must contain a single formater '{statistic}'.
    
    derivative: str (optional)
        Method for the computation of the potential derivatives.

            * FFT (Default): Compute the derivatives using FFT and IFFT.
            * gradient     : Use forth order finite differences.
    
    returns:

    kappa: array
        Returns the convergence map computed from the reconstructed
        lensing potential.

    gamma, gamma1, gamma2 maps are saved directly on fout.
    '''
    # Getting the convergence map
    kappa = fits.getdata(fname)
    # Zero padding the Map in place
    kappa = np.pad(kappa, (kappa.shape[0]//2, kappa.shape[1]//2))

    # Setting the FFT for inverting the Laplacian
    kx     = 2 * np.pi * fftfreq(kappa.shape[0], 1)
    ky     = 2 * np.pi * fftfreq(kappa.shape[1], 1)
    # Setting the zeroth mode to infinity to avoid divergence
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    KX[kx == 0, ky == 0] = np.inf
    KY[kx == 0, ky == 0] = np.inf

    # Getting the potential
    potential = lensing_potential( kappa, KX, KY).real

    if derivative == "gradient":

        # Computing derivatives
        d1_x, d1_y   = gradientO4(potential)
        d2_xx, d2_xy = gradientO4( d1_x )
        d2_yx, d2_yy = gradientO4( d1_y )
        d2_xx, d2_yy  = laplacian_O3(potential)

        # Recomputing kappa (for redundancy check), gamma1, gamma2
        kappa     = 1/2 * ( d2_xx + d2_yy )
        kappa     = unpad(kappa, 2)
        gamma1    = 1/2 * ( d2_xx - d2_yy )
        gamma2    = d2_xy
        gamma     = np.sqrt(gamma1**2 + gamma2**2)
        potential = unpad(potential, 2)

    elif derivative == "FFT":

        # Recomputing kappa (for redundancy check), gamma1, gamma2
        KX, KY = np.meshgrid(kx, ky, indexing='ij')
        kappa  = convergence_fft(potential, KX, KY).real
        kappa  = unpad(kappa, 2)
        gamma1 = shear1_fft(potential, KX, KY).real
        gamma1 = unpad(gamma1, 2)
        gamma2 = shear2_fft(potential, KX, KY).real
        gamma2 = unpad(gamma2, 2)
        gamma  = np.sqrt(gamma1**2 + gamma2**2)
        potential = unpad(potential, 2)

    else:

        raise NotImplementedError("derivative method '{}' is not implemented!".format(derivative))

    header = fits.getheader(fname)

    if fout is None:

        # Checking if it was run by PowerBornApp or if it has kappa in the name
        fout = fname.replace("kappaBApp", "gamma1") if "kappaBApp" in fname else fname.replace("kappa", "gamma1")
        if fout == fname:

            raise ValueError("fout can not be None if 'kappa' or 'kappaBApp' is not in fname.")

        fits.writeto(fout, gamma1.astype(np.float32), header=header)
        fout = fname.replace("kappaBApp", "gamma2") if "kappaBApp" in fname else fname.replace("kappa", "gamma2")
        fits.writeto(fout, gamma2.astype(np.float32), header=header)
        fout = fname.replace("kappaBApp", "gamma") if "kappaBApp" in fname else fname.replace("kappa", "gamma")
        fits.writeto(fout, gamma.astype(np.float32), header=header)
        fout = fname.replace("kappaBApp", "phi") if "kappaBApp" in fname else fname.replace("kappa", "phi")
        fits.writeto(fout, potential.astype(np.float32), header=header)

    else:

        try:

            # Checking if it was run by PowerBornApp or if it has kappa in the name
            fits.writeto(fout.format(statistic="gamma1"), gamma1.astype(np.float32), header=header)
            fits.writeto(fout.format(statistic="gamma2"), gamma2.astype(np.float32), header=header)
            fits.writeto(fout.format(statistic="gamma"), gamma.astype(np.float32), header=header)
            fits.writeto(fout.format(statistic="phi"), potential.astype(np.float32), header=header)

        except KeyError:

            raise ValueError("fout should have a single formater '{statistic}'!")

    return kappa
