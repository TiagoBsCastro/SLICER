import numpy as np
from numpy.fft import fft, ifft, fft2, ifft2, fftfreq
from astropy.io import fits
from derivatives import laplacian_O3, gradientO4

##################################################  Utils  ########################################################

# Laplacian in 2D
def convergence_fft(potential, KX, KY, zero_mode=None):

    kappa_k = np.fft.fft2(potential) * (KX**2+KY**2)

    if zero_mode is None:

        kappa = - np.fft.ifft2( kappa_k, potential.shape ).real

        return kappa - kappa[0,0]

    else:

        kappa_k[0, 0] = - zero_mode

        return - np.fft.ifft2( kappa_k, potential.shape ).real

# Shear direct terms
def shear1_fft(potential, KX, KY, zero_mode=None):

    shear1_k = np.fft.fft2(potential) * (KX**2-KY**2)/2.0

    if zero_mode is None:

        gamma1 = - np.fft.ifft2( shear1_k, potential.shape ).real
        return gamma1 - gamma1[0,0]

    else:

        shear1_k[0, 0] = - zero_mode

        return - np.fft.ifft2( shear1_k, potential.shape ).real

# Shear cross terms
def shear2_fft(potential, KX, KY, zero_mode=None):

    shear2_k = KX*KY*fft2(potential)

    if zero_mode is None:

        gamma2 = - ifft2( shear2_k, potential.shape).real

        return gamma2 - gamma2[0,0]

    else:

        shear2_k[0,0] = - zero_mode

        return - ifft2( shear2_k, potential.shape).real

def lensing_potential (kappa, KX, KY):

    # Going to Fourier Space
    with np.errstate(divide='ignore', invalid='ignore'):

        phi_k = np.fft.fft2(kappa) / (KX**2+KY**2)

    # Changing the monopole
    phi_k[0,0] = 0.0

    return - np.fft.ifft2( phi_k, kappa.shape ).real

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
    if (kappa.shape[0] % 2) and (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2+1), (kappa.shape[1]//2, kappa.shape[1]//2+1) ) )

    elif (kappa.shape[0] % 2) and not (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2+1), (kappa.shape[1]//2, kappa.shape[1]//2) ) )

    elif not (kappa.shape[0] % 2) and (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2), (kappa.shape[1]//2, kappa.shape[1]//2+1) ) )

    elif not (kappa.shape[0] % 2) and not (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2), (kappa.shape[1]//2, kappa.shape[1]//2) ) )

    # Setting the FFT for inverting the Laplacian
    kx     = 2 * np.pi * fftfreq(kappa.shape[0], 1)
    ky     = 2 * np.pi * fftfreq(kappa.shape[1], 1)
    # Setting the zeroth mode to infinity to avoid divergence
    KX, KY = np.meshgrid(kx, ky, indexing='ij')

    # Getting the potential
    potential = lensing_potential( kappa, KX, KY )

    if derivative == "gradient":

        # Computing derivatives
        d1_x, d1_y   = gradientO4(potential)
        d2_xx, d2_xy = gradientO4( d1_x )
        d2_yx, d2_yy = gradientO4( d1_y )
        d2_xx, d2_yy = laplacian_O3(potential)

        # Recomputing kappa (for redundancy check), gamma1, gamma2
        kappa     = 1/2 * ( d2_xx + d2_yy )
        gamma1    = 1/2 * ( d2_xx - d2_yy )
        gamma2    = d2_xy
        gamma     = np.sqrt(gamma1**2 + gamma2**2)

    elif derivative == "FFT":

        # Recomputing kappa (for redundancy check), gamma1, gamma2
        kappa  = convergence_fft(potential, KX, KY)
        gamma1 = shear1_fft(potential, KX, KY)
        gamma2 = shear2_fft(potential, KX, KY)
        gamma  = np.sqrt(gamma1**2 + gamma2**2)

    else:

        raise NotImplementedError("derivative method '{}' is not implemented!".format(derivative))

    # Unpadding the maps
    potential = unpad(potential, 2)
    kappa     = unpad(kappa, 2)
    gamma1    = unpad(gamma1, 2)
    gamma2    = unpad(gamma2, 2)
    gamma     = unpad(gamma, 2)

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
