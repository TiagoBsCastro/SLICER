import numpy as np
from scipy.stats import binned_statistic
from astropy.io import fits
from derivatives import laplacian_O3, gradientO4

##################################################  Utils  ########################################################
# Laplacian in 2D
def convergence_fft(potential, KX, KY):

    return - np.fft.ifft2( np.fft.fft2(potential) * (KX**2+KY**2)/2.0, potential.shape ).real

# Shear direct terms
def shear1_fft(potential, KX, KY):

    return - np.fft.ifft2( np.fft.fft2(potential) * (KX**2-KY**2)/2.0, potential.shape ).real

# Shear cross terms
def shear2_fft(potential, KX, KY):

    return - np.fft.ifft2( KX*KY*np.fft.fft2(potential), potential.shape).real

def lensing_potential (kappa, KX, KY):

    # Going to Fourier Space
    with np.errstate(divide='ignore', invalid='ignore'):

        phi_k = 2.0 * np.fft.fft2(kappa) / (KX**2+KY**2)

    # Changing the monopole
    phi_k[0,0] = 0.0

    return - np.fft.ifft2( phi_k, kappa.shape ).real

def unpad (array, n):

    return array[array.shape[0]//(2*n):-array.shape[0]//(2*n), array.shape[1]//(2*n):-array.shape[1]//(2*n)]

def PS (field, KX, KY, n=50):

    # Taking the fourier transform
    K  = np.sqrt(KX**2 + KY**2).flatten()
    P  = np.fft.fft2(field).flatten()
    P *= P.conj(); P = P.real

    # Binning; taking the third bin due to 0-padding
    kmin = np.min([ KX[2,0], KY[0,2] ])
    kmax = np.max([ KX[1:, 0].max(), KY[0, 1:].max() ])
    bins = np.geomspace(kmin, kmax, n)

    P = binned_statistic(K, P, bins=bins, statistic='mean').statistic
    K = binned_statistic(K, K, bins=bins, statistic='mean').statistic

    # Normalization
    bin_width1 = 2*np.pi/KX[1,0]/KX.shape[0]
    bin_width2 = 2*np.pi/KY[0,1]/KY.shape[1]
    N = bin_width1*bin_width2/KX.shape[0]/KX.shape[1]

    return np.transpose( [K, N*P] )

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
    size  = fits.getheader(fname)['ANGLE'] * np.pi / 180.0

    # Zero padding the Map in place
    if (kappa.shape[0] % 2) and (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2+1), (kappa.shape[1]//2, kappa.shape[1]//2+1) ) )

    elif (kappa.shape[0] % 2) and not (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2+1), (kappa.shape[1]//2, kappa.shape[1]//2) ) )

    elif not (kappa.shape[0] % 2) and (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2), (kappa.shape[1]//2, kappa.shape[1]//2+1) ) )

    elif not (kappa.shape[0] % 2) and not (kappa.shape[1] % 2):

        kappa = np.pad(kappa, ( (kappa.shape[0]//2, kappa.shape[0]//2), (kappa.shape[1]//2, kappa.shape[1]//2) ) )

    size *= 2

    # Setting the FFT for inverting the Laplacian
    kx = 2 * np.pi * np.fft.fftfreq(kappa.shape[0], size/kappa.shape[0])
    ky = 2 * np.pi * np.fft.fftfreq(kappa.shape[1], size/kappa.shape[1])
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
        gamma2    = 1/2 * ( d2_xy + d2_yx )

    elif derivative == "FFT":

        # Recomputing kappa (for redundancy check), gamma1, gamma2
        kappa  = convergence_fft(potential, KX, KY)
        gamma1 = shear1_fft(potential, KX, KY)
        gamma2 = shear2_fft(potential, KX, KY)

    else:

        raise NotImplementedError("derivative method '{}' is not implemented!".format(derivative))

    Pk  = PS(kappa, KX, KY); Nk = kappa.var()
    Pg1 = PS(gamma1, KX, KY); Ng1 = gamma1.var()
    Pg2 = PS(gamma2, KX, KY); Ng2 = gamma2.var()

    potential = unpad(potential, 2)
    kappa     = unpad(kappa, 2); kappa -= kappa.mean(); Nk /= kappa.var()
    gamma1    = unpad(gamma1, 2); gamma1 -= gamma1.mean(); Ng1 /= gamma1.var()
    gamma2    = unpad(gamma2, 2); gamma2 -= gamma2.mean(); Ng2 /= gamma2.var()
    gamma     = np.sqrt(gamma1**2 + gamma2**2);

    Pk[:,1] /= Nk; Pg1[:,1] /= Ng1; Pg2[:,1] /= Ng2;
    Pg  = np.array( [Pg1[:,0], Pg1[:, 1]+Pg2[:, 1]] )

    header = fits.getheader(fname)

    if fout is None:

        # Checking if it was run by PowerBornApp or if it has kappa in the name
        fout = fname.replace(".fits", ".txt")
        np.savetxt(fout, Pk)
        fout = fname.replace("kappaBApp", "gamma1") if "kappaBApp" in fname else fname.replace("kappa", "gamma1")
        if fout == fname:

            raise ValueError("fout can not be None if 'kappa' or 'kappaBApp' is not in fname.")

        fits.writeto(fout, gamma1.astype(np.float32), header=header)
        fout = fout.replace(".fits", ".txt")
        np.savetxt(fout, Pg1)
        fout = fname.replace("kappaBApp", "gamma2") if "kappaBApp" in fname else fname.replace("kappa", "gamma2")
        fits.writeto(fout, gamma2.astype(np.float32), header=header)
        fout = fout.replace(".fits", ".txt")
        np.savetxt(fout, Pg2)
        fout = fname.replace("kappaBApp", "gamma") if "kappaBApp" in fname else fname.replace("kappa", "gamma")
        fits.writeto(fout, gamma.astype(np.float32), header=header)
        fout = fout.replace(".fits", ".txt")
        np.savetxt(fout, Pg)
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
