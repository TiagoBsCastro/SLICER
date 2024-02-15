import numpy as np
import pyfftw
from scipy.stats import binned_statistic
from astropy.io import fits
from derivatives import laplacian_O3, gradientO4
from copy import deepcopy

##################################################  Utils  ########################################################
# Laplacian in 2D
def statistic_fft(phi, KX, KY, threads=1, input_space='real', statistic=None):

    if (statistic is None) or (statistic not in ['kappa', 'shear1', 'shear2']):

        raise RuntimeError('Statistic should be `kappa`, `shear1`, or `shear2`')

    # Prepare the FFT object for forward and inverse FFT
    if input_space=='real':

        # Allocate aligned arrays for output, and the inverse FFT output
        fft_input    = pyfftw.empty_aligned((phi.shape[0], phi.shape[1]), dtype='float32')
        fft_input[:] = phi
        fft_output   = pyfftw.empty_aligned((phi.shape[0], phi.shape[1]//2 + 1), dtype='complex64')
        fft_object   = pyfftw.FFTW(fft_input, fft_output, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_ESTIMATE',), threads=threads)

        # Perform the forward FFT
        phi_k = fft_object() 

    elif input_space=='fourier':

        # Allocate aligned arrays for output, and the inverse FFT output
        fft_input     = pyfftw.empty_aligned((phi.shape[0], 2*(phi.shape[1] - 1)), dtype='float32')
        fft_output    = pyfftw.empty_aligned(phi.shape, dtype='complex64')
        fft_output[:] = phi
        phi_k         = fft_output[:]

    else:

        raise RuntimeError(f'Unknown input space {input_space}. It should be either real or fourier.')

    if statistic == 'kappa':

        phi_k *= (KX**2 + KY**2) / 2.0

    elif statistic == 'shear1':

        phi_k *= (KX**2 - KY**2) / 2.0

    elif statistic == 'shear2':

        phi_k *= KX * KY

    ifft_object = pyfftw.FFTW(fft_output, fft_input, axes=(0, 1), direction='FFTW_BACKWARD', flags=('FFTW_ESTIMATE',), threads=threads)
    # Perform the inverse FFT
    return  -ifft_object()

def lensing_potential (kappa, KX, KY, threads=1, return_phi_k=False):

    # Allocate aligned arrays for output, and the inverse FFT output
    fft_input  = pyfftw.empty_aligned(kappa.shape, dtype='float32')
    fft_output = pyfftw.empty_aligned((kappa.shape[0], kappa.shape[1]//2 + 1), dtype='complex64')
    fft_input[:] = kappa
    
    # Prepare the FFT object for forward and inverse FFT
    fft_object  = pyfftw.FFTW(fft_input, fft_output, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_ESTIMATE',), threads=threads)
    ifft_object = pyfftw.FFTW(fft_output, fft_input, axes=(0, 1), direction='FFTW_BACKWARD', flags=('FFTW_ESTIMATE',), threads=threads)

    # Perform the forward FFT
    phi_k = fft_object()

    # Perform operations in Fourier space
    with np.errstate(divide='ignore', invalid='ignore'):
        phi_k /= (KX**2 + KY**2)
    
    # Handling the monopole term explicitly
    phi_k[0,0] = 0.0

    if return_phi_k:
        # This is needed to break the alias between phi_k and fft_output
        # that will be destroyed during the inverse FFT.
        phi_k = deepcopy(phi_k)

    # Inverse FFT
    ifft_object()
    # Perform the inverse FFT
    if return_phi_k:
        return -2 * ifft_object.output_array, -2 * phi_k
    else:
        return -2 * ifft_object.output_array

def unpad (array, n):

    return array[array.shape[0]//(2*n):-array.shape[0]//(2*n), array.shape[1]//(2*n):-array.shape[1]//(2*n)]

def PS (field, size, threads=1):

    from pyfftw.interfaces import numpy_fft as fft

    # Taking the fourier transform
    kx = 2 * np.pi * np.fft.fftfreq(field.shape[0], size/field.shape[0])
    ky = 2 * np.pi * np.fft.rfftfreq(field.shape[1], size/field.shape[1])
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    K = np.sqrt(KX**2 + KY**2).flatten()

    kf   = np.min([ KX[1,0], KY[0,1] ])
    bins = np.array([i*kf for i in range(field.shape[0])])

    P = fft.rfft2( field , threads=threads); P = P.flatten()
    P = P ** 2

    P = binned_statistic(K, P, bins=bins, statistic='mean').statistic
    K = binned_statistic(K, K, bins=bins, statistic='mean').statistic

    # Normalization
    bin_width1 = 2*np.pi/KX[1,0]/KX.shape[0]
    bin_width2 = 2*np.pi/KY[0,1]/KY.shape[1]
    N = bin_width1*bin_width2/KX.shape[0]/KX.shape[1]

    return np.transpose( [K, N*P] )

####################################################################################################################

def smr(fname, zero_padding=False, fout=None, derivative="FFT", threads=1, overwrite=False):
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

    threads: int (optional)
        Number of threads to be used by pyfftw

    overwrite: bool (optional)
        Whether overwite files in case the fout name already exists.

    returns:

    kappa: array
        Returns the convergence map computed from the reconstructed
        lensing potential.

    gamma, gamma1, gamma2 maps are saved directly on fout.
    '''
    # Getting the convergence map info
    nx   = fits.getheader(fname)['NAXIS1']
    ny   = fits.getheader(fname)['NAXIS2']
    size = fits.getheader(fname)['ANGLE'] * np.pi / 180.0

    # Zero padding the Map in place
    if zero_padding:

        size *= 2
        nx   *= 2
        ny   *= 2
        kappa = fits.getdata(fname)

        if (nx % 2) and (ny % 2):

            kappa = np.pad(kappa, ( (nx//2, nx//2+1), (ny//2, ny//2+1) ) )

        elif (nx % 2) and not (ny % 2):

            kappa = np.pad(kappa, ( (nx//2, nx//2+1), (ny//2, ny//2) ) )

        elif not (nx % 2) and (ny % 2):

            kappa = np.pad(kappa, ( (nx//2, nx//2), (ny//2, ny//2+1) ) )

        elif not (nx % 2) and not (ny % 2):

            kappa = np.pad(kappa, ( (nx//2, nx//2), (ny//2, ny//2) ) )

    else:

        kappa = fits.getdata(fname)

    # Setting the FFT for inverting the Laplacian
    kx = 2 * np.pi * np.fft.fftfreq(nx, size/nx).astype(np.float32)
    ky = 2 * np.pi * np.fft.rfftfreq(ny, size/ny).astype(np.float32)
    KX, KY = np.meshgrid(kx, ky, indexing='ij')

    # Getting the potential
    potential, phi_k = lensing_potential( kappa, KX, KY, threads=threads, return_phi_k=True)

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
        kappa  = statistic_fft(phi_k, KX, KY, threads=threads, input_space='fourier', statistic='kappa')
        gamma1 = statistic_fft(phi_k, KX, KY, threads=threads, input_space='fourier', statistic='shear1')
        gamma2 = statistic_fft(phi_k, KX, KY, threads=threads, input_space='fourier', statistic='shear2')

    else:

        raise NotImplementedError("derivative method '{}' is not implemented!".format(derivative))

    if zero_padding:

        potential = unpad(potential, 2)
        kappa     = unpad(kappa, 2); kappa -= kappa.mean()
        gamma1    = unpad(gamma1, 2); gamma1 -= gamma1.mean()
        gamma2    = unpad(gamma2, 2); gamma2 -= gamma2.mean()
        size     /= 2.0

    header = fits.getheader(fname)
    gamma  = np.sqrt(gamma1**2 + gamma2**2);

    if fout is None:

        # Checking if it was run by PowerBornApp or if it has kappa in the name
        fout = fname.replace(".fits", ".txt")
        fout = fname.replace("kappaBApp", "gamma1") if "kappaBApp" in fname else fname.replace("kappa", "gamma1")
        if fout == fname:

            raise ValueError("fout can not be None if 'kappa' or 'kappaBApp' is not in fname.")

        fits.writeto(fout, gamma1.astype(np.float32), header=header, overwrite=overwrite)
        fout = fname.replace("kappaBApp", "gamma2") if "kappaBApp" in fname else fname.replace("kappa", "gamma2")
        fits.writeto(fout, gamma2.astype(np.float32), header=header, overwrite=overwrite)
        fout = fname.replace("kappaBApp", "gamma") if "kappaBApp" in fname else fname.replace("kappa", "gamma")
        fits.writeto(fout, gamma.astype(np.float32), header=header, overwrite=overwrite)
        fout = fname.replace("kappaBApp", "phi") if "kappaBApp" in fname else fname.replace("kappa", "phi")
        fits.writeto(fout, potential.astype(np.float32), header=header, overwrite=overwrite)

    else:

        try:

            # Checking if it was run by PowerBornApp or if it has kappa in the name
            fits.writeto(fout.format(statistic="gamma1"), gamma1.astype(np.float32), header=header, overwrite=overwrite)
            fits.writeto(fout.format(statistic="gamma2"), gamma2.astype(np.float32), header=header, overwrite=overwrite)
            fits.writeto(fout.format(statistic="gamma"), gamma.astype(np.float32), header=header, overwrite=overwrite)
            fits.writeto(fout.format(statistic="phi"), potential.astype(np.float32), header=header, overwrite=overwrite)

        except KeyError:

            raise ValueError("fout should have a single formater '{statistic}'!")

    return kappa
