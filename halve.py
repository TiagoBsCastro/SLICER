from astropy.io import fits
import numpy as np

fname = "gadget.016.plane_10000_0.fits"
res   = 10000

for i in range(4):

	orig   = fits.getdata(fname)
	header = fits.getheader(fname)

	npix = int(np.sqrt(orig.size)/2)

	degr = np.zeros((npix, npix), dtype=np.float32)


	for i in range(npix):

	    for j in range(npix):

                 degr[i,j]  = (orig[2*i,2*j] + orig[2*i+1, 2*j] + orig[2*i, 2*j+1] + orig[2*i+1,2*j+1])

	header['NAXIS1'] = header['NAXIS1']/2
	header['NAXIS2'] = header['NAXIS2']/2
	fname = fname.replace(".fits", "_halved.fits")
	fname = fname.replace(str(res), str(res/2))
        res  /= 2
	fits.writeto(fname, degr, header, overwrite=True)
        
