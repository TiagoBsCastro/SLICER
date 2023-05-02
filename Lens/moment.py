import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from scipy.stats import moment

NMOMENT = 8

# Ref.
files = glob.glob("ref/*/lens_maps/maps_1024_*/gadget_z1.1137_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

# Om02
files = glob.glob("Om02/lens_maps/maps_1024_*/gadget_z1.0162_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

# Om04
files = glob.glob("Om04/lens_maps/maps_1024_*/gadget_z1.0398_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

# h06
files = glob.glob("h06/lens_maps/maps_1024_*/gadget_z1.1137_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

# h08
files = glob.glob("h08/lens_maps/maps_1024_*/gadget_z1.1137_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

# s807
files = glob.glob("s807/lens_maps/maps_1024_*/gadget_z1.1137_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

# s809
files = glob.glob("s809/lens_maps/maps_1024_*/gadget_z1.1137_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

quit()
# I have to fix SLICER to work with wCDM
# w-07
files = glob.glob("w-07/lens_maps/maps_1024_*/gadget_z1.056_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

# w-1p3
files = glob.glob("w-1p3/lens_maps/maps_1024_*/gadget_z1.0349_kappa.fits")
mean  = np.mean([fits.getdata(f).mean() for f in files])
m = np.zeros(NMOMENT-1)
for f in files:

    kappa = fits.getdata(f)
    for i in range(2, NMOMENT+1):

        m[i-2] += np.sum( (kappa - mean)**i )

m = m/len(files)/kappa.size
print(m)

