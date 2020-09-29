from astropy.cosmology import FlatLambdaCDM
from string import Template

################################ Lensing parameters #####################################

# Map index
imap = 9
# The index of the source redshifts
source_redshifts = [15, 23, 29, 36]

################################ Cosmology Pars #########################################
Om0    = 0.301
H0     = 68.2
ns     = 0.973
sigma8 = 0.799
Ob0    = 0.048

pert_params     =  {'flat': True, 'H0': H0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': sigma8, 'ns': ns, 'relspecies':False}

################################ Base Directory #########################################
baseDirectory   = "/beegfs/tcastro/Alice/hr/"
############################# Mass Maps Directory #######################################
massMapsNames   =  Template(baseDirectory+'density_maps/maps_512_${mapIndex}/gadget.*fits')
############################# Kappa Maps Directory ######################################
kappaMapsNames  = Template(baseDirectory+'lens_maps/maps_512_${mapIndex}/gadget_z${zs}_kappa.fits')
############################# Delta Maps Directory ######################################
deltaMapsNames  = Template(baseDirectory+'lens_maps/maps_512_${mapIndex}/gadget_z${zs}_delta.fits')
################################# Planes Files ##########################################
planesFileNames = Template(baseDirectory+'density_maps/maps_512_${mapIndex}/planes_list_${mapIndex}.txt')
