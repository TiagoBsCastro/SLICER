import sys
from string import Template

################################ Lensing parameters #####################################

# Map index
#imap = 9
imap = sys.argv[1]
# The index of the source redshifts
source_redshifts = 'All'

################################ Cosmology Pars #########################################
Om0    = 0.3000
H0     = 70.000
ns     = 0.9661
sigma8 = 0.9000
Ob0    = 0.0494
w0     = -1.300
wa     = 0.0000

pert_params     =  {'flat': True, 'H0': H0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': sigma8, 'ns': ns, 'relspecies':False, 'de_model':'w0wa', 'w0':w0, 'wa':wa}

################################ Base Directory #########################################
baseDirectory   = "w-1p3/"
############################# Mass Maps Directory #######################################
massMapsNames   =  Template(baseDirectory+'density_maps/maps_1024_${mapIndex}/gadget.*fits')
############################# Kappa Maps Directory ######################################
kappaMapsNames  = Template(baseDirectory+'lens_maps/maps_1024_${mapIndex}/gadget_z${zs}_kappa.fits')
############################# Delta Maps Directory ######################################
deltaMapsNames  = Template(baseDirectory+'lens_maps/maps_1024_${mapIndex}/gadget_z${zs}_delta.fits')
################################# Planes Files ##########################################
planesFileNames = Template(baseDirectory+'density_maps/maps_1024_${mapIndex}/planes_list_${mapIndex}.txt')
