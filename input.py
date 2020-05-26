from astropy.cosmology import FlatLambdaCDM
from string import Template

imap = 9

################################ Cosmology Pars #########################################
Om0    = 0.272
H0     = 70.4
ns     = 0.963
sigma8 = 0.809
Ob0    = 0.168 * 0.272

pert_params     =  {'flat': True, 'H0': H0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': sigma8, 'ns': ns, 'relspecies':False}

############################# Mass Maps Directory #######################################
massMapsNames   =  Template('maps_625_${mapIndex}/gadget.*fits')
############################# Kappa Maps Directory ######################################
kappaMapsNames  = Template('maps_625_${mapIndex}/gadget_z${zs}_kappa.fits')
############################# Delta Maps Directory ######################################
deltaMapsNames  = Template('maps_625_${mapIndex}/gadget_z${zs}_delta.fits')
################################# Planes Files ##########################################
planesFileNames = Template('maps_625_${mapIndex}/planes_list_${mapIndex}.txt')
