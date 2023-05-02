import plot_style
import matplotlib.pyplot as plt
import numpy as np

################### INPUT #################
aplc = [0.731, 0.626, 0.576, 0.512, 0.452, 0.399, 0.337]
abox = [0.73631, 0.62708, 0.57870, 0.51717, 0.45849, 0.40322, 0.33523]
###########################################

for a1, a2 in zip(aplc, abox):

    box = np.loadtxt("hlist_{:.5f}.txt".format(a2)) 
    plc = np.load("piccolo-C0-hmf-a={:.3f}.npy".format(a1)) 
    plt.plot(box[:, 0], box[:, 5]/plc[:, 1], label="Box/PLC")
    plt.xscale('log') 
    plt.legend() 
    plt.xlabel(r"$M_{\rm Vir.}\,[\,M_\odot/h]$")
    plt.ylabel(r"${\rm Res.}\,dn/dlogM$")
    plt.show()
