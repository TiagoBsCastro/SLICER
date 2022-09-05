# This script creates the InputParams.ini for SLICER
from string import Template
from random import seed
from random import randint
import sys

# Setting the random seed
seed(3395)
run   = sys.argv[1]
index = int(sys.argv[2])
r1 = randint(1, 9999) * index
r2 = randint(1, 9999) * index
r3 = randint(1, 9999) * index

params = Template("""##### 1. Number of Map Pixels ##
1024
##### 2. Source Redshift #######
1.0
##### 3. Field of View #########
2.0
##### 4. File with Snapshots ###
snaps.txt
##### 5. Snapshots Directory ###
${run}
##### 6. PLC Sim. Name #########
gadget
##### 7. Seed for Pos. Center ##
$r1
##### 8. Seed for Pos. Reflec. #
$r2
##### 9. Seed for Axis Sel. ####
$r3
##### 10. Part. in Planes ######
0
##### 11. PLC Directory ########
${run}/density_maps/maps_1024_${index}/
##### 12. PLC Suffix ###########
${index}
##### 13. Part. Degradation ####
0""")

print(params.substitute(run=run, index=index, r1=r1, r2=r2, r3=r3))
