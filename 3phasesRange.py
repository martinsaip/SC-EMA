__author__ = 'Aydin & Titrian'

from tools.gibbs_plot import plotYoungsModulus
from tools.gibbs_plot import plot3phase

from polycrystal import cubic
from polycrystal import tetragonal
from polycrystal import trigonal
from polycrystal import hexagonal
from polycrystal import orthorombic
import numpy as np

import sys
import os
import numpy as np

tmpFolder = str(sys.argv[1])

# This works for Windows
file = open("../public/"+tmpFolder+"/3phase/crystalStack.txt","r")
crysStack = file.read()
crysList = crysStack.splitlines()
file.close()

crys1 = eval(crysList[0])
crys2 = eval(crysList[1])
crys3 = eval(crysList[2])

crystalName1 = crys1.crystalname
crystalName2 = crys2.crystalname
crystalName3 = crys3.crystalname
#conc1 = np.linspace(1,0,11)
#conc2 = np.linspace(0,1,11)



bulkList = []
youngsList = []
shearList = []
poissonList = []
kList = []
conc1 = []
#myks = []
myvals = []

steps = 10
for k in range(steps+1):
    conc1 = np.round(1 - k*1 / steps, 8)    # Thanks to Python 3, we do not need to explicitly state "float"
    for l in range(k+1):
        conc2 = np.round(1 - conc1 - l*1 / steps, 8)    # Thanks to Python 3, we do not need to explicitly state "float"
        conc3 = np.round(1 - conc1 - conc2, 8)
        kList.append([conc1, conc2, conc3])
#        print("concentrations of elements A,B,C: %s %s %s" % (conc1, conc2, conc3) + '<br>')

        #        myks.append([conc1, conc2, conc3])
        crys1.setConc( conc1 )
        crys2.setConc( conc2 )
        crys3.setConc( conc3 )

        my_crystal = crys1 + crys2 + crys3

        myDict = my_crystal.getElasticDict()

        my_crystal.setK0andMue0()
        my_crystal.setYoungsMod()
        my_crystal.setpoissonratio()

        bulkList.append(my_crystal.K0 * 100)
        shearList.append(my_crystal.mue0 * 100)
        youngsList.append(my_crystal.E * 100)
        poissonList.append(my_crystal.v)

#write to separate table
f = open("../public/"+tmpFolder+"/3phase/dataTable.txt","w")
f.write("#conc.Phase A | #conc.Phase B | #conc.Phase C | Bulk Mod (GPa) | Shear Mod (GPa) | Youngs Mod (GPa) | poiss. ratio\n")
#steps = 10
for k in range(len(kList)):
    line = '{0:.1f}\t\t{1:.1f}\t\t{2:.1f}\t\t{3:3.4f}          {4:3.4f}           {5:3.4f}            {6:3.4f}\n'.\
    format(float(kList[k][0]), float(kList[k][1]), float(kList[k][2]),
            float(bulkList[k]),float(shearList[k]),float(youngsList[k]), float(poissonList[k]))
    f.write(line)

f.close()

plot3phase(concList = kList, valList = bulkList, labelA = crystalName1, labelB = crystalName2, labelC = crystalName3,
    axeslabelx = 'Bulk Modulus (GPa)', saveFig = "True", type = "bulkMod", tmpFolder = tmpFolder)
plot3phase(concList = kList, valList = shearList, labelA = crystalName1, labelB = crystalName2, labelC = crystalName3,
    axeslabelx = 'Shear modulus (GPa)', saveFig = "True", type = "shearMod", tmpFolder = tmpFolder)
plot3phase(concList = kList, valList = youngsList, labelA = crystalName1, labelB = crystalName2, labelC = crystalName3,
    axeslabelx = 'Youngs modulus (GPa)', saveFig = "True", type = "youngsMod", tmpFolder = tmpFolder)
