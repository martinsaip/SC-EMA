__author__ = 'Aydin & Titrian'

from tools.gibbs_plot import plotYoungsModulus
from tools.gibbs_plot import plot2phase

from polycrystal import cubic
from polycrystal import tetragonal
from polycrystal import trigonal
from polycrystal import hexagonal
from polycrystal import orthorombic
import numpy as np

import sys
import os

tmpFolder = str(sys.argv[1]);

file = open("../public/"+tmpFolder+"/2phase/crystalStack.txt","r")
crysStack = file.read()
crysList = crysStack.splitlines()
file.close()

crys1 = eval(crysList[0])
crys2 = eval(crysList[1])
#conc1 = float(concList[0])
#conc2 = float(concList[1])
crystalName1 = crys1.crystalname
crystalName2 = crys2.crystalname

conc1 = np.linspace(1,0,11)
conc2 = np.linspace(0,1,11)


bulkList = []
youngsList = []
shearList = []
poissonList = []

for k in range(len(conc1)):
    crys1.setConc(conc1[k])
    crys2.setConc(conc2[k])
#    print("concentration phase II: %s" % (conc2[k]) )+ '<br>')
    my_crystal = crys1 + crys2

    myDict = my_crystal.getElasticDict()

    my_crystal.setK0andMue0()
    my_crystal.setYoungsMod()
    my_crystal.setpoissonratio()

    bulkList.append(my_crystal.K0 * 100)
    shearList.append(my_crystal.mue0 * 100)
    youngsList.append(my_crystal.E * 100)
    poissonList.append(my_crystal.v)

#    print(conc2[k], bulkList [k])

#write to separate table
f = open("../public/"+tmpFolder+"/2phase/dataTable.txt","w")
f.write("#conc. Phase II | Bulk Mod (GPa) | Shear Mod (GPa) | Youngs Mod (GPa) | poiss. ratio\n")
for k in range(len(conc2)):
#    line = "%f" % conc1[k]
    line = '{0:f}\t  {1:f}\t   {2:f}\t     {3:f}\t\t{4:f}\n'.format(conc2[k], float(bulkList[k]), float(shearList[k]),
        float(youngsList[k]), float(poissonList[k]))
    f.write(line)
f.close()

plot2phase(concList = conc2, valList = bulkList, labelA = crystalName1, labelB = crystalName2, axesLabelY = "Bulk modulus (GPa)",saveFig= "True", type = "bulkMod", tmp = tmpFolder)
plot2phase(concList = conc2, valList = shearList, labelA = crystalName1, labelB = crystalName2, axesLabelY = "Shear modulus (GPa)",saveFig= "True", type = "shearMod", tmp = tmpFolder)
plot2phase(concList = conc2, valList = youngsList, labelA = crystalName1, labelB = crystalName2, axesLabelY = "Youngs modulus (GPa)",saveFig= "True", type = "youngsMod", tmp = tmpFolder)

#print("The OUTPUT homogenized polycrystalline BULK MODULUS ="    + '<br>','%.2f' % (my_crystal.K0   * 100)  + ' (GPa) ' + '<br>')
#print("The OUTPUT homogenized polycrystalline SHEAR MODULUS ="   + '<br>','%.2f' % (my_crystal.mue0 * 100)  + ' (GPa) ' + '<br>')
#print("The OUTPUT homogenized polycrystalline YOUNG'S MODULUS = "+ '<br>','%.2f' % (my_crystal.E    * 100)  + ' (GPa) ' + '<br>')
#print("The OUTPUT homogenized polycrystalline POISSON RATIO = "  + '<br>','%.4f' % (my_crystal.v)  +  '<br>')
