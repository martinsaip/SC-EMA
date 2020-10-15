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
import numpy as np

inputVal = sys.argv[1]
conc2 = float(inputVal)
conc1 = 1 - conc2

file = open("..\\2phases\\crystalStack.txt","r")
crysStack = file.read()
crysList = crysStack.splitlines()
file.close()

#file2 = open("..\\2phases\\concentrationStack.txt","r")
#concStack = file2.read()
#concList = concStack.splitlines()
#file2.close()
#
our_string = sys.argv[1]

crys1 = eval(crysList[0])
crys2 = eval(crysList[1])

#conc1 = float(concList[0])
#conc2 = float(concList[1])
crystalName1 = crys1.crystalname
crystalName2 = crys2.crystalname

#conc1 = np.linspace(1.0,0.0,11)
#conc2 = np.linspace(0.0,1.0,11)

#for el in conc1:
#    print el

#bulkList = []
#youngsList = []
#shearList = []
#poissonList = []

#kList = []
#conc1 = []
#for k in range(len(conc1)):
crys1.setConc(conc1)
crys2.setConc(conc2)
my_crystal = crys1 + crys2
my_crystal.setK0andMue0()

if conc1 == 1:
    fileNew = open("..\\2phases\\bulkStack.txt", "w")
else:
    fileNew = open("..\\2phases\\bulkStack.txt", "a")
fileNew.write(str(conc2) + " "  + str(my_crystal.K0 * 100) + "\n")
fileNew.close()

#
#steps = 10
#
#for k in range(steps+1):
#    conc = 1.0-k*1.0/steps
#    kList.append([conc, 1-conc])
#    #print conc
#    print( "concentration phase A: %s" % (conc) )
#
#    crys1.setConc(1 - k*1./steps)
#    crys2.setConc(k*1./steps)
##    mynew = crys1 + crys2
#
#    conc1.append(k*1.0/steps)
#    print conc1
#    print bulkList

#for k in range(len(conc1)):
#    crys1.setConc(conc1[k])
#    crys2.setConc(conc2[k])
#
#    my_crystal = crys1 + crys2
#
#    myDict = my_crystal.getElasticDict()
#
#    my_crystal.setK0andMue0()
#    my_crystal.setYoungsMod()
#    my_crystal.setpoissonratio()
#
#    bulkList.append(my_crystal.K0 * 100)
#    shearList.append(my_crystal.mue0 * 100)
#    youngsList.append(my_crystal.E * 100)
#    poissonList.append(my_crystal.v)

#for k in range(len(conc1)):
#    print conc2[k], bulkList[k]
#    print "concentration phase A:" , conc2[k]
#
#plot2phase(concList = conc2, valList = bulkList, labelA = crystalName1, labelB = crystalName2,
#    axesLabelY = "Bulk modulus (GPa)")
#
#plot2phase(concList = conc2, valList = shearList, labelA = crystalName1, labelB = crystalName2,
#    axesLabelY = "Shear modulus (GPa)")
#
#plot2phase(concList = conc2, valList = youngsList, labelA = crystalName1, labelB = crystalName2,
#    axesLabelY = "Youngs modulus (GPa)")
#
##myDict = my_crystal.getElasticDict()
#
##plotYoungsModulus(my_crystal.crystalname, myDict, saveFig= "True" )
#
#print "Bulk Modulus       = ",'%.2f' % (my_crystal.K0   * 100)        + ' (GPa) ' + '<br>'
#print "Shear Modulus      = ",'%.2f' % (my_crystal.mue0 * 100)        + ' (GPa) ' + '<br>'
#print "Young's Modulus    = ",'%.2f' % (my_crystal.E    * 100)        + ' (GPa) ' + '<br>'
#print "Poisson Ratio      = ",'%.4f' % (my_crystal.v)    +  '<br>'
