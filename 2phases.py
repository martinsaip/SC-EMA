__author__ = 'Aydin'

from tools.gibbs_plot import plotYoungsModulus
from polycrystal import cubic
from polycrystal import tetragonal
from polycrystal import trigonal
from polycrystal import hexagonal
from polycrystal import orthorombic
import numpy as np

import sys
import os

tmpFolder = str(sys.argv[1])


file = open("../public/"+tmpFolder+"/2phase/crystalStack.txt","r")
crysStack = file.read()
crysList = crysStack.splitlines()
file.close()

file2 = open("../public/"+tmpFolder+"/2phase/concentrationStack.txt","r")
concStack = file2.read()
concList = concStack.splitlines()
file2.close()

crys1 = eval(crysList[0])
crys2 = eval(crysList[1])
conc1 = float(concList[0])
conc2 = float(concList[1])

crys1.setConc(conc1)
crys2.setConc(conc2)

my_crystal = crys1 + crys2
my_crystal.setK0andMue0()
my_crystal.setYoungsMod()
my_crystal.setpoissonratio()

#our_string = sys.argv[1]

#myDict = my_crystal.getElasticDict()

#plotYoungsModulus(my_crystal.crystalname, myDict, saveFig= "True" )

print("The OUTPUT homogenized polycrystalline BULK MODULUS ="+ '<br>','%.2f' % (my_crystal.K0   * 100)        + ' (GPa) ' + '<br>')
print("The OUTPUT homogenized polycrystalline SHEAR MODULUS ="+ '<br>','%.2f' % (my_crystal.mue0 * 100)        + ' (GPa) ' + '<br>')
print("The OUTPUT homogenized polycrystalline YOUNG'S MODULUS = "+ '<br>','%.2f' % (my_crystal.E    * 100)        + ' (GPa) ' + '<br>')
print("The OUTPUT homogenized polycrystalline POISSON RATIO = "+ '<br>','%.4f' % (my_crystal.v)    +  '<br>')
