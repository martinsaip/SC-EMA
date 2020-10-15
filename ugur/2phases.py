__author__ = 'Aydin & Titrian'

from tools.gibbs_plot import plotYoungsModulus
from polycrystal import cubic
from polycrystal import tetragonal
from polycrystal import trigonal
from polycrystal import hexagonal
from polycrystal import orthorombic
import numpy as np

import sys
import os

#this works for windows
file = open("..\\2phases\\crystalStack.txt","r")
crysStack = file.read()
crysList = crysStack.splitlines()
file.close()
file2 = open("..\\2phases\\concentrationStack.txt","r")
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
#my_crystal.areamoduli()

#our_string = sys.argv[1]

#myDict = my_crystal.getElasticDict()

#plotYoungsModulus(my_crystal.crystalname, myDict, saveFig= "True" )

print("Bulk Modulus       = ",'%.2f' % (my_crystal.K0   * 100)        + ' (GPa) ' + '<br>')
print("Shear Modulus      = ",'%.2f' % (my_crystal.mue0 * 100)        + ' (GPa) ' + '<br>')
print("Young's Modulus    = ",'%.2f' % (my_crystal.E    * 100)        + ' (GPa) ' + '<br>')
print("Poisson Ratio      = ",'%.4f' % (my_crystal.v)    +  '<br>')
