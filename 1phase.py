#__author__ = 'Aydin & Titrian'
import sys
import os
import tempfile
import time


os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

from tools.gibbs_plot import plotYoungsModulus
from polycrystal import cubic
from polycrystal import tetragonal
from polycrystal import trigonal
from polycrystal import hexagonal
from polycrystal import orthorombic
import numpy as np

our_string = sys.argv[1]
currentUserID = str(sys.argv[2])

#myarray = our_string.split()
#high = int(myarray.pop())
#low  = int(myarray.pop())
#tmpFolder = str(myarray.pop())
#pnts = int(myarray.pop())
#strStruc = ''.join(myarray)

#myIP = tmpFolder.split(".")[0]
#print myIP
#myid = tmpFolder.split(".")[-1]

#1. read the log-file
#try:
#   f = open("../public/logfile.txt",'r+')
#except IOError, e:
#   print e
#
#readList = f.read()
#
#if len(readList) != 0:
#   userIDs = []
#   userTime = []
#   readList = readList.split()
#   for k in len(readList):
   
#   if currentUserID in userIDS:
#	print "yayy found"
#else:
#   f.write('{0:s} {1:s}\n'.format(currentUserID, str(1)))
#
#f.close()
f = open("../public/logfile","a")

my_crystal = eval(our_string)
my_crystal.setK0andMue0()
my_crystal.setYoungsMod()
my_crystal.setpoissonratio()
my_crystal.areamoduli()

myDict = my_crystal.getElasticDict()

#plotYoungsModulus(my_crystal.crystalname, myDict, saveFig= "True" )

crystalname = my_crystal.crystalname

print("Based on input elastic constants, the mechanical stability of your system was checked and the system was found stable (for detail see the &quot;about&quot; page)." + '<br>')
print("The OUTPUT homogenized polycrystalline BULK MODULUS ="+ '<br>','%.2f' % (my_crystal.K0   * 100)        + ' (GPa) ' + '<br>')
print("The OUTPUT homogenized polycrystalline SHEAR MODULUS ="+ '<br>','%.2f' % (my_crystal.mue0 * 100)        + ' (GPa) ' + '<br>')
print("The OUTPUT homogenized polycrystalline YOUNG'S MODULUS = "+ '<br>','%.2f' % (my_crystal.E    * 100)        + ' (GPa) ' + '<br>')
print("The OUTPUT homogenized polycrystalline POISSON RATIO = "+ '<br>','%.4f' % (my_crystal.v)    +  '<br>')
print("If using the REUSS homogenization method, the OUTPUT homogenized polycrystalline BULK MODULUS = "+ '<br>','%.2f' % (my_crystal.reuss_bulk  * 100) + ' (GPa) ' + '<br>')
print("If using the VOIGT homogenization method, the OUTPUT homogenized polycrystalline BULK MODULUS = "+ '<br>','%.2f' % (my_crystal.voigt_bulk  * 100) + ' (GPa) ' + '<br>')
print("If using the REUSS homogenization method, the OUTPUT homogenized polycrystalline SHEAR MODULUS = "+ '<br>','%.2f' % (my_crystal.reuss_shear * 100) + ' (GPa) ' + '<br>')
print("If using the VOIGT homogenization method, the OUTPUT homogenized polycrystalline SHEAR MODULUS = "+ '<br>','%.2f' % (my_crystal.voigt_shear * 100) + ' (GPa) ' + '<br>')

if crystalname == 'cubic': 
    print("Specifically for cubic SINGLE CRYSTALS, the AREA MODUL OF ELASTICITY for the (001) plane is A(001) = "+ '<br>','%.2f' % (my_crystal.A1   * 100)        + ' (GPa) ' + '<br>')
    print("Specifically for cubic SINGLE CRYSTALS, the AREA MODUL OF ELASTICITY for the (110) plane is A(110) = "+ '<br>','%.2f' % (my_crystal.A2   * 100)        + ' (GPa) ' + '<br>')
    print("Specifically for cubic SINGLE CRYSTALS, the AREA MODUL OF ELASTICITY for the (111) plane is A(111) = "+ '<br>','%.2f' % (my_crystal.A3   * 100)        + ' (GPa) ' + '<br>')

f.write("someone used the run button\n")
f.close()
