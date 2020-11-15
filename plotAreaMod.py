from tools.gibbs_plot import plotYoungsModulus, plotAreaModulus
from polycrystal import cubic
from polycrystal import hexagonal
from polycrystal import tetragonal
from polycrystal import trigonal
from polycrystal import orthorombic
#from time import time
import time

import sys
import os
import getpass

mydate = (time.strftime("%d/%m/%Y"))
mytime = (time.strftime("%H:%M:%S"))

logged_user = getpass.getuser()

ourString = sys.argv[1]
myarray = ourString.split()
high = int(myarray.pop())
low  = int(myarray.pop())
tmpFolder = str(myarray.pop())
pnts = int(myarray.pop())
strStruc = ''.join(myarray)
myIP = tmpFolder.split(".")[0]
myid = tmpFolder.split(".")[-1]

if (high == 0) and (low == 0):
    plt = None
elif (low >= high):
    print("warning", "please check your plot range!")
    exit()
else:
    plt = [low, high]
if pnts > 200:
    print("warning", " too many data points! Please use a value below 200 points.")
    exit()

myStruc = eval(strStruc)

myDict = myStruc.getElasticDict()
f = open("../public/logfile", "a")
time1 = time.time()
plotAreaModulus(myStruc.crystalname, myDict, saveFig = "True", plotRange = plt, tmpFolder = tmpFolder, plotPoints  = pnts)
print("Area mod calculated.")
time2 = time.time()
#f.write("plotAreaModulus " + str(time2-time1) + "\n")
f.write("plotAreaModulus  :: runtime: " + str(time2-time1)+".. Date: " + str(mydate)+ "..Time: " +str(mytime)+ " from:" +str(myIP) + "_" + str(myid) + "\n")
f.close()
