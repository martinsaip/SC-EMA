import os, sys, tempfile

os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
try:
   import numpy as np
   import matplotlib
   matplotlib.use("agg")
   import matplotlib.pyplot as plt
   from mpl_toolkits.mplot3d import Axes3D
   from matplotlib import cm, colors
#except RuntimeError, e:
except:
   print("RuntimeError", e)
   sys.exit(1)
var = sys.argv[1]

print(int(var)*10)
