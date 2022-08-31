__author__ = 'aydin'
import sys
import os
import tempfile

os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

import getpass

logged_user = getpass.getuser()

import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from time import time
#from mpl_toolkits.mplot3d import Axes3D

a = [1.,0.]
b = [np.cos(60./180 * np.pi), np.sin(60./180. * np.pi)]

global S

def pos( v ):
    """
    The function returns the x and y coordinate
    as an array element. E.g., [x,y]
    in the triangular plot.
    v is the given concentration
    """
    global a, b
    return [v[1] * a[0] + v[2]* b[0], v[1] * a[1] + v[2] * b[1]]

def retX( v ):
    """
    v is the concentration value, e.g., v = [1.0, 0.0, 0.0].
    The function returns the x coordinate in the triangular
    for the given concentration (argument) v.
    """
    global a, b
    return (v[1] * a[0] + v[2] * b[0])

def retY( v ):
    """
    The function returns the y coord in the triangular
    for the given concentration (arg) v.
    v is the concentration value, e.g., v = [1.0, 0.0, 0.0].
    """
    global a, b
    return  (v[1] * a[1] + v[2] *b[1])

# TODO: this is an arbitrary function for test reasons
# and must be deleted later
def func( v ):
    """
    this is an arbitrary function for test reasons
    x^2 + y^2 + sin(z)
    x,y,z = concentrations of first, second, and third element
    returns the x,y tuple of the triangle plot
    """
    global a, b
    return v[0]**3 + v[1]**4 + np.sin(1-v[0]-v[1])**2

def plot3phase(concList = None, valList = None, labelA = None, labelB = None, labelC = None, axeslabelx = None,
               saveFig = "False", type = None, tmpFolder = None):
    """
    This function plots a triangular (Gibbs diagram) for 3 phases.
    It shows the Youngs modulus with respect to the concentration.
    """
    typeList = ["bulkMod", "shearMod", "youngsMod"]

    if type not in typeList:
        print("type is missing: Must be bulkMod, shearMod or youngsMod")
        sys.exit()

    if len(concList) != len(valList):
        print("error in length")

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    x = [ retX(concList[k]) for k in range(len(concList))]
    y = [ retY(concList[k]) for k in range(len(concList))]
    z = [ float(valList[k]) for k in range(len(valList))]
#    print len(x)
#    print len(y)
#    print len(z)

    plt.plot(x,y, 'kx')
    plt.figtext(0.05, 0.04, labelA)
    plt.figtext(0.75, 0.04, labelB)
    plt.figtext(0.42, 0.92, labelC)

    plt.tricontourf(x,y,z, 45, cmap = plt.cm.jet)
    plt.colorbar()

    plt.xlabel(axeslabelx)
#    plt.show()
#    plt.draw()

    if saveFig == "True":
        if type == "bulkMod":
            if tmpFolder == None:
                plt.savefig("../public/3phase/Bulk_modulus.png")
            else:
                plt.savefig("../public/"+tmpFolder+"/3phase/Bulk_modulus.png")
        elif type == "shearMod":
            if tmpFolder == None:
                plt.savefig("../public/3phase/Shear_modulus.png")
            else:
                plt.savefig("../public/"+tmpFolder+"/3phase/Shear_modulus.png")
        elif type == 'youngsMod':
            if tmpFolder == None:
                plt.savefig("../public/3phase/Youngs_modulus.png")
            else:
                plt.savefig("../public/"+tmpFolder+"/3phase/Youngs_modulus.png")

    else:
        plt.show()
        plt.draw()


def vecProd(pos = None): # TODO: hard coded here, change S value
    global S
    a = pos[0]
    b = pos[1]
    c = pos[2]
    myvec = np.matrix([a**2, b**2, c**2, b*c, c*a, a*b])
    tmyvec = np.transpose(myvec)
    return float(1./(np.transpose(S * tmyvec ) * tmyvec))

def modifyVecs( x = None, y = None, z = None):
    newX = np.copy(x)
    newY = np.copy(y)
    newZ = np.copy(z)
    for i in range(len(x)):
        for j in range(len(x)):
            val = [x[i][j], y[i][j], z[i][j]]
            newVal = np.dot(vecProd(val), val )
            newX[i][j] = newVal[0]
            newY[i][j] = newVal[1]
            newZ[i][j] = newVal[2]
    return newX, newY, newZ

def sphere(dim = 0):
#    zvalues = np.linspace(-1, 1, dim+1)
    anglesAlt = np.linspace(np.pi, 0, dim+1)
    anglesAz = np.linspace(2*np.pi, 0, dim+1)
    z = np.matrix([np.cos(anglesAlt) for _ in range(dim+1)])

    xmat = np.matrix([[np.sin(valAlt)* np.cos(valAz) for valAz in anglesAz] for valAlt in anglesAlt])
    ymat = np.matrix([[np.sin(valAlt)*np.sin(valAz) for valAz in anglesAz] for valAlt in anglesAlt])
    zmat = np.transpose(z)

    return xmat,ymat,zmat

def getColors(x = None, y = None, z = None):
    colList = np.copy(x)
    valList = []
    for i in range(len(x)):
        for j in range(len(x)):
            val = np.array([x[i][j], y[i][j], z[i][j]])
            valList.append(np.linalg.norm(val))
        #            colList[i][j] = np.linalg.norm(val)
    return valList
def prod(pos = None):
    global S
    global indexMat
    global kroneckerMat
    indexMat = np.matrix([[0,5,4],[5,1,3],[4,3,2]])
    kroneckerMat = np.matrix([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    a = pos[0]
    b = pos[1]
    c = pos[2]

#    CC = [a,b,c]
#    sum = 0
#    for i in range(3):
#        for j in range(3):
#            for k in range(3):
#                for l in range(3):
#                    pref = 1.0
#                    if indexMat[i,j] > 2: pref *= 0.5
#                    if indexMat[k,l] > 2: pref *= 0.5
#                    sum = sum + pref * (S[indexMat[i,j],indexMat[k,l]]*(kroneckerMat[i,j] - CC[i]*CC[j])*(kroneckerMat[k,l] - CC[k]*CC[l]))
#    return 1./sum

    #Testing for cubic case
    myvec = np.matrix([1-a**2, 1-b**2, 1-c**2, -b*c, -c*a,  -a*b])
    tmyvec = np.transpose(myvec)
    return float(1./(np.transpose(S * tmyvec ) * tmyvec))

def modify( x = None, y = None, z = None):
    newX = np.copy(x)
    newY = np.copy(y)
    newZ = np.copy(z)
    for i in range(len(x)):
        for j in range(len(x)):
            val = [x[i][j], y[i][j], z[i][j]]
            newVal = np.dot(prod(val), val)
            newX[i][j] = newVal[0]
            newY[i][j] = newVal[1]
            newZ[i][j] = newVal[2]
    return newX, newY, newZ

#def plotAreaModulus(crystalType = None, vals = None, file = None):
def plotAreaModulus(crystalType = None, vals = None, saveFig = "False", plotRange = None, tmpFolder = None, plotPoints = None):
    timeA = time()
    time1 = time()
    global S
    if crystalType == 'cubic':
        c11 = vals['C11']; c12 = vals['C12']; c44 = vals['C44']
        mat = np.matrix([[c11,c12,c12,0,0,0],[c12,c11,c12,0,0,0],[c12,c12,c11,0,0,0],
                         [0,0,0,c44,0,0],[0,0,0,0,c44,0],[0,0,0,0,0,c44]])
    elif crystalType == 'hexagonal':
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']
        c33 = vals['C33']; c55 = vals['C55']
        c66 = 0.5 * (c11 - c12)
        mat = np.matrix([[c11,c12,c13, 0, 0, 0],[c12,c11,c13, 0, 0, 0],[c13,c13,c33, 0, 0, 0],
                         [0, 0, 0,c55, 0,  0],[0, 0, 0, 0,c55,  0],[0, 0, 0, 0, 0,c66]])
    elif crystalType == "tetragonal":
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']
        c33 = vals['C33']; c44 = vals['C44']; c66 = vals['C66']
        mat = np.matrix([[c11,c12,c13, 0, 0, 0],[c12,c11,c13, 0, 0, 0],[c13,c13,c33, 0, 0, 0],
                         [0, 0, 0,c44, 0, 0],[0, 0, 0, 0,c44, 0],[0, 0, 0, 0, 0,c66]])
    elif crystalType == "trigonal":
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']
        c14 = vals['C14']; c33 = vals['C33']; c44 = vals['C44']
        c66 = 0.5* (c11 - c12)
        mat = np.matrix([[c11, c12, c13, c14, 0, 0],[c12, c11, c13,-c14, 0, 0],[c13, c13, c33, 0, 0, 0],
                         [c14,-c14, 0, c44, 0, 0], [ 0, 0, 0, 0, c44, c14/2.],[ 0, 0, 0, 0, c14/2., c66]])
    elif crystalType == "orthorombic":
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']; c22 = vals['C22']
        c23 = vals['C23']; c33 = vals['C33']; c44 = vals['C44']
        c55 = vals['C55']; c66 = vals['C66'];
        mat = np.matrix([[c11, c12, c13, 0, 0, 0],[c12, c22, c23, 0, 0, 0],[c13, c23, c33, 0, 0, 0],
                         [0, 0, 0, c44, 0, 0],[0, 0, 0, 0, c55, 0],[0, 0, 0, 0, 0, c66]])
    else:
        print("nothing to do")
    S = np.linalg.inv(mat)
    time2 = time()
    #print "I inverted the matrix", time2 - time1

    if plotPoints == None:
        plotPoints = 50
    [x,y,z] = sphere(plotPoints)
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
#    time1 = time()
    time1 = time()
    [x,y,z] = modify(x,y,z)
 #   time2 = time()

#    time2 = time()
#    print time2 - time1

#    if plotRange == None:
##        lowerRange = np.min([np.min(x),np.min(y),np.min(z)])*1.1
 #       upperRange = np.max([np.max(x),np.max(y),np.max(z)])*1.1
 #   else:
 ##       lowerRange = np.float64(-plotRange[1])#np.float64(-plotRange[1])
 #       upperRange = np.float64(plotRange[1]) 
    #lowerRange = np.min([np.min(x),np.min(y),np.min(z)])*1.1
    #upperRange = np.max([np.max(x),np.max(y),np.max(z)])*1.1
 #   cls = getColors(x,y,z)

 #   if plotRange == None:
 #       MaxVal = max(cls)
 #       MinVal = min(cls)
 #   else:
 #       MinVal = np.float64(plotRange[0])
 #       MaxVal = np.float64(plotRange[1])

#    cls_old = getColors(x,y,z)
#    MaxVal = max(cls)
#    MinVal = min(cls)
    #    print min(cls), max(cls)

    #low = lowerRange
    #high = upperRange
#
#    MaxValColor = np.float64(high)
#    MinValColor = np.float64(low)
#    MinValCLS = np.float64(low)
#    MaxValCLS = np.float64(high)

#    cls = (cls - min(cls))/(max(cls)-min(cls))#/np.max(cls)

    #cls = (cls - MinValCLS)/(MaxValCLS - MinValCLS)

    #cls = np.reshape(cls,(len(x),len(x)))

    #fig = plt.figure()#facecolor=None)
    #ax = Axes3D(fig)
    #ax = fig.gca(projection='3d')

    #ax.set_axis_bgcolor("#bdb76b")
    #time1 = time()

    #surf = ax.plot_surface(x,y,z, rstride=1, cstride=1, shade = True, antialiased = True,facecolors = cm.gist_rainbow(cls))

    #time2 = time()

    #axesX = ax.axes.set_xlabel("[100]")
    #axesY = ax.axes.set_ylabel("[010]")
    #axesZ = ax.axes.set_zlabel("[001]")

    #    m = cm.ScalarMappable(cmap=cm.RdYlBu)
    #m = cm.ScalarMappable(cmap = cm.gist_rainbow)
    ##    cls = cls * (MaxValColor - MinValColor) + MinValColor
    #cls = cls * (MaxVal - MinVal) + MinVal

    #mycls = np.ndarray(shape  = (len(cls),len(cls)))

#        for k in range(25):
#            mycls[k,:] = k * (0.1 - 0.0) / (25 - 1)
#            mycls[k,:] = k * 1.0 /(len(cls)-1)
#        for k in range(25,51):
#            mycls[k,:] = k * (1.0 - 0.0) / (51 - 1)
#    for k in range(len(x)):
#        for l in range(len(x)):
#            mycls[k][l] = (k * 1.0 + l * 0.1 /(len(x) - 1))  /(len(x) - 1)
    #if plotRange == None:
    #    cls = cls * (MaxVal - MinVal) + MinVal
    #else:
    ##    mycls = np.ndarray(shape = (len(cls),len(cls)))
    #    for k in range(len(cls)):
    #       mycls[k,:] = k * 1.0 / (len(cls)-1)
    #    cls = mycls * (MaxVal - MinVal) + MinVal

#    mycls = mycls * (MaxValColor - MinValColor) + MinValColor

    #m.set_array(cls)
#    m.set_array(mycls)
    #fig.colorbar(m, shrink = 0.8)

    #ax.set_autoscalez_on(False)
    #ax.set_autoscalex_on(False)
    #ax.set_autoscaley_on(False)
    #    print MinVal, MaxVal
    #print "I am still here"
    #ax.autoscale_view( scalex = False, scaley = False, scalez = False) # for matplotlib version 1.0.xx

#    ax.set_xlim([lowerRange,upperRange])
#    ax.set_ylim([lowerRange,upperRange])
#    ax.set_zlim([lowerRange,upperRange])
    #ax.set_xlim([-high, high])
    #ax.set_ylim([-high, high])
    #ax.set_zlim([-high, high])
    #lowerRange = -high
    #upperRange = high
    #ax.set_xlim3d([lowerRange,upperRange])
    #ax.set_ylim3d([lowerRange,upperRange])
    #ax.set_zlim3d([lowerRange,upperRange])

    #if file == None:
    #    plt.show()
    #else:
    #    plt.savefig(file)
    #time2 = time()
    #print "I am saving", time2 - time1

    #if saveFig == "True":
#	try:
#	   os.remove("../public/" + str(tmpFolder) + "/" + str(crystalType) + "_area.png")
#	except OSError:
#	   pass
#	plt.savefig("../public/" + str(tmpFolder) + "/" + str(crystalType) + "_area.png")
#    else:
#        plt.show()

    # scale arragement (lifang calculation -200 & 200 )
    if plotRange == None:
        lowerRange = np.min([np.min(x),np.min(y),np.min(z)])*1.1
        upperRange = np.max([np.max(x),np.max(y),np.max(z)])*1.1
    else:
        lowerRange = np.float64(-plotRange[1])
        upperRange = np.float64(plotRange[1])

#    cls_old = getColors(x,y,z)
    cls = getColors(x,y,z)

#    a = [np.float64(6000.),np.float64(12000.)]
#    a = [np.float64(50),np.float64(230)]

#    MaxVal = 230.0 #max(cls)
#    MinVal = 50.0 #min(cls)
    if plotRange == None:
        MaxVal = max(cls)
        MinVal = min(cls)
    else:
        MinVal = np.float64(plotRange[0])
        MaxVal = np.float64(plotRange[1])

    #cls = (cls - min(cls))/(max(cls)-min(cls))#/np.max(cls)
    cls = (cls - MinVal)/(MaxVal - MinVal)
#    cls = (cls - a[0])/((a[1]-a[0]))#/np.max(cls)
#    cls = (cls - a[0]) / a[1]
    cls = np.reshape(cls,(len(x),len(x)))

    fig = plt.figure()#facecolor=None)
    ax = Axes3D(fig)
    ax = fig.gca(projection='3d')

    #ax.set_axis_bgcolor("#bdb76b")
    surf = ax.plot_surface(x,y,z, rstride=1, cstride=1, shade = True, antialiased = True,facecolors = cm.gist_rainbow(cls))
    axesX = ax.axes.set_xlabel("[100]")
    axesY = ax.axes.set_ylabel("[010]")
    axesZ = ax.axes.set_zlabel("[001]")

    m = cm.ScalarMappable(cmap=cm.gist_rainbow)
    
    if plotRange == None:
        cls = cls * (MaxVal - MinVal) + MinVal
    else:
        mycls = np.ndarray(shape = (len(cls),len(cls)))
        for k in range(len(cls)):
           mycls[k,:] = k * 1.0 / (len(cls)-1)
        cls = mycls * (MaxVal - MinVal) + MinVal
    m.set_array(cls) 
    fig.colorbar(m, shrink = 0.8)

    #ax.set_autoscalez_on(False)
    #ax.set_autoscalex_on(False)
    #ax.set_autoscaley_on(False)

    ax.autoscale_view( scalex = False, scaley = False, scalez = False) # for matplotlib version 1.0.xx

#    print MinVal, MaxVal
    ax.set_xlim3d([lowerRange,upperRange])
    ax.set_ylim3d([lowerRange,upperRange])
    ax.set_zlim3d([lowerRange,upperRange])

    if saveFig == "True":
        try:
           os.remove("../public/" + str(tmpFolder) + "/" + str(crystalType) + "_area.png")
        except OSError:
           pass
        plt.savefig("../public/" + str(tmpFolder) + "/" + str(crystalType) + "_area.png")
    else:
        plt.show()

def plotYoungsModulus(crystalType = None, vals = None, saveFig = "False", plotRange = None, tmpFolder = None, plotPoints = None):
    global S
    if crystalType == 'cubic':
        c11 = vals['C11']; c12 = vals['C12']; c44 = vals['C44']
        mat = np.matrix([[c11,c12,c12,0,0,0],[c12,c11,c12,0,0,0],[c12,c12,c11,0,0,0],
            [0,0,0,c44,0,0],[0,0,0,0,c44,0],[0,0,0,0,0,c44]])
    elif crystalType == 'hexagonal':
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']
        c33 = vals['C33']; c55 = vals['C55']
        c66 = 0.5 * (c11 - c12)
        mat = np.matrix([[c11,c12,c13, 0, 0, 0],[c12,c11,c13, 0, 0, 0],[c13,c13,c33, 0, 0, 0],
                         [0, 0, 0,c55, 0,  0],[0, 0, 0, 0,c55,  0],[0, 0, 0, 0, 0,c66]])
    elif crystalType == "tetragonal":
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']
        c33 = vals['C33']; c44 = vals['C44']; c66 = vals['C66']
        mat = np.matrix([[c11,c12,c13, 0, 0, 0],[c12,c11,c13, 0, 0, 0],[c13,c13,c33, 0, 0, 0],
                         [0, 0, 0,c44, 0, 0],[0, 0, 0, 0,c44, 0],[0, 0, 0, 0, 0,c66]])
    elif crystalType == "trigonal":
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']
        c14 = vals['C14']; c33 = vals['C33']; c44 = vals['C44']
        c66 = 0.5* (c11 - c12)
        mat = np.matrix([[c11, c12, c13, c14, 0, 0],[c12, c11, c13,-c14, 0, 0],[c13, c13, c33, 0, 0, 0],
                         [c14,-c14, 0, c44, 0, 0], [ 0, 0, 0, 0, c44, c14/2.],[ 0, 0, 0, 0, c14/2., c66]])
    elif crystalType == "orthorombic":
        c11 = vals['C11']; c12 = vals['C12']; c13 = vals['C13']; c22 = vals['C22']
        c23 = vals['C23']; c33 = vals['C33']; c44 = vals['C44']
        c55 = vals['C55']; c66 = vals['C66'];
        mat = np.matrix([[c11, c12, c13, 0, 0, 0],[c12, c22, c23, 0, 0, 0],[c13, c23, c33, 0, 0, 0],
                         [0, 0, 0, c44, 0, 0],[0, 0, 0, 0, c55, 0],[0, 0, 0, 0, 0, c66]])
    else:
        print("nothing to do")
    S = np.linalg.inv(mat)

    if plotPoints == None:
        plotPoints = 50

    [x,y,z] = sphere(plotPoints)
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    [x,y,z] = modifyVecs(x,y,z)

    # scale arragement (lifang calculation -200 & 200 )
    if plotRange == None:
        lowerRange = np.min([np.min(x),np.min(y),np.min(z)])*1.1
        upperRange = np.max([np.max(x),np.max(y),np.max(z)])*1.1
    else:
        lowerRange = np.float64(-plotRange[1])
        upperRange = np.float64(plotRange[1])

#    cls_old = getColors(x,y,z)
    cls = getColors(x,y,z)

#    a = [np.float64(6000.),np.float64(12000.)]
#    a = [np.float64(50),np.float64(230)]

#    MaxVal = 230.0 #max(cls)
#    MinVal = 50.0 #min(cls)
    if plotRange == None:
        MaxVal = max(cls)
        MinVal = min(cls)
    else:
        MinVal = np.float64(plotRange[0])
        MaxVal = np.float64(plotRange[1])

    #cls = (cls - min(cls))/(max(cls)-min(cls))#/np.max(cls)
    cls = (cls - MinVal)/(MaxVal - MinVal)
#    cls = (cls - a[0])/((a[1]-a[0]))#/np.max(cls)
#    cls = (cls - a[0]) / a[1]
    cls = np.reshape(cls,(len(x),len(x)))

    fig = plt.figure()#facecolor=None)
    ax = Axes3D(fig)
    ax = fig.gca(projection='3d')

    #ax.set_axis_bgcolor("#bdb76b")
    surf = ax.plot_surface(x,y,z, rstride=1, cstride=1, shade = True, antialiased = True,facecolors = cm.gist_rainbow(cls))
    axesX = ax.axes.set_xlabel("[100]")
    axesY = ax.axes.set_ylabel("[010]")
    axesZ = ax.axes.set_zlabel("[001]")

    m = cm.ScalarMappable(cmap=cm.gist_rainbow)
    
    if plotRange == None:
    	cls = cls * (MaxVal - MinVal) + MinVal
    else:
        mycls = np.ndarray(shape = (len(cls),len(cls)))
        for k in range(len(cls)):
           mycls[k,:] = k * 1.0 / (len(cls)-1)
        cls = mycls * (MaxVal - MinVal) + MinVal
    m.set_array(cls) 
    fig.colorbar(m, shrink = 0.8)

    #ax.set_autoscalez_on(False)
    #ax.set_autoscalex_on(False)
    #ax.set_autoscaley_on(False)

    ax.autoscale_view( scalex = False, scaley = False, scalez = False) # for matplotlib version 1.0.xx

#    print MinVal, MaxVal
    ax.set_xlim3d([lowerRange,upperRange])
    ax.set_ylim3d([lowerRange,upperRange])
    ax.set_zlim3d([lowerRange,upperRange])


    if saveFig == "True":
        try:
           os.remove("../public/" + str(tmpFolder) + "/" + str(crystalType) + ".png")
        except OSError:
           pass
        plt.savefig("../public/" + str(tmpFolder) + "/" + str(crystalType) + ".png")
    else:
        plt.show()

def plotHashinShtriktman():
    pass

def plotVoigtAndReuss(concList = None, valList = None, crystalType = None, mynew = None):
    if len(concList) != len(valList):
        print("error in length")
        sys.exit()

    if concList == 0:
#        crystalType = crys_1
#        self.pt1Dbulk = (mynew.voigt_bulk)
#        self.pt1Ubulk = (mynew.reuss_bulk)
#        self.pt1Dshear = mynew.voigt_shear
#        self.pt1Ushear = mynew.reuss_shear
        pass
    if concList == 1:
#        crystalType = crys_2
#        self.pt2Dbulk = (mynew.voigt_bulk)
#        self.pt2Ubulk = (mynew.reuss_bulk)
#        self.pt2Dshear = mynew.voigt_shear
#        self.pt2Ushear = mynew.reuss_shear
        pass

def plot2phase(concList = None, valList = None, labelA = None, labelB = None,
               axesLabelY = None, saveFig = "False", type = None, tmp = None):
    typeList = ["bulkMod", "shearMod", "youngsMod"]

    if type not in typeList:
        print("type is missing: Must be bulkMod, shearMod or youngsMod")
        sys.exit()

    if len(concList) != len(valList):
        print("error in length")
        sys.exit()

    axeslabelx = "x(" + str(labelB) +")"
    #axesLabely = "Youngs modulus (GPa)"

    plt.plot(concList, valList)
    plt.figtext(0.10,0.03, labelA)
    plt.figtext(0.89, 0.03, labelB)
    plt.xlabel(axeslabelx)
    plt.ylabel(axesLabelY)
#    plt.plot(0,pt1Dbulk,'ro')
#    plt.plot(0,pt1Ubulk,'bo')
#    plt.plot(1,pt2Dbulk,'ro')
#    plt.plot(1,pt2Ubulk,'bo')
    
    saveName = ""
    if type == "bulkMod": saveName = "Bulk_modulus.png"
    elif type == "shearMod": saveName = "Shear_modulus.png"
    elif type == "youngsMod": saveName = "Youngs_modulus.png"

    if saveFig == "True":
        try:
            if tmp == None:
                plt.savefig("../public/2phase/"+saveName)
            else:
                plt.savefig("../public/"+tmp+"/2phase/"+saveName)
        except IOError as e:
            print("Error: ", e)
    else:
        plt.show()

    plt.clf()


#if __name__ == '__main__':
#    vec = [] # concentration triplet
#    for k in range(11):
#        for l in range(11):
#            for m in range(11):
#                if k + l + m == 10:
#                    vec.append([k/10.,l/10.,m/10.])

    #TODO: this must be changed! z values must be calculated

#    z = [ func(vec[k]) for k in range(len(vec))]

#    k = np.linspace(0,10,30)
#    l = np.ones(30)
#    plot2phase(k,l)
#    plot3phase(vec, z, type = "bulkMod")
