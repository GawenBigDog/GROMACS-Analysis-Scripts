"""
Probability Density Function
"""

import numpy as np
import sys, os

from MDAnalysis import *
from MDAnalysis.lib.distances import *
import MDAnalysis.units
from MDAnalysis.tests.datafiles import GRO, XTC
import seaborn
import imageio

try:
    import matplotlib
    
    matplotlib.use('Qt4Agg')  # Qt4Agg for plots and agg for no plots
    import pylab
    
    have_matplotlib = True
except ImportError:
    have_matplotlib = False

import matplotlib.pyplot as plt
    
if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write("Usage: %s <gro file> <xtc file> <title>\n")
        sys.exit()

dmin, dmax = 0.0, 100.0
nbins = 101
binWidth = (dmax-dmin)/nbins
        
groFile = sys.argv[1]
xtcFile = sys.argv[2]
title = sys.argv[3]

try:
    sys.argv[4]
except IndexError:
    skipFlag = False
else:
    skipFlag = True
    skip = sys.argv[4]

universe = Universe(groFile, xtcFile)

pdf = np.zeros(nbins)

distance=np.linspace(dmin,dmax,nbins)
volume=np.linspace(dmin,dmax,nbins)

# loop over all frames
frames=0
boxvolume=0
timeStepCounter = 0

# Volume calculation is expensive.  Finish it here
volume[0] = 0
i=1
while i < nbins:
    volume[i] = (4/3)*np.pi*(np.power(distance[i],3) - (4/3)*np.power(distance[i-1],3))
    i += 1

print "##################################"

for timestep in universe.trajectory:
    if (skipFlag) :
        if (timeStepCounter%int(skip) != 0) :
            timeStepCounter += 1
            continue

    pdf = np.zeros(nbins)
    i=1
    while i < nbins:
        selection = "name ROH and sphlayer "+str(distance[i-1])+" "+str(distance[i])+" (name BB and bynum 1938:3672)"
        #selection = "name ROH and around "+str(distance[i-1])+" (name BB and bynum 1:3672)"
        group = universe.select_atoms(selection)
        print "Frame     : "+str(timeStepCounter)
        print "Selection : "+str(selection)
        pdf[i] = group.n_residues
        i += 1
        print "##################################"
    frames = 1
    boxvolume = timestep.volume
    normGroup = universe.select_atoms("name ROH").n_residues
    normVolume = boxvolume/frames
    normDensity = normGroup/normVolume

    # Calculate the normalization factor
    normFactor = np.zeros(nbins)
    normFactor[0] = 0
    i=1
    while i < nbins:
        normFactor[i] = normDensity * volume[i]
        i += 1

    pdf /= frames

    i=1
    while i < nbins:
        pdf[i] = pdf[i] / normFactor[i]
        i += 1
    
    if np.isnan(pdf[0]):
        pdf[0] = 0.0
    pdf = pdf / np.sum(pdf)
    #pdf = (count - np.min(count))/(np.max(count)-np.min(count))

    print "##################################"
    print "Distance Array :"
    print distance
    print "##################################"
    print "PDF Array :"
    print pdf
    print "##################################"
    print np.sum(pdf)

    np.savetxt(str(title)+"-pdf"+str(timeStepCounter)+".dat", (distance, pdf))

    plt.bar(distance, pdf, color='k')
    plt.title(title)
    plt.xlabel(r"distance $r$ in nm")
    plt.ylabel(r"probability")
    plt.savefig(str(title)+"-pdf"+str(timeStepCounter)+".png")

    plt.clf()

    timeStepCounter += 1

images = []
i=0
while (i <= timeStepCounter):
    filename=str(title)+"-pdf"+str(timeStepCounter)+".png"
    images.append(imageio.imread(filename))
    i += 1
imageio.mimsave(str(title)+"-pdf"+str(timeStepCounter)+".gif", images)
