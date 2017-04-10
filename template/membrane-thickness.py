#!/usr/bin/env python
"""
membrane-thickness.py

Details : compute the local thickness of the membrane
Use-case: create a membrane thickness map
Output  : a 2D membrane map
"""

import sys
import numpy as np
from MDAnalysis import *
from MDAnalysis.analysis.align import *
from MDAnalysis.tests.datafiles import GRO, XTC

# Defines files and arrays here
N = 40
nNodes = N * N

membrane = Universe('PO4.gro', 'PO4.xtc')

centroid_list = []
upper_z_nodes = []
lower_z_nodes = []

xmax = -np.inf
xmin = np.inf
ymax = -np.inf
ymin = np.inf

# Process the trajectories
nFrames = 0
for timestep in membrane.trajectory:
     centroid = membrane.select_atoms("name PO4").center_of_geometry()
     centroid_list.append([centroid[0],centroid[1],centroid[2]])
     for atom in membrane.atoms:
          if ( atom.position[0] - centroid[0] > xmax ):
               xmax = atom.position[0] - centroid[0]
          if ( atom.position[0] - centroid[0] < xmin ):
               xmin = atom.position[0] - centroid[0]
          if ( atom.position[1] - centroid[1] > ymax ):
               ymax = atom.position[1] - centroid[1]
          if ( atom.position[1] - centroid[1] < ymin ):
               ymin = atom.position[1] - centroid[1]
     nFrames += 1

print 'Number of frames: '+str(nFrames)
     
# Create mesh grid here
xi, yi = np.linspace(xmin, xmax, N + int(N/10)), np.linspace(ymin, ymax, N + int(N/10))
xi = xi[int(N/20):-int(N/20)]
yi = yi[int(N/20):-int(N/20)]
xi, yi = np.meshgrid(xi, yi)

xi = np.reshape(xi,nNodes)
yi = np.reshape(yi,nNodes)

# *_node[frame][nodes][X, Y, Z, Count] 
lower_nodes = [[[0,0,0,0] for i in range(nNodes)] for j in range(nFrames)]
upper_nodes = [[[0,0,0,0] for i in range(nNodes)] for j in range(nFrames)]

k = 0
while k < nFrames:
     i = 0
     while i < nNodes:
          lower_nodes[k][i][0] = xi[i]
          lower_nodes[k][i][1] = yi[i]
          upper_nodes[k][i][0] = xi[i]
          upper_nodes[k][i][1] = yi[i]
          i += 1
     k += 1

# Process the upper and lower leaflets
k = 0
for timestep in membrane.trajectory:
     print timestep.time
     for atom in membrane.atoms:
          mindist = np.inf
          minIndx = 0
          
          x = atom.position[0] - centroid_list[k][0]
          y = atom.position[1] - centroid_list[k][1]
          z = atom.position[2] - centroid_list[k][2]

          i = 0
          while (i < nNodes):
               dist = (lower_nodes[k][i][0] - x) ** 2 \
                      + (lower_nodes[k][i][1] - y) ** 2
               if (dist < mindist):
                    mindist = dist
                    minIndx = i
               i += 1
          if ( z < 0 ):
               lower_nodes[k][minIndx][2] += z
               lower_nodes[k][minIndx][3] += 1
          else:
               upper_nodes[k][minIndx][2] += z
               upper_nodes[k][minIndx][3] += 1
     k += 1

# Average local values and calculate thickness
k = 0
while ( k < nFrames ):
     f0 = open(str(k)+'.dat', 'w')
     i = 0
     while ( i < nNodes):
          if (lower_nodes[k][i][3] != 0): 
               lower_nodes[k][i][2] = lower_nodes[k][i][2]/lower_nodes[k][i][3]
          if (upper_nodes[k][i][3] != 0):
               upper_nodes[k][i][2] = upper_nodes[k][i][2]/upper_nodes[k][i][3]
               
          thickness = upper_nodes[k][i][2] - lower_nodes[k][i][2]
          f0.write('{:6.2f} {:10.4f} {:10.4f} {:16.8f}\n'\
                   .format(i, lower_nodes[k][i][0], lower_nodes[k][i][1], thickness))
          i += 1
     f0.close()
     k += 1


          
     
