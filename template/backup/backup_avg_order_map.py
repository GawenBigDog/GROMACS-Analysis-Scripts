#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.patches as ptch
import scipy
import scipy.interpolate

### SOME TOOLS
def closest_node(coord, nodes):
  nodes = np.asarray(nodes)
  dist_2 = np.sum((nodes - coord)**2, axis=1)
  return np.argmin(dist_2)

def neighbor_node(coord, nodes, tol):
  nodes = np.asarray(nodes)
  dist_2 = np.sum((nodes - coord)**2, axis=1)
  neighborhoodIndex = []
  tol_2 = tol**2
  for i, x in enumerate(dist_2):
    if x <= tol_2:
      neighborhoodIndex.append(i)
  if np.size(neighborhoodIndex) == 0:
    neighborhoodIndex.append(np.argmin(dist_2))
  return neighborhoodIndex

def find(lst, a):
  result = []
  for i, x in enumerate(lst):
    if x == a:
      result.append(i)
  return result

def gauss(x,p):
  """
  Return the gauss function N(x), with mean p[0] and std p[1].
  Normalized such that N(x=p[0]) = 1.
  """
  return np.exp((-(x - p[0])**2) / (2 * p[1]**2))
        
### REAL STUFF

if len(sys.argv) < 2:
  sys.stderr.write("Usage: python %s <.dat file>\n" % sys.argv[0])
  exit(0)
  
filename = sys.argv[1]
try:
  sys.argv[2]
except IndexError:
  nodes = 10
else:
  nodes = int(sys.argv[2])

#tol = 2
  
f0 = open(filename,'r')
lines = f0.readlines()[1:]

data = [] # [frame, [x, y], orderP2]

for line in lines:
  data.append([int(line.split()[1]),[float(line.split()[3]),float(line.split()[4])],float(line.split()[6])])

f0.close()

temp_arr = [item[1] for item in data]
pos_arr = np.asarray(temp_arr)
temp_arr = [item[0] for item in data]
fra_arr = np.asarray(temp_arr)

# pos_arr[index][x or y]

xi, yi = np.linspace(np.min(pos_arr[:,0]), np.max(pos_arr[:,0]), nodes), np.linspace(np.min(pos_arr[:,1]), np.max(pos_arr[:,1]), nodes)
xi, yi = np.meshgrid(xi, yi)

meshGrid = []

for i in range(len(xi)):
  for j in range(len(yi)):
    meshGrid.append([xi[i][j], yi[i][j]])

gridOrderP2 = [0.] * (nodes**2)

mes_arr = np.asarray(meshGrid)

# Loop over all the frames
for i in range(np.max(fra_arr)):
  indexList = []

  index = closest_node(pos_arr[i], mes_arr)

  temp_arr = [item[1] for item in data if item[0]==i]
  pos_arr = np.asarray(temp_arr)

  temp_arr = [item[2] for item in data if item[0]==i]
  ord_arr = np.asarray(temp_arr)

  # Find the best matching unit

  for i in range(len(pos_arr)):
    index = closest_node(pos_arr[i], mes_arr)
    indexList.append(index)

#  for i in range(len(pos_arr)):
#    index = neighbor_node(pos_arr[i], mes_arr, tol)
#    indexList.append(index)
  
  counter = [0.] * (nodes**2)
  for j in range(len(indexList)):
    #for k in range(len(indexList[j])):
    index = indexList[j]
    gridOrderP2[index] += ord_arr[j]
    counter[index] += 1

  for j in range(len(gridOrderP2)):
    if (counter[j] != 0):
      gridOrderP2[j] /= counter[j]
    
gridOrderP2 /= np.max(fra_arr)

# Interpolate
rbf = scipy.interpolate.Rbf(xi, yi, gridOrderP2, function='linear')
orderP2_i = rbf(xi, yi)

plt.imshow(orderP2_i, vmin=0, vmax=1, origin='lower',
                                 extent=[np.min(xi), np.max(xi), np.min(yi), np.max(yi)])
plt.scatter(xi, yi, c=gridOrderP2)
plt.colorbar()
plt.clim(0,1)
plt.show()
  
