#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import scipy.interpolate

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <# of Files> <# Nodes in one direction>\n")
        sys.exit()

nFiles = int(sys.argv[1])
N = int(sys.argv[2])

nNodes = N * N

x = [0] * nNodes
y = [0] * nNodes
thickness = [0] * nNodes

k = 0
while k <= nFiles:

    f = open(str(k)+'.dat','r')
    lines = f.readlines()
    
    i=0
    for line in lines:
        x[i] += float(line.split()[1])
        y[i] += float(line.split()[2])
        thickness[i] += float(line.split()[3])
        i += 1
    f.close()
    k += 1

i = 0
while i < nNodes:
    x[i] = x[i] / nFiles
    y[i] = y[i] / nFiles
    thickness[i] = thickness[i] / nFiles
    i += 1

fig = plt.figure()
ax = fig.add_subplot(111)

nrows, ncols = np.sqrt(np.size(x)),np.sqrt(np.size(x))
x = np.array(x).reshape((nrows, ncols)).astype(np.float)
y = np.array(y).reshape((nrows, ncols)).astype(np.float)
thickness = np.array(thickness).reshape((nrows, ncols)).astype(np.float)

rbf = scipy.interpolate.Rbf(x, y, thickness, function='linear', smooth=1)
thickness_i = rbf(x, y)
p = plt.imshow(thickness_i,
               vmin=np.min(thickness_i), vmax=np.max(thickness_i),
               origin='lower',
               cmap = cm.jet,
               extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
               interpolation='nearest')

#thickness = np.ma.masked_equal(thickness,0)
#p = plt.pcolor(x, y, thickness,
#           cmap=cm.jet,
#           vmin=np.min(thickness), vmax=np.max(thickness))
cb = fig.colorbar(p)
xlim = plt.xlim([np.min(x),np.max(x)])
ylim = plt.ylim([np.min(y),np.max(y)])

plt.xlabel('Membrane Node', fontsize=15)
plt.ylabel('Membrane Node', fontsize=15)
    
#plt.savefig(title+'.png')
#plt.savefig(title+'.eps', format='eps', dpi=1200)

#sns.set_palette("husl")

plt.show()
