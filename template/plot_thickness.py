#!/usr/bin/env python

import sys, os
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.patches as ptch
import seaborn

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: %s <datafile> <nres> <resolution> <title>\n")
        sys.exit()
        
datafile = sys.argv[1]
nres = sys.argv[2]
resolution = sys.argv[3]
title = sys.argv[4]

mij = [[0.]*int(nres) for i in range(int(nres))]

f = open(datafile,'r')
lines = f.readlines()
i=0
for f0 in lines:
    j=0
    while j < int(nres):
        mij[i][j] = float(f0.split()[j])
        j+=1
    print int(i)
    i+=1
    if i == int(nres):
        break
f.close()
        
levels = np.linspace(1.0,5.0,resolution)

nres_1 = int(nres)-1

fig = plt.figure()
ax1 = fig.add_subplot(111)
xlim = plt.xlim([0,nres_1])
ylim = plt.ylim([0,nres_1])
CS = plt.contourf(mij,levels=levels,cmap = cm.get_cmap('jet', len(levels)-1))
cbar = fig.colorbar(CS,ticks=[-1.,0,1.])
xdivider1 = plt.plot([0, nres_1], [nres_1, nres_1], color='k', linestyle='-', linewidth=2)
xdivider2 = plt.plot([0, nres_1], [0, 0], color='k', linestyle='-', linewidth=2)
ydivider1 = plt.plot([0, 0], [0, nres_1], color='k', linestyle='-', linewidth=2)
ydivider2 = plt.plot([nres_1, nres_1], [0, nres_1], color='k', linestyle='-', linewidth=2)

plt.title(title, fontsize=20)
plt.xlabel('Membrane Node', fontsize=15)
plt.ylabel('Membrane Node', fontsize=15)
    
#plt.savefig(title+'.png')
#plt.savefig(title+'.eps', format='eps', dpi=1200)

plt.show()
