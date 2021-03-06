#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import scipy.interpolate
#import seaborn as sns

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <datafile>\n")
        sys.exit()
        
datafile = sys.argv[1]

x = []
y = []
thickness = []

f = open(datafile,'r')
lines = f.readlines()
i=0
print len(lines)
for f0 in lines:
    x.append(f0.split()[1])
    y.append(f0.split()[2])
    thickness.append(f0.split()[3])
f.close()
        
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
               interpolation='none')

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
