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
orderP2 = []

f = open(datafile,'r')
lines = f.readlines()
i=0
print len(lines)
for line in lines:
    x.append(line.split()[1])
    y.append(line.split()[2])
    orderP2.append(line.split()[3])
f.close()
        
fig = plt.figure()
ax = fig.add_subplot(111)

print np.size(x)

nrows, ncols = np.sqrt(np.size(x)),np.sqrt(np.size(x))
x = np.array(x).reshape((nrows, ncols)).astype(np.float)
y = np.array(y).reshape((nrows, ncols)).astype(np.float)
orderP2 = np.array(orderP2).reshape((nrows, ncols)).astype(np.float)

rbf = scipy.interpolate.Rbf(x, y, orderP2, function='linear', smooth=1)
orderP2_i = rbf(x, y)
p = plt.imshow(orderP2_i,
               vmin=np.min(orderP2_i), vmax=np.max(orderP2_i),
               origin='lower',
               cmap = cm.jet,
               extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
               interpolation='none')

#orderP2 = np.ma.masked_equal(orderP2,0)
#p = plt.pcolor(x, y, orderP2,
#           cmap=cm.jet,
#           vmin=np.min(orderP2), vmax=np.max(orderP2))
cb = fig.colorbar(p)
xlim = plt.xlim([np.min(x),np.max(x)])
ylim = plt.ylim([np.min(y),np.max(y)])

plt.xlabel('Membrane Node', fontsize=15)
plt.ylabel('Membrane Node', fontsize=15)
    
#plt.savefig(title+'.png')
#plt.savefig(title+'.eps', format='eps', dpi=1200)

#sns.set_palette("husl")

plt.show()
