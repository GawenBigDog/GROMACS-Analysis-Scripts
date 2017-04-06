#!/usr/bin/env python

import sys, os
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.patches as ptch
import scipy
import scipy.interpolate
#import seaborn

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <.dat file> <frame>\n")
        sys.exit()

datafile=sys.argv[1]

try:
    sys.argv[2]
except IndexError:
    frame = 0
else:
    frame = sys.argv[2]

f = open(datafile,'r')
lines = f.readlines()[1:]

x = []
y = []
orderP2 = []

for line in lines:
    if float(line.split()[1]) == frame:
        x.append(float(line.split()[3]))
        y.append(float(line.split()[4]))
        orderP2.append(float(line.split()[6]))
f.close()

fig = plt.figure()

# Set up a regular grid of interpolation points
xi, yi = np.linspace(np.min(x), np.max(x), 100), np.linspace(np.min(y), np.max(y), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpolate
rbf = scipy.interpolate.Rbf(x, y, orderP2, function='linear')
orderP2_i = rbf(xi, yi)

plt.imshow(orderP2_i, vmin=0, vmax=1, origin='lower',
                      extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
plt.scatter(x, y, c=orderP2)
plt.colorbar()
plt.show()
