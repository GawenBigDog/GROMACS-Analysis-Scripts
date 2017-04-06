#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <.dat file> <nBins>\n")
        sys.exit()

data=sys.argv[1]
nBins=int(sys.argv[2])

f0 = open(data,'r')
lines = f0.readlines()
f0.close()

distance = np.zeros(int(nBins))
pdf = np.zeros(int(nBins))

print distance

i=0
while i < nBins- 1:
    print "Distance: "+str(lines[0].split()[i])
    print "PDF: "+str(lines[1].split()[i])
    distance[i] = float(lines[0].split()[i])
    pdf[i] = float(lines[1].split()[i])
    i += 1

plt.bar(distance, pdf, color='k')
plt.xlabel(r"distance $r$ in nm")
plt.ylabel(r"probability")

plt.show()
