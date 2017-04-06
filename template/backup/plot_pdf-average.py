#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: %s <filename> <total # of files> <nBins>\n")
        sys.exit()

baseName=sys.argv[1]
totalCount=int(sys.argv[2])
nBins=int(sys.argv[3])

distance = np.zeros(int(nBins))
pdf = np.zeros(int(nBins))

i=0
while i <= totalCount:
    f0 = open(str(baseName)+str(i)+'.dat','r')
    lines = f0.readlines()
    f0.close()

    j=0
    while j < nBins - 1:
        distance[j] += float(lines[0].split()[j])
        pdf[j] += float(lines[1].split()[j])
        j += 1
    i += 1

distance /= totalCount
pdf /= totalCount
    
plt.bar(distance, pdf, color='k')
plt.xlabel(r"distance $r$ in $\AA$")
plt.ylabel(r"probability")

plt.show()
