#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <filename> <total # of files>\n")
        sys.exit()

baseName=sys.argv[1]
totalCount=int(sys.argv[2])

i=1
while i <= totalCount:
    f0 = open(str(baseName)+str(i)+'.dat','r')
    lines = f0.readlines()
    f0.close()

    if (i == 1):
        nframes=int(len(lines))
        distance = np.zeros(nframes)
        pdf = np.zeros(nframes)
        
    j=0
    for line in lines:
        distance[j] += float(line.split()[0])
        pdf[j] += float(line.split()[2])
        j += 1
    i += 1

distance /= totalCount
pdf /= totalCount
    
plt.bar(distance, pdf, color='k')
plt.xlabel(r"distance $r$ in $\AA$")
plt.ylabel(r"probability")

plt.show()
