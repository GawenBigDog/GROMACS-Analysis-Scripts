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
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <.dat file>\n")
        sys.exit()

data=sys.argv[1]

f0 = open(data,'r')
lines = f0.readlines()
f0.close()

nframes=int(len(lines)) 

distance = np.zeros(nframes)
pdf = np.zeros(nframes)

i=0
for line in lines:
    distance[i] = float(line.split()[0])
    pdf[i] = float(line.split()[2])
    i += 1

plt.bar(distance, pdf, color='k')
plt.xlabel(r"distance $r$ in nm")
plt.ylabel(r"probability")

plt.show()
