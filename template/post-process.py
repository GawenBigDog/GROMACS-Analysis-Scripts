#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
import scipy.interpolate

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <Filename> <# Nodes in one direction>\n")
        sys.exit()

# Defines files and arrays here
filename = sys.argv[1]
N = int(sys.argv[2])
nNodes = N * N

f0 = open(filename,'r')
lines = f0.readlines()[1:]

nFrames = 0

# Find the extreme values and the centroids of the system
for line in lines:
    if int(line.split()[1]) > nFrames:
        nFrames = int(line.split()[1])

nFrames += 1
nRes = int(len(lines)) / nFrames

f0.close()

# data[frame][residue][X, Y, Z, OrderP2]
data = [[[0,0,0,0] for i in range(nRes)] for j in range(nFrames)]

# renumber_array[frame][residue number]
renumber_array = [0 for i in range(nFrames)]

# centroids[frame][X, Y]
centroids = [[0,0,0] for i in range(nFrames)]

xmax = -np.inf
xmin = np.inf
ymax = -np.inf
ymin = np.inf

# Reread the file to file in the data array and renumber residues
for line in lines:
    frame = int(line.split()[1])
    data[frame][renumber_array[frame]][0] = 10 * float(line.split()[3])
    data[frame][renumber_array[frame]][1] = 10 * float(line.split()[4])
    data[frame][renumber_array[frame]][2] = 10 * float(line.split()[5])
    data[frame][renumber_array[frame]][3] = float(line.split()[6])
    centroids[frame][0] += data[frame][renumber_array[frame]][0]
    centroids[frame][1] += data[frame][renumber_array[frame]][1]
    centroids[frame][2] += data[frame][renumber_array[frame]][2]
    renumber_array[frame] += 1

# Center the plot at (0,0)
frame = 0
while frame < nFrames:
    centroids[frame][0] = centroids[frame][0] / nRes
    centroids[frame][1] = centroids[frame][1] / nRes
    centroids[frame][2] = centroids[frame][2] / nRes

    resnum=0
    while resnum < nRes:
        data[frame][resnum][0] = data[frame][resnum][0] - centroids[frame][0]
        data[frame][resnum][1] = data[frame][resnum][1] - centroids[frame][1]
        data[frame][resnum][2] = data[frame][resnum][2] - centroids[frame][2]
        x = data[frame][resnum][0]
        y = data[frame][resnum][1]
        if x > xmax:
            xmax = x
        if x < xmin:
            xmin = x
        if y > ymax:
            ymax = y
        if y < ymin:
            ymin = y
        resnum += 1
    frame += 1
    
# Create mesh grid here
xi, yi = np.linspace(xmin, xmax, N + int(N/10)), np.linspace(ymin, ymax, N + int(N/10))
xi = xi[int(N/20):-int(N/20)]
yi = yi[int(N/20):-int(N/20)]
xi, yi = np.meshgrid(xi, yi)

xi = np.reshape(xi,nNodes)
yi = np.reshape(yi,nNodes)

# *_node[frame][nodes][X, Y, orderP2, Count] 
lower_nodes = [[[0,0,0,0] for i in range(nNodes)] for j in range(nFrames)]
upper_nodes = [[[0,0,0,0] for i in range(nNodes)] for j in range(nFrames)]

frame = 0
while frame < nFrames:
    i = 0
    while i < nNodes:
        # Create the array of nodes here
        lower_nodes[frame][i][0] = xi[i]
        lower_nodes[frame][i][1] = yi[i]
        upper_nodes[frame][i][0] = xi[i]
        upper_nodes[frame][i][1] = yi[i]
        i += 1
    frame += 1

# Process the upper and lower leaflets
frame = 0
while frame < nFrames:
    print frame
    resnum = 0
    while resnum < nRes:
        mindist = np.inf
        minIndx = 0
        
        i = 0
        while (i < nNodes):
            dist = (lower_nodes[frame][i][0] - data[frame][resnum][0]) ** 2 \
                   + (lower_nodes[frame][i][1] - data[frame][resnum][1]) ** 2
            if (dist < mindist):
                mindist = dist
                minIndx = i
            i += 1
        if ( data[frame][resnum][2] < 0 ):
            lower_nodes[frame][minIndx][2] += data[frame][resnum][3]
            lower_nodes[frame][minIndx][3] += 1
        else:
            upper_nodes[frame][minIndx][2] += data[frame][resnum][3]
            upper_nodes[frame][minIndx][3] += 1

        resnum += 1
    frame += 1

# Average local values
frame = 0
while ( frame < nFrames ):
    f0 = open(str(frame)+'-upper.dat', 'w')
    f1 = open(str(frame)+'-lower.dat', 'w')
    f2 = open(str(frame)+'-average.dat', 'w')
    i = 0
    while ( i < nNodes):
        if (lower_nodes[frame][i][3] != 0):
            lower_nodes[frame][i][2] = lower_nodes[frame][i][2]/lower_nodes[frame][i][3]
        if (upper_nodes[frame][i][3] != 0):
            upper_nodes[frame][i][2] = upper_nodes[frame][i][2]/upper_nodes[frame][i][3]

        interp_lower_nodes = 0
        interp_upper_nodes = 0
        if lower_nodes[frame][i][2] == 0:
            count = 0
            if ( i + 1 ) % (N) != 0: # Checframe Right
                if lower_nodes[frame][i + 1][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i + 1][2]
                    count += 1
            else:
                if lower_nodes[frame][i - N + 1][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i - N + 1][2]
                    count += 1
            if i == 1 or i == N * N - N + 1 or ( ( i - 1 ) % (N - 1) != 0 and i != 0): # Check Left
                if lower_nodes[frame][i - 1][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i - 1][2]
                    count += 1
            else:
                if lower_nodes[frame][i + N - 1][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i + N - 1][2]
                    count += 1
            if ( i + N ) < nNodes: # Checframe Up
                if lower_nodes[frame][i + N][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i + N][2]
                    count += 1
            else:
                if lower_nodes[frame][i + N - nNodes][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i + N - nNodes][2]
                    count += 1
            if ( i - N ) >= 0: # Check Down
                if lower_nodes[frame][i - N][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i - N][2]
                    count += 1
            else:
                if lower_nodes[frame][i - N + nNodes][2] != 0:
                    interp_lower_nodes += lower_nodes[frame][i - N + nNodes][2]
                    count += 1
            if count != 0:
                interp_lower_nodes = interp_lower_nodes / count
        else:
            interp_lower_nodes = lower_nodes[frame][i][2]

        if upper_nodes[frame][i][2] == 0:
            count = 0
            if ( i + 1 ) % (N) != 0: # Check Right
                if upper_nodes[frame][i + 1][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i + 1][2]
                    count += 1
            else:
                if upper_nodes[frame][i - N + 1][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i - N + 1][2]
                    count += 1
            if i == 1 or i == N * N - N + 1 or ( ( i - 1 ) % (N - 1) != 0 and i != 0): # Check Left
                if upper_nodes[frame][i - 1][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i - 1][2]
                    count += 1
            else:
                if upper_nodes[frame][i + N - 1][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i + N - 1][2]
                    count += 1
            if ( i + N ) < nNodes: # Check Up
                if upper_nodes[frame][i + N][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i + N][2]
                    count += 1
            else:
                if upper_nodes[frame][i + N - nNodes][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i + N - nNodes][2]
                    count += 1
            if ( i - N ) >= 0: # Check Down
                if upper_nodes[frame][i - N][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i - N][2]
                    count += 1
            else:
                if upper_nodes[frame][i - N + nNodes][2] != 0:
                    interp_upper_nodes += upper_nodes[frame][i - N + nNodes][2]
                    count += 1
            if count!= 0:
                interp_upper_nodes = interp_upper_nodes / count
        else:
            interp_upper_nodes = upper_nodes[frame][i][2]
            
        f0.write('{:6.2f} {:10.4f} {:10.4f} {:16.8f}\n'
                 .format(i, lower_nodes[frame][i][0], lower_nodes[frame][i][1], interp_lower_nodes))
        f1.write('{:6.2f} {:10.4f} {:10.4f} {:16.8f}\n'
                 .format(i, upper_nodes[frame][i][0], upper_nodes[frame][i][1], interp_upper_nodes))
        average_node = ( interp_upper_nodes + interp_lower_nodes ) / 2
        f2.write('{:6.2f} {:10.4f} {:10.4f} {:16.8f}\n'
                 .format(i, upper_nodes[frame][i][0], upper_nodes[frame][i][1], average_node))
        i += 1
    f0.close()
    f1.close()
    f2.close()
    frame += 1
                                
