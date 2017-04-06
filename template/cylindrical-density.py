#!/usr/bin/env python
"""
cylindrical-density.py

Details : compute the density of atoms ("targets") around a centering selection
Use-case: Lipid packing around a membrane protein
Output  : a 2D histogram in r,z 

Adapted from
Lightweight Object-Oriented Structure library (LOOS)
Tod D. Romo and Alan Grossfield
Department of Biochemistry and Biophysics
School of Medicine & Dentistry, University of Rochester
http://grossfieldlab.github.io/loos/
"""

import sys
import numpy as np

from MDAnalysis import *
from MDAnalysis.analysis.align import *
from MDAnalysis.tests.datafiles import TPR, XTC

header = " ".join(sys.argv)
print "# ", header

tprFile = sys.argv[1]
xtcFile = sys.argv[2]
pdbFile = sys.argv[3]

selection = sys.argv[4]
target_selection = sys.argv[5]

byFrameCheck = sys.argv[6]
if (byFrameCheck == 'True'):
    byFrame = True
else:
    byFrame = False
    
zmin = float(sys.argv[7])
zmax = float(sys.argv[8])
znum_bins = int(sys.argv[9])
rmin = float(sys.argv[10])
rmax = float(sys.argv[11])
rnum_bins = int(sys.argv[12])
title = sys.argv[13]

universe = Universe(tprFile, xtcFile)
ref = Universe(tprFile, pdbFile)

rms_fit_trj(universe, ref, filename='fit.xtc', select=selection)

fit_u = Universe(tprFile, 'fit.xtc')

zbin_width = (zmax - zmin) / znum_bins

rbin_width = (rmax - rmin) / rnum_bins
rmin2 = rmin*rmin
rmax2 = rmax*rmax

hist = np.zeros([rnum_bins, znum_bins])

nframes = 0

for timestep in fit_u.trajectory:

    print timestep.time
    nframes += 1

    if (byFrame):
        hist = np.zeros([rnum_bins, znum_bins])
    
    centroid = fit_u.select_atoms(selection).center_of_geometry()

    target = fit_u.select_atoms(target_selection)

    for atom in target.atoms:
        x = atom.position[0] - centroid[0]
        y = atom.position[1] - centroid[1]
        z = atom.position[2] - centroid[2]

        r2 = x*x + y*y

        if (zmin < z < zmax) and (rmin2 < r2 < rmax2):
            r = np.sqrt(r2)

            rbin = int((r - rmin)/ rbin_width)
            zbin = int((z - zmin)/ zbin_width)
            
            hist[rbin, zbin] += 1.0
    if(byFrame):
        f0 = open('pdf-'+str(title)+'-'+str(nframes)+'.dat', 'w') 
        for i in range(rnum_bins):
            rinner = rmin + i*rbin_width
            router = rinner + rbin_width
            rval = rmin + (i+0.5)*rbin_width
            norm = np.pi * (router*router - rinner*rinner)
            for j in range(znum_bins):
                zval = zmin + (j+0.5)*zbin_width
                f0.write('{:5.2f} {:5.2f} {:10.8f}\n'.format(rval, zval, hist[i,j] / norm))
        f0.close()

        
if not byFrame:
    hist /= nframes
    f0 = open('pdf-'+str(title)+'-all.dat', 'w')
    for i in range(rnum_bins):
        rinner = rmin + i*rbin_width
        router = rinner + rbin_width
        rval = rmin + (i+0.5)*rbin_width
        norm = np.pi * (router*router - rinner*rinner)
        for j in range(znum_bins):
            zval = zmin + (j+0.5)*zbin_width
            f0.write('{:5.2f} {:5.2f} {:10.8f}\n'.format(rval, zval, hist[i,j] / norm))
    f0.close()

