#!/bin/bash

tprFile=TPRFILE
xtcFile=XTCFILE
pdbFile=PDBFILE
selection=""
target_selection=""
byFrame=True
zmin=-50
zmax=50
znum_bins=1
rmin=0
rmax=50
rnum_bins=50
title=""

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title
