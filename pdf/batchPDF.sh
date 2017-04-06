#!/bin/bash

tprFile=../data/topol_noWAT.tpr
xtcFile=../data/0mff-first500ns_noPBC.xtc
pdbFile=../data/ref_noPBC.pdb
selection="name BB and bynum 0:3683"
target_selection="name ROH"
byFrame=False
zmin=-50
zmax=50
znum_bins=1
rmin=0
rmax=50
rnum_bins=50
title="gsec"

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title

selection="name BB and bynum 3684:3791"
title="app1"

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title

selection="name BB and bynum 3792:3899"
title="app2"

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title

selection="name BB and bynum 3900:4007"
title="app3"

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title

selection="name BB and bynum 4008:4115"
title="app4"

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title

selection="name BB and bynum 4116:4223"
title="app5"

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title

selection="name BB and bynum 4224:4331"
title="app6"

python cylindrical-density.py $tprFile $xtcFile $pdbFile "$selection" "$target_selection" $byFrame $zmin $zmax $znum_bins $rmin $rmax $rnum_bins $title
