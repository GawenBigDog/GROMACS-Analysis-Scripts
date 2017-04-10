#!/bin/bash

xtcFile=../data/1mff-last500ns_noPBC.xtc
tprFile=../data/topol_noWAT.tpr

./filter.sh $xtcFile $tprFile 1

python leaflet_finder.py ../data/ref_noPBC.pdb "name PO4" l0 l1

gmx trjconv -f $xtcFile -s $tprFile -n l0.ndx -o l0.xtc

gmx trjconv -f $xtcFile -s $tprFile -n l1.ndx -o l1.xtc

gmx trjconv -f $xtcFile -s $tprFile -o prot.xtc

gmx make_ndx -f $tprFile -o PO4.ndx

gmx trjconv -f $xtcFile -s $tprFile -n PO4.ndx -o PO4.xtc
gmx convert-tpr -s $tprFile -o PO4.tpr -n PO4.ndx
echo "0" | gmx trjconv -f PO4.xtc -o PO4.gro -s PO4.tpr -dt 9999999999
