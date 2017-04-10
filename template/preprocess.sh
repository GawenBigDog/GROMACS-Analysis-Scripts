#!/bin/bash

xtcFile=../data/1mff-first500ns_noPBC.xtc
tprFile=../data/topol_noWAT.tpr

gmx make_ndx -f $tprFile -o PO4.ndx <<EOF
a PO4
q
EOF

echo "PO4" | gmx trjconv -f $xtcFile -s $tprFile -n PO4.ndx -o PO4.xtc
echo "PO4" | gmx convert-tpr -s $tprFile -o PO4.tpr -n PO4.ndx
echo "0" | gmx trjconv -f PO4.xtc -o PO4.gro -s PO4.tpr -dt 9999999999
