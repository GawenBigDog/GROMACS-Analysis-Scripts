#!/bin/bash

# Specify Parameters here
ndxFile=./mindist.ndx
xtcFile=XTCFILE
tprFile=TPRFILE
proteinGroup=PROTGRP
imax=MAXPROT

# Create a new index file with the appropriate groups
gmx make_ndx -f $tprFile -o $ndxFile <<EOF
del 13-81
a 1-3672
a 3673-3780
a 3781-3888
a 3889-3996
a 3999-4104
a 4105-4212
a 4213-4320
q
EOF

# Perform minimum distance calculation with mindist Gromacs analysis tool
i=1
while [ $i -le $imax ]; do
    compGroup=`expr $proteinGroup + $i` 
    gmx mindist -f $xtcFile -s $tprFile -n $ndxFile -od mindist.xvg <<EOF
$proteinGroup
$compGroup
EOF
    
    mv mindist.xvg app${i}_mindist.xvg

    i=`expr $i + 1`

done
