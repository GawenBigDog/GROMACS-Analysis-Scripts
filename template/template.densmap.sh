#!/bin/bash

xtcFile=XTCFILE
ndxFile=NDXFILE

gmx densmap -f $xtcFile -n $ndxFile

gmx xpm2ps -f densmap.xpm -o densmap.eps -rainbow blue #-gradient 1 0 0
