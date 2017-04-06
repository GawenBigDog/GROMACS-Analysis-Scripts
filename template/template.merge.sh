#!/bin/bash

filename=FILENAME
ndxFile=NDXFILE

gmx trjcat -f *.xtc -o ${filename}.xtc

gmx trjconv -s topol.tpr -f ${filename}.xtc -o ${filename}_noPBC.xtc -pbc mol -ur compact -n ${ndxFile}
