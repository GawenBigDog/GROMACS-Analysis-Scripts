#!/bin/bash

filename=FILENAME
inFormat=INFORMAT
outFormat=OUTFORMAT
ndxFile=NDXFILE

gmx trjconv -f $filename.$inFormat -s topol.tpr -o $filename.$outFormat -dt 10000
gmx trjconv -f $filename.$inFormat -s topol.tpr -o ref.$outFormat -dt 999999999

gmx convert-tpr -s topol.tpr -o topol_noWAT.tpr -n $ndxFile 
gmx trjconv -f ${filename}_noPBC.xtc -s topol_noWAT.tpr -o ref_noPBC.$outFormat -dt 999999999