#!/bin/bash

### BEGIN ADJUSTABLE PARAMETERS ###

# General name of the system to analyze
allAtomName="0mff-first500ns"
filename=${allAtomName}"_noPBC"
ndxName="0mff+gsec+app"

start=0     # Start time (ps) for order parameter analysis
end=5000000 # Stop time (ps) for order parameter analysis

protGrp=13 # Group specification for Protein
maxProt=6  # Number of proteins in the system

nFrames=50 # Number of frames for the thickness analysis 

### END ADJUSTABLE PARAMETERS ###

# General directory specification
dataDirectory="../data/"
template="template/"

dDirectory="densmap/"
mDirectory="mindist/"
oDirectory="order/"
pDirectory="pdf/"
tDirectory="thickness/"

# General trajectory file specification with no solvent
xtcFile=${filename}".xtc"
tprFile="topol.tpr"
tprFileAlt="topol_noWAT.tpr"
groFile=${filename}".gro"
pdbFile=${filename}".pdb"
ndxFile=${ndxName}".ndx"
refFile="ref_noPBC.pdb"

# General trajectory file specification with solvent
allAtomPDB=${allAtomName}".pdb"

# File extension specification for file conversion
inFormat="xtc"
outFormat="pdb"

# Generate all analysis scripts
sed -e "s~FILENAME~${allAtomName}~" ${template}template.convert_xtc.sh | \
    sed -e "s~INFORMAT~$inFormat~" | \
    sed -e "s~OUTFORMAT~$outFormat~" | \
    sed -e "s~NDXFILE~${ndxFile}~" \
    > data/convert_xtc.sh

sed -e "s~XTCFILE~${dataDirectory}${xtcFile}~" ${template}template.densmap.sh | \
    sed -e "s~NDXFILE~${dataDirectory}${ndxFile}~" \
    > ${dDirectory}densmap.sh

sed -e "s~FILENAME~${allAtomName}~" ${template}template.merge.sh | \
    sed -e "s~NDXFILE~${ndxFile}~" \
    > data/merge.sh

sed -e "s~XTCFILE~${dataDirectory}${xtcFile}~" ${template}template.mindist.sh | \
    sed -e "s~TPRFILE~${dataDirectory}${tprFileAlt}~" | \
    sed -e "s~PROTGRP~${protGrp}~" | \
    sed -e "s~MAXPROT~${maxProt}~" \
    > ${mDirectory}mindist.sh

sed -e "s~XTCFILE~${dataDirectory}${xtcFile}~" ${template}template.batchOrder.sh | \
    sed -e "s~TPRFILE~${dataDirectory}${tprFile}~" | \
    sed -e "s~STARTTIME~${start}~" | \
    sed -e "s~ENDTIME~${end}~" \
    > ${oDirectory}batchOrder.sh

sed -e "s~NFRAMES~${nFrames}~" ${template}template.gmat_thickness.in | \
    sed -e "s~PDBFILE~${dataDirectory}${allAtomPDB}~" \
    > ${tDirectory}gmat_thickness.in

sed -e "s~XTCFILE~${dataDirectory}${xtcFile}~" ${template}template.run-pdf.sh | \
    sed -e "s~TPRFILE~${dataDirectory}${tprFileAlt}~" | \
    sed -e "s~PDBFILE~${dataDirectory}${refFile}~" \
    > ${pDirectory}run-pdf.sh

chmod +x data/convert_xtc.sh
chmod +x ${dDirectory}densmap.sh
chmod +x data/merge.sh
chmod +x ${mDirectory}mindist.sh
chmod +x ${oDirectory}batchOrder.sh
chmod +x ${pDirectory}run-pdf.sh

echo
echo "To properly run all the scripts, first run the following: "
echo "     (1) data/merge.sh"
echo "     (2) data/convert_xtc.sh"
echo
