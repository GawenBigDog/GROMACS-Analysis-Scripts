#!/bin/bash

i=0
imax=56

lipidnames="PIPC POPC DPPC PAPC PEPC PGPC PUPC PRPC PIPE POPE DPPE PAPE PQPE PGPE PUPE PRPE POPA PGPA PGPS PRPS POPI PVPI PAPI PPC POSM PVSM DPSM PGSM PNSM DPCE"

rm tmp tmp2
rm orderParameter_all.dat

for lipid in $lipidnames; do
    cat order-$lipid.dat >> tmp
done

grep -v "{" tmp > tmp2

echo "   {LIPID}    {FRAME}   {RESNUM}        {X}        {Y}        {Z} {ORDERPARAM}" > orderParameter_all.dat

cat tmp2 >> orderParameter_all.dat

rm tmp tmp2
	  
