#!/bin/bash

xtcFile=XTCFILE
tprFile=TPRFILE
start=STARTTIME
end=ENDTIME

python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PIPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 POPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 DPPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PAPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PEPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PGPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PUPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PRPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PIPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 POPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 DPPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PAPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PQPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PGPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PUPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PRPE
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 POPA
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PGPA
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PGPS
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PRPS
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 POPI
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PVPI
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PAPI
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1  PPC
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 POSM
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PVSM
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 DPSM
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PGSM
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 PNSM
python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 DPCE
#python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 DPG3
#python order_map.py $xtcFile $tprFile $start $end 1 0 0 1 CHOL

rm *#
