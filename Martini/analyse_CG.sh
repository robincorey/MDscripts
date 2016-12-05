#!/bin/bash
# analyse CG stuff
# ugly as fuck
CURRENT_DIR=`pwd`
GMX5=/usr/local/gromacs_5/bin
ARRAY=(unbiased/)

#echo -e Protein'\n'Protein | g_rms -s unbiased.tpr -f unbiased.pbc.xtc -o CG_7000.xvg >& ca_rms &

echo -e aPO1 "|" aPO2 '\n' q | $GMX5/make_ndx -f  ${ARRAY[i]%/}.tpr -o  ${ARRAY[i]%/}_phos.ndx
$GMX5/g_select -sf selection_0.5.head.dat -f ${ARRAY[i]%/}.pbc.xtc  -s ${ARRAY[i]%/}.tpr -selrpos atom -os ../${ARRAY[i]%/}_0.5_head.xvg -n  ${ARRAY[i]%/}_phos.ndx >& gselhead &
$GMX5/g_select -sf selection_1.head.dat -f ${ARRAY[i]%/}.pbc.xtc  -s ${ARRAY[i]%/}.tpr -selrpos atom -os ../${ARRAY[i]%/}_1_head.xvg -n  ${ARRAY[i]%/}_phos.ndx >& gsel_1_head &
$GMX5/g_select -sf selection_2.head.dat -f ${ARRAY[i]%/}.pbc.xtc  -s ${ARRAY[i]%/}.tpr -selrpos atom -os ../${ARRAY[i]%/}_2_head.xvg -n  ${ARRAY[i]%/}_phos.ndx >& gsel_2_head &
grep "0.000" ${ARRAY[i]%/}_2_head.xvg | awk '{print $2}' > ${ARRAY[i]%/}_2head.num
sed "1 i\\${ARRAY[i]%/}" ${ARRAY[i]%/}_2head.num -i
grep "0.000" ${ARRAY[i]%/}_1_head.xvg | awk '{print $2}' > ${ARRAY[i]%/}_1head.num
sed "1 i\\${ARRAY[i]%/}" ${ARRAY[i]%/}_1head.num -i
grep "0.000" ${ARRAY[i]%/}_0.5_head.xvg | awk '{print $2}' > ${ARRAY[i]%/}_0.5head.num
sed "1 i\\${ARRAY[i]%/}" ${ARRAY[i]%/}_0.5head.num -i
paste -d , *_0.5head.num > 0.5.head.csv
paste -d , *_1head.num > 1.head.csv
paste -d , *_2head.num > 2.head.csv


