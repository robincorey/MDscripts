#!/bin/bash

#ARRAY=(192.small 192_reseed.1.small 192_reseed.2.small 192_reseed.3.small 192_reseed.4.small unbiased_BC3)
ARRAY=(all_CG)
ARRAYLEN=`echo ${#ARRAY[@]} - 1 | bc`
CURRENTDIR=`pwd`
TPR=unbiased.tpr
GRO=equilibration_2.gro

for i in `seq 0 $ARRAYLEN`; do

#echo -e a1843-2983 '\n' q | make_ndx -f ${GRO} -o YEG.ndx
# make protein whole
echo System | $GMX5/trjconv -s ${TPR} -f ${ARRAY[i]}.xtc -o ${ARRAY[i]}.1.xtc -pbc nojump -n YEG.ndx -nice -19
# center on prot
echo -e a_1843-2983 '\n' System | $GMX5/trjconv -s ${TPR} -f ${ARRAY[i]}.1.xtc -o ${ARRAY[i]}.2.xtc -center -n YEG.ndx -nice -19
# fit trans
echo -e a_1843-2983'\n' System | $GMX5/trjconv -s ${TPR} -f ${ARRAY[i]}.2.xtc -fit rotxy+trans -o ${ARRAY[i]}.3.xtc -n YEG.ndx -nice -19
# pbc to remake sim box
echo System | $GMX5/trjconv -s ${TPR} -f ${ARRAY[i]}.3.xtc  -o ${ARRAY[i]}.final.xtc -pbc mol -n YEG.ndx -nice -19
# make a nice input
#echo System | $GMX5/editconf -f ${TPR} -o ${ARRAY[i]}_input.gro -n YEG.ndx
if [[ -f densmap.xpm ]]; then
	rm densmap.xpm
fi 
echo a_1843-2983 | $GMX5/gmx densmap -f ${ARRAY[i]}.final.xtc -s ${TPR} -n YEG.ndx
mv densmap.xpm ${ARRAY[i]}_protein_densmap.xpm
echo PO1_PO2_GL0 | $GMX5/gmx densmap -f ${ARRAY[i]}.final.xtc -s ${TPR} -n CL.ndx
mv densmap.xpm ${ARRAY[i]}_CL_densmap.xpm
#$GMX5/gmx xpm2ps -f 192.small.trimmed_CL_densmap.xpm -f2 192.small.trimmed_protein_densmap.xpm -o combined -rainbow red -xpm combined -diag none -combine add
$GMX5/gmx xpm2ps -f ${ARRAY[i]}_CL_densmap.xpm -o ${ARRAY[i]}_CL_densmap.xpm -xpm ${ARRAY[i]}_CL_densmap.xpm
done
