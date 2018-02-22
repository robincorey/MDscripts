#!/bin/bash

XTC=$1
TPR=$2
XTC2=${XTC%.*}

echo -e aBB '\n' q | gmx_avx make_ndx -f ${TPR} -o BB.ndx
echo -e Protein '\n' System | gmx_avx trjconv -f ${XTC2}.gro -s ${TPR} -o ${XTC2}.pdb -center -pbc res 
#echo -e BB '\n' BB | gmx rms -f md-center.xtc -s md.tpr -o rmsd.xvg 
