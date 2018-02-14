#!/bin/bash

XTC=$1
TPR=$2
XTC2=${XTC%.*}

#echo -e Protein '\n' System | gmx_sse trjconv -f ${XTC2}.xtc -s ${TPR} -center -o ${XTC2}.1.xtc -pbc nojump -nice -5 >& trjone
#echo trj1 done
#echo System | gmx_sse trjconv -f ${XTC2}.1.xtc -s ${TPR} -o ${XTC2}.pbc.xtc -pbc mol -nice -5 >& trjtwo
#echo trj2 done
#rm -f ${XTC2}.1.xtc

echo -e Protein '\n' System | gmx_sse trjconv -f ${XTC} -s ${TPR} -o ${XTC2}.pbc.xtc -center -pbc res -nice -5 -skip 500

