#!/bin/bash

XTC=$1
TPR=$2
XTC2=${XTC%.*}

echo -e Protein '\n' System | gmx_sse trjconv -f ${XTC2}.xtc -s ${TPR} -center -o ${XTC2}.1.xtc -pbc nojump >& trjone
echo trj1 done
echo System | gmx_sse trjconv -f ${XTC2}.1.xtc -s ${TPR} -o ${XTC2}.pbc.xtc -pbc mol >& trjtwo
echo trj2 done
rm -f ${XTC2}.1.xtc
