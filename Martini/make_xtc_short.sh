#!/bin/bash

XTC=$1
TPR=$2

xtc_name=`echo $XTC | cut -f1 -d'.'`

echo 0 | gmx_sse trjconv -f ${XTC} -s ${TPR} -skip 100 -o ${xtc_name}.pbc.xtc -pbc mol -nice -19
