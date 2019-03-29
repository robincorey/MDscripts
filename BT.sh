#!/bin/bash

XTC=$1
TPR=$2
P=$3

name=${TPR::-4}

echo -e a$P '\n' q | gmx_sse make_ndx -f ${TPR} -o P.ndx
echo $P | gmx_sse density -f ${XTC} -s ${TPR} -n P.ndx -o BT_$name -d Z -sl 200 -nice -19

peak1=`cat BT_dynamic.xvg | grep -v \# | grep -v \@ | sort -k2 -h | awk '{print $1}' | tail -n 1`
peak2=`cat BT_dynamic.xvg | grep -v \# | grep -v \@ | sort -k2 -h | awk '{print $1}' | tail -n 2 | head -n 1`

diff=`echo "scale=4; $peak1 - $peak2" | bc`
if [[ $diff < 2 ]]; then
        peak2=`cat BT_dynamic.xvg | grep -v \# | grep -v \@ | sort -k2 -h | awk '{print $1}' | tail -n 3 | head -n 1`
        diff=`echo "scale=4; $peak1 - $peak2" | bc`
fi

echo BT is $diff > BT.value.xvg

$CG/PLOT_BT.py  BT_$name.xvg

