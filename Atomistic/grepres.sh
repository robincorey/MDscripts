#!/bin/bash

file=$1
file2=${file//.xvg}

ARRAY=(PHE-1270 THR-1271 HIS-1272 GLN-1273 ARG-1274 GLU-1275 LEU-1276)

for (( i=0; i<${#ARRAY[@]}; i++ )); do
	grep ${ARRAY[i]} $file > ${file2}_${ARRAY[i]}.xvg
	python PLOT_rama_hist.py ${file2}_${ARRAY[i]}.xvg
done

