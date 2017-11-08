#!/bin/bash

file=$1
file2=${file//.xvg}

ARRAY=(LEU-1261 GLU-1262 ARG-1263 GLN-1264 HIS-1265 THR-1266 PHE-1267)

for (( i=0; i<${#ARRAY[@]}; i++ )); do
	grep ${ARRAY[i]} $file > ${file2}_${ARRAY[i]}.xvg
	python PLOT_rama_hist.py ${file2}_${ARRAY[i]}.xvg
done

