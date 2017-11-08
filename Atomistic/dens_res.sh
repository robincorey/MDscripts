#!/bin/bash

ARRAY=(PHE-1270 THR-1271 HIS-1272 GLN-1273 ARG-1274 GLU-1275 LEU-1276)
ARRAY2=(20333 20353 20367 20384 20401 20425 20440)

for (( i=0; i<${#ARRAY[@]}; i++ )); do
	num=`grep -c ${ARRAY2[i]} $1`
	echo ${ARRAY[i]} $num 
done

