#!/bin/bash

rm -f *plot *err

for i in {1..5}
do
	for j in {200..200..25}
	do 
		grep total ${j}_*out*_$i | awk '{print $6}' >> ${i}.plot
		grep total ${j}_*out*_$i | awk '{print $8}' >> ${i}.err	
	done
done

paste -d ' ' *plot > plot.xvg
paste -d ' ' *err > err.xvg

python PLOT_FEP_data.py
