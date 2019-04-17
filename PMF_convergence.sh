#!/bin/bash

for i in `seq 50000 50000 800000`
do
        j=`echo "$i+200000" | bc`
	echo $j
	mkdir -p convergence
        if [[ ! -f convergence/${i}_pmf.xvg ]]
	then
		gmx wham -it t.dat -if f.dat -o convergence/${i}_data.xvg -hist convergence/${i}_hist.xvg -temp 323 -b 200000 -e $j -bsprof convergence/${i}_bsprof.xvg -nBootstrap 200 -bsres convergence/${i}_pmf.xvg >& convergence/out_${i} 
	fi
done

$CG/PLOT_PMF_equil.py > equil_data.txt
$CG/PLOT_equil_timecourse.py 
