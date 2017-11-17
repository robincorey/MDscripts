#!/bin/bash

NUC=(ATP ADP)

CD=`pwd`

for (( f=0; f<${#NUC[@]}; f++ )); do
	rm -f ${NUC[f]}_rama_half*.xvg 
	cd $CD
	mkdir -p ${NUC[f]}_rama
	cd ${NUC[f]}_rama
	XTC="../${NUC[f]}/${NUC[f]}.pbc.xtc"
	TPR="../${NUC[f]}/md_tom_a*p.tpr"
	echo -e r1260-1268 '\n' r1269-1277 '\n' q  | $GMX51 make_ndx -f $TPR -o subs.ndx  >& out2
	sed 's/r_1260-1268/half1/g' subs.ndx -i
	sed 's/r_1269-1277/half2/g' subs.ndx -i
	# find contacts
	HALF=(half1 half2)
	for (( k=0; k<${#HALF[@]}; k++ )); do
		TIMES=(0 1 2 3 4 5 10 15 20 30 50 75 100 200 300 400 500 600 700 800 900 1000)
		for (( j=0; j<${#TIMES[@]}; j++ )); do
			time1=`echo "${TIMES[j]} * 1000" | bc`
			time2=`echo "(${TIMES[j]} * 1000) + 1000" | bc`
			timeac=`echo "(${TIMES[j]} * 1000) + 500" | bc`
			$GMX51 rama -s ${TPR} -f ${XTC} -o rama_${TIMES[j]}_${HALF[k]}.xvg >& rama_${TIMES[j]}_${HALF[k]} 
			forbidden=`grep "^6\|^7\|^8\|^9\|^10" rama_${TIMES[j]}_${HALF[k]}.xvg | awk '{print $2}' | grep "^10\|^11\|^12\|^13\|^14\|^15\|^16\|^17" | wc -l`	
			echo ${timeac} $forbidden >> ${NUC[f]}_rama_${HALF[k]}.xvg
		echo 
		done
	done
done

