#!/bin/bash

NUC=(ATP) # )ADP)

CD=`pwd`

for (( f=0; f<${#NUC[@]}; f++ )); do 
	cd $CD
	cd ${NUC[f]}
	if [ ${NUC[f]} == 'ATP' ]; then
		TIMES=(70 80 110)
	elif [ ${NUC[f]} == 'ADP' ]; then
		TIMES=(50 65 80)
	else
		echo NUC loop failed
		exit 0
	fi

for i in {0..4}; do
	XTC="../${NUC[f]}/RQH_full_Data_Archer/RQH_full_${TIMES[2]}/RQH_full_${TIMES[2]}_${i}.xtc"
	TPR="../${NUC[f]}/RQH_full_Data_Archer/RQH_full_${TIMES[2]}/RQH_full_${TIMES[2]}_${i}.tpr"
	# define y
	echo -e r757-1174 '\n' q | $GMX51 make_ndx -f $TPR -o Y.ndx >& out1
	# deine subs
	echo -e r1260-1268 '\n' r1269-1277 '\n' q  | $GMX51 make_ndx -f $TPR -o Y_subs.ndx -n Y.ndx >& out2
	sed 's/r_757-1174/SecY/g' Y_subs.ndx -i
	sed 's/r_1260-1268/half1/g' Y_subs.ndx -i
	sed 's/r_1269-1277/half2/g' Y_subs.ndx -i
	# find contacts
	HALF=(half1 half2)
	for (( k=0; k<${#HALF[@]}; k++ )); do
		echo "group ${HALF[k]} and within 0.5 of group SecY" > gselect_${HALF[k]}_${TIMES[2]}_${i}
		$GMX51 "select" -f ${XTC} -s ${TPR} -sf gselect_${HALF[k]}_${TIMES[2]}_${i} -on contact_${HALF[k]}_${TIMES[2]}_${i} -selrpos whole_res_com -rmpbc -n Y_subs -b 50000 -e 100000 >& out_${HALF[k]}_${TIMES[2]}_${i}
#		TIMES=(50 55 60 65 70 85 90 95 100)
#		for (( j=0; j<${#TIMES[@]}; j++ )); do
#			time1=`echo "${TIMES[j]} * 1000" | bc`
#			time2=`echo "(${TIMES[j]} * 1000) + 5000" | bc`
#			$GMX51 "select" -f ${XTC} -s ${TPR} -sf gselect_${HALF[k]} -on contact_${j}_${HALF[k]} -selrpos whole_res_com -rmpbc -n Y_subs -b $time1 -e $time1 >& out_a_${j}_${HALF[k]}
		done
	done
done
