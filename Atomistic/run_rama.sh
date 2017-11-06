#!/bin/bash

for i in {0..0}; do
	XTC="../../ATP/RQH_full_Data_Archer/RQH_full_110/RQH_full_110_${i}.pbc.xtc"
	TPR="../../ATP/RQH_full_Data_Archer/RQH_full_110/RQH_full_110_${i}.tpr"
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
		echo "group ${HALF[k]} and within 0.5 of group SecY" > gselect_${HALF[k]}
		TIMES=(50 55 60 65 70 85 90 95 100)
		for (( j=0; j<${#TIMES[@]}; j++ )); do
			time1=`echo "${TIMES[j]} * 1000" | bc`
			time2=`echo "(${TIMES[j]} * 1000) + 5000" | bc`
			$GMX51 "select" -f ${XTC} -s ${TPR} -sf gselect_${HALF[k]} -on contact_${j}_${HALF[k]} -selrpos whole_res_com -rmpbc -n Y_subs -b $time1 -e $time1 >& out_a_${j}_${HALF[k]}
			$GMX51 trjconv -f ${XTC} -s ${TPR} -n contact_${j}_${HALF[k]}.ndx -o contact_${j}_${HALF[k]}.xtc -b $time1 -e $time2 >& out_b_${j}_${HALF[k]}
			$GMX51 rama -s ${TPR} -f contact_${j}_${HALF[k]}.xtc -o rama_${j}_${HALF[k]}.xvg >& rama_${j}_${HALF[k]} &
		done
	done
done
