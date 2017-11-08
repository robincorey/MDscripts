#!/bin/bash

NUC=(ATP ADP)

CD=`pwd`

for (( f=0; f<${#NUC[@]}; f++ )); do 
	cd $CD
	mkdir -p ${NUC[f]}_distance
	cd ${NUC[f]}_distance
	if [ ${NUC[f]} == 'ATP' ]; then
		TIMES=(70 80 110)
	elif [ ${NUC[f]} == 'ADP' ]; then
		TIMES=(50 65 80)
	else
		echo NUC loop failed
		exit 0
	fi
	for (( j=0; j<${#TIMES[@]}; j++ )); do
		for i in {0..4}; do
			if [ ${NUC[f]} == 'ATP' ]; then
				XTC="../../${NUC[f]}/RQH_full_Data_Archer/RQH_full_${TIMES[j]}/RQH_full_${TIMES[j]}_${i}.xtc"
				TPR="../../${NUC[f]}/RQH_full_Data_Archer/RQH_full_${TIMES[j]}/RQH_full_${TIMES[j]}_${i}.tpr"
			 elif [ ${NUC[f]} == 'ADP' ]; then
				XTC="../../${NUC[f]}/RQH_full_Data_BC3_Ian/RQH_full_${NUC[f]}_${TIMES[j]}/RQH_full_${NUC[f]}_${TIMES[j]}_${i}.secondpart.xtc"
				TPR="../../${NUC[f]}/RQH_full_Data_BC3_Ian/RQH_full_${NUC[f]}_${TIMES[j]}/RQH_full_${NUC[f]}_${TIMES[j]}_${i}.secondpart.tpr"
			fi
			echo -e r757-1174 '\n' q | $GMX51 make_ndx -f $TPR -o Y.ndx >& out1
			sed 's/r_757-1174/SecY/g' Y.ndx -i
			for k in {1260..1277}; do
				echo -e "r${k}"'\n'q | $GMX51 make_ndx -f $TPR -o ${k}.ndx -n Y.ndx >& out2
				echo -e "r_${k}"'\n'SecY | $GMX51 mindist -f $XTC -s $TPR -n ${k}.ndx -od ${NUC[f]}_${TIMES[j]}_${i}_${k}.xvg -nice -19 -tu us >& out3
				echo ${k} >> timespent_${NUC[f]}_${TIMES[j]}_${i}.xvg
				num=`egrep -v "\#|\@" ${NUC[f]}_${TIMES[j]}_${i}_${k}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^0|^1|^2|^3" | wc -l`
				tot=`egrep -v "\#|\@" ${NUC[f]}_${TIMES[j]}_${i}_${k}.xvg | wc -l`
				occ=`echo "scale=3 ; $num / $tot * 100" | bc`
				echo $occ >> timespent_${NUC[f]}_${TIMES[j]}_${i}.xvg
			done
		done
	done
done
