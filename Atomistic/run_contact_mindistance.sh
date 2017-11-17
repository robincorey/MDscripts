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
			rm -f ${NUC[f]}_${TIMES[j]}_${i}.xvg
			if [ ${NUC[f]} == 'ATP' ]; then
				XTC="../../${NUC[f]}/RQH_full_Data_Archer/RQH_full_${TIMES[j]}/RQH_full_${TIMES[j]}_${i}.xtc"
				TPR="../../${NUC[f]}/RQH_full_Data_Archer/RQH_full_${TIMES[j]}/RQH_full_${TIMES[j]}_${i}.tpr"
			 elif [ ${NUC[f]} == 'ADP' ]; then
				XTC="../../${NUC[f]}/RQH_full_Data_BC3_Ian/RQH_full_${NUC[f]}_${TIMES[j]}/RQH_full_${NUC[f]}_${TIMES[j]}_${i}.secondpart.xtc"
				TPR="../../${NUC[f]}/RQH_full_Data_BC3_Ian/RQH_full_${NUC[f]}_${TIMES[j]}/RQH_full_${NUC[f]}_${TIMES[j]}_${i}.secondpart.tpr"
			fi
			echo -e "r1260 & aCA" '\n' "r1268 & aCA" '\n' "r1269 & aCA" '\n' "r1277 & aCA" '\n' q | $GMX51 make_ndx -f $TPR -o dist_${i}.ndx >& out1
			echo -e "r_1260_&_CA"'\n'"r_1268_&_CA" | $GMX51 mindist -f $XTC -s $TPR -n dist_${i} -od ${NUC[f]}_${TIMES[j]}_${i}_${k}_dist_half1.xvg -nice -19 -tu us >& out3
			echo -e "r_1269_&_CA"'\n'"r_1277_&_CA" | $GMX51 mindist -f $XTC -s $TPR -n dist_${i} -od ${NUC[f]}_${TIMES[j]}_${i}_${k}_dist_half2.xvg -nice -19 -tu us >& out4
			ave1=`tail -n 2500 ${NUC[f]}_${TIMES[j]}_${i}_${k}_dist_half1.xvg | awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }'`
			ave2=`tail -n 2500 ${NUC[f]}_${TIMES[j]}_${i}_${k}_dist_half2.xvg | awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }'`
			echo ${k} $ave1 $ave2 >> ${NUC[f]}_${TIMES[j]}_${i}.xvg
		done
	done
done
