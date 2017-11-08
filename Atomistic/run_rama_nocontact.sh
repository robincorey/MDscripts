#!/bin/bash

NUC=(ADP)
CD=`pwd`
for (( f=0; f<${#NUC[@]}; f++ )); do 
	cd $CD
	if [ ${NUC[f]} == 'ATP' ]; then
		TIMES=(110) ##70 80 110)
	elif [ ${NUC[f]} == 'ADP' ]; then
		TIMES=(50 65 80)
	else
		echo NUC loop failed
		exit 0
	fi
	for (( d=0; d<${#TIMES[@]}; d++ )); do
		for i in {0..4}; do
			mkdir -p "${NUC[f]}_${TIMES[d]}_${i}_rama"
			WRITE="${NUC[f]}_${TIMES[d]}_${i}_rama"
			XTC="../${NUC[f]}/RQH_full_Data_BC3_Ian/RQH_full_${NUC[f]}_${TIMES[d]}/RQH_full_${NUC[f]}_${TIMES[d]}_${i}.secondpart.xtc"
			TPR="../${NUC[f]}/RQH_full_Data_BC3_Ian/RQH_full_${NUC[f]}_${TIMES[d]}/RQH_full_${NUC[f]}_${TIMES[d]}_${i}.secondpart.tpr"
			echo -e r757-1174 '\n' q | $GMX51 make_ndx -f $TPR -o ${WRITE}/Y.ndx >& ${WRITE}/out1
			echo -e r1260-1268 '\n' r1269-1277 '\n' q  | $GMX51 make_ndx -f $TPR -o ${WRITE}/Y_subs.ndx -n ${WRITE}/Y.ndx >& ${WRITE}/out2
			sed 's/r_757-1174/SecY/g' ${WRITE}/Y_subs.ndx -i
			sed 's/r_1260-1268/half1/g' ${WRITE}/Y_subs.ndx -i
			sed 's/r_1269-1277/half2/g' ${WRITE}/Y_subs.ndx -i
			HALF=(half1 half2)
			for (( k=0; k<${#HALF[@]}; k++ )); do
				echo -e ${HALF[k]} '\n' q  | $GMX51 trjconv -f ${XTC} -s ${TPR} -n ${WRITE}/Y_subs.ndx -o ${WRITE}/${HALF[k]}.xtc >& ${WRITE}/out_trj_${HALF[k]}
				echo -e ${HALF[k]} '\n' q  | $GMX51 convert-tpr -s ${TPR} -o ${WRITE}/${HALF[k]}.tpr -n ${WRITE}/Y_subs.ndx >& ${WRITE}/out_ndx_${HALF[k]}
				$GMX51 rama -s ${WRITE}/${HALF[k]}.tpr -f ${WRITE}/${HALF[k]}.xtc -o ${WRITE}/${HALF[k]}.xvg >& ${WRITE}/out_rama_${HALF[k]} &
			done
		done
	done
done
