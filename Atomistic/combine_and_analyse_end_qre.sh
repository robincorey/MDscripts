#!/bin/bash

CD=`pwd`
DESCIN=RQH_full
NUC=(ATP)
ARRAY=(70 80 110)

for (( k=0; k<${#NUC[@]}; k++ )); do
mkdir -p ${NUC[k]}_contactres
cd ${NUC[k]}_contactres
rm DSSP_averages*txt
	for (( i=0; i<${#ARRAY[@]}; i++ )); do
		for j in {0..4}; do
			SYS=${DESCIN}_${ARRAY[i]}
			TPR=../${SYS}/${SYS}_${j}.tpr
			echo -e r1270 '\n' r1273-1275 '\n' q | $GMX51 make_ndx -f ${TPR} -o subs_${j}.ndx >& ndx_${ARRAY[i]}_${j}
			echo  r_1270 | $GMX51 do_dssp -s ${TPR} -f ../${SYS}/${SYS}_${j}.pbc.xtc -sc ${SYS}_${j}.dssp.contact1.xvg -n subs_${j} -nice -19 >& dssp1_${ARRAY[i]}_${j}
			echo  r_1273-1275 | $GMX51 do_dssp -s ${TPR} -f ../${SYS}/${SYS}_${j}.pbc.xtc -sc ${SYS}_${j}.dssp.contact2.xvg -n subs_${j} -nice -19 >& dssp2_${ARRAY[i]}_${j}
			C=(contact1 contact2)
			for (( c=0; c<${#C[@]}; c++ )); do
				ave1=0
				ave2=0
				cat  ${SYS}_${j}.dssp.${C[c]}.xvg | grep -v \# | grep -v \@ >  ${SYS}_${j}.dssp.${C[c]}_data.xvg	
				if ! grep -q 3-Helix ${SYS}_${j}.dssp.${C[c]}_data.xvg; then 
					ave1=`awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' ${SYS}_${j}.dssp.${C[c]}_data.xvg`
				else
					helixcolumn=`grep "3-Helix" ${SYS}_${j}.dssp.${C[c]}_data.xvg | awk '{print $2}' | tr -d "s"`
					helixcolumn2=`echo "$helixcolumn + 4" | bc`
					ave1a=`awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' ${SYS}_${j}.dssp.${C[c]}_data.xvg`
					ave1b=`awk '{ sum += $8; n++ } END { if (n > 0) print sum / n; }' ${SYS}_${j}.dssp.${C[c]}_data.xvg`
					echo "$ave1a + $ave1b"
					ave1=`echo "$ave1a + $ave1b" | bc`
				fi
				echo ${ARRAY[i]}_${j} $ave1 $ave2 >> DSSP_averages_${contact}.txt	
			done
		done
	done
rm -f *#*
done

