#!/bin/bash

ARRAY=(192.small 192_reseed.1.small 192_reseed.2.small 192_reseed.3.small 192_reseed.4.small unbiased_BC3)
ARRAYLEN=`echo ${#ARRAY[@]} - 1 | bc`
CURRENTDIR=`pwd`
TPR=unbiased.tpr
GRO=minimization.gro

##############################

mkdir -p Lys_kinetics
cd Lys_kinetics

SECY=(88 154)
for i in `seq 0 ${ARRAYLEN}`; do
        echo ${ARRAY[i]} started
	for (( j=0; j<${#SECY[@]}; j++ )); do
		rm -f kinetics_${SECY[j]}_${ARRAY[i]}.xvg	
		grep -v \# ../SecY_${SECY[j]}_${ARRAY[i]}.xvg | grep -v \@ | awk '{print $1,$2}' > SecY_${SECY[j]}_${ARRAY[i]}.data.xvg
		sed "s/e-/*10^-/g" SecY_${SECY[j]}_${ARRAY[i]}.data.xvg -i
		sed "s/e+/*10^/g" SecY_${SECY[j]}_${ARRAY[i]}.data.xvg -i
		while read line; do
			temps=`echo $line | awk '{print $1}'`
			dist=`echo $line | awk '{print $2}'`
			dist2=`echo "scale = 4; $dist * 1000" | bc `
			dist3=`echo $dist2 | cut -f1 -d"."`
			if [[ ${dist3} -lt "1000" ]]; then
				echo $temps $dist3 >>  kinetics_${SECY[j]}_${ARRAY[i]}.xvg
			fi
		done < SecY_${SECY[j]}_${ARRAY[i]}.data.xvg 		
	done
done
