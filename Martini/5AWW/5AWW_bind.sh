#!/bin/bash
#not easy

CL=(CL NoCL)
CD=`pwd`

for (( j=0; j<${#CL[@]}; j++ )); do
	cd $CD
	cd insane_${CL[j]}
	for i in {0..15}; do
		tail -n 2500 Sec_dist_${i}.xvg | awk -F , '{print $2}' > Sec_dist_${i}_tail.xvg 
		#sed -n '1~4p' Sec_dist_${i}_tail.xvg nice code, but ignroing
		# 0.8 seems a good cutoff
		sed '/e+00/d' Sec_dist_${i}_tail.xvg | sed '/e+01/d' |  sed '/e+02/d' | sed '/^9/d'  > Sec_dist_${i}_${CL[j]}.xvg
	done
	cat  Sec_dist_*_${CL[j]}.xvg  > ${CL[j]}_dimer.xvg
done
