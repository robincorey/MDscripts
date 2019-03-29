#!/bin/bash

XTC=$1 # centred on protein
TPR=$2
XTC2=${XTC%.*}

echo -e aROH '\n' q | gmx_avx make_ndx -f ${TPR} -o ROH.ndx
echo -e ROH '\n' ROH | gmx_avx trjconv -f ${XTC} -n ROH.ndx -s ${TPR} -o ROH.xtc -center
echo ROH | gmx_avx convert-tpr -s $TPR -o ROH.tpr -n ROH.ndx
echo System | gmx_avx traj -f ROH.xtc -s ROH -ox ROH -z yes -x no -y no -tu ns -nice -15

centre=`tail -n 1 ROH.xvg | awk '{sum=0; for(i=2; i<=NF; i++){sum+=$i}; sum/=NF; print sum}'`

wc=`awk '{print NF}' ROH.xvg | tail -n 1`
for (( i=1; i<${wc}; i++ ))
do 
	awk -v var=$i '{print $var}' ROH.xvg | if grep -q "^[4-7]" && grep -q "^[8-9]"; then echo $i ; fi
done
