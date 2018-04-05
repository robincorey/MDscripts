#!/bin/bash

XTC=$1 # centred on protein
TPR=$2
XTC2=${XTC%.*}

echo -e aPO4 '\n' q | gmx_avx make_ndx -f ${TPR} -o PO4.ndx
#echo -e PO4 '\n' PO4 | gmx_avx trjconv -f ${XTC} -n PO4.ndx -s ${TPR} -o PO4.xtc -center
##echo PO4 | gmx_avx convert-tpr -s $TPR -o PO4.tpr -n PO4.ndx
echo PO4 | gmx_avx traj -f ${XTC} -s ${TPR} -ox PO4 -z yes -x no -y no -tu ns -nice -15 -n PO4

centre=`tail -n 1 PO4.xvg | awk '{sum=0; for(i=2; i<=NF; i++){sum+=$i}; sum/=NF; print sum}'`

wc=`awk '{print NF}' PO4.xvg | tail -n 1`
for (( i=1; i<${wc}; i++ ))
do 
	awk -v var=$i '{print $var}' PO4.xvg | if grep -q "^[4-7]" && grep -q "^[8-9]"; then echo $i ; fi
done
