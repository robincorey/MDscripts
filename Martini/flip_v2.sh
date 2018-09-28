#!/bin/bash

MD=$1
TPR=$2

make_PO4 () { 
echo -e aPO4 '\n' q | gmx_sse make_ndx -f $TPR -o PO4
echo PO4 '\n' System | gmx_sse trjconv -f $MD -s $TPR -center -n PO4 -o ${MD::-4}.center.xtc
echo PO4 | gmx_sse traj -s $TPR -f  ${MD::-4}.center.xtc -ox PO4_pre.xvg -nox -noy -z -tu us -n PO4 -b 0.1
grep -v \# PO4_pre.xvg | grep -v \@ > PO4.xvg
}

make_PO4

centre=`tail -n 1 PO4.xvg | awk '{sum=0; for (i=2; i<=NF; i++) { sum+= $i } print sum/NF}'`
c1=`echo "$centre - 1.5" | bc`
c2=`echo "$centre + 1.5" | bc`
rm -f res.txt
awk -v var="$c1" -v var2="$c2" '{for(i=2;i<=NF;i++){if($i<var){col[i]++};if($i>var2){col1[i]++}}} END{for(j in col){if(col[j]>=1 && col1[j]>=1){print j-1}}}' PO4.xvg > res.txt
if [ ! -s res.txt ] 
then
	echo "##############################
No flipping residues
##############################"
	exit 0
fi
rm -f up.xvg down.xvg
top=`tail -n 1 PO4.xvg | awk -v var="$centre" '{sum=0; for (i=2; i<=NF; i++) {if($i<var) sum+= $i } print sum/(NF/2)}'`
bottom=`tail -n 1 PO4.xvg | awk -v var="$centre" '{sum=0; for (i=2; i<=NF; i++) {if($i>var) sum+= $i } print sum/(NF/2)}'`
starttime=`head -n 1 PO4.xvg | awk '{print $1}'`
endtime=`tail -n 1 PO4.xvg | awk '{print $1}'`
high=`echo "$top - $centre" | bc`
low=`echo "$bottom - $centre" | bc`

while IFS= read -r line
do
	col=`echo "$line + 1" | bc`
	begin=`awk -v var=$col '{print $var}' PO4.xvg | head -n 20 | tail -n 1`
	echo $col $begin
	if [[ "$begin" > "$centre" ]]
	then
		echo down
		echo $line >> down.xvg
	else
		echo up
		echo $line >> up.xvg
	fi
done < res.txt

for j in up down
do
	num=`wc -l ${j}.xvg | awk '{print $1}'`
	sed "s/##CENT##/$centre/g" $CG/PLOT_flip.py > PLOT_flip_${j}.py
	sed "s/##end##/$endtime/g" PLOT_flip_${j}.py -i
	sed "s/##start##/$starttime/g" PLOT_flip_${j}.py -i
	sed "s/##high##/$high/g" PLOT_flip_${j}.py -i
	sed "s/##low##/$low/g" PLOT_flip_${j}.py -i
	sed "s/##num##/$num/g" PLOT_flip_${j}.py -i
	sed "s/##data##/${j}.xvg/g" PLOT_flip_${j}.py -i
	sed "s/##w##/${j}/g" PLOT_flip_${j}.py -i
	python PLOT_flip_${j}.py
done
