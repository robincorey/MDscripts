#!/bin/bash

CG=/sansom/s137/bioc1535/Desktop/CG_KIT
LAMBDA=40

rm -f *_data
rm -f times.xvg

check_
prepare_files () {
for i in `seq 0 $LAMBDA`
do 
	for j in {1..5}
	do
		if [[ -f "MD_${i}_${j}.part0001.xvg" ]]
		then 
			cat MD_${i}_${j}.part*xvg > MD_${i}_${j}.all.xvg
		fi
	done
done
}

full_analysis () {
for g in {25..200..25}
do
echo $g >> times.xvg
for i in {1..5}
do
	for j in `seq 1 $i`
	do
		printf "*$j.xvg " >> ${i}_${g}_data
	done
	vals=`cat ${i}_${g}_data`
	end=`echo "($g * 1000) + 50000" | bc`
	gmx_sse bar -f $vals -o ${i}_bar -oi ${i}_barint -temp 323 -b 50000 -e $end >& ${g}_ns_out_$i &
done
done
}

short_analysis () {

rm -f *plot *err

for i in {1..5}
do
	for j in {25..200..25}
	do 
		grep total ${j}_*out*_$i | awk '{print $6}' >> ${i}.plot
		grep total ${j}_*out*_$i | awk '{print $8}' >> ${i}.err	
	done
done
}

plot_fep () {
paste -d ' ' *plot > plot.xvg
paste -d ' ' *err > err.xvg
python $CG/PLOT_FEP_data.py
}

prepare_files
