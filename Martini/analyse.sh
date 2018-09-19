#!/bin/bash

rm -f *_data
rm -f times.xvg

for g in {200..200..25}
do
echo $g >> times.xvg

for i in {1..5}
do
	#for j in {0..20}
	#do
	#cat dhdl_${j}_$i.xvg
	for j in `seq 1 $i`
	do
		printf "*$j.xvg " >> ${i}_${g}_data
	done
	vals=`cat ${i}_${g}_data`
	#echo $vals
	end=`echo "($g * 1000) + 50000" | bc`
	#echo $end
	gmx_sse bar -f $vals -o ${i}_bar -oi ${i}_barint  -temp 323 -b 50000 -e $end >& ${g}_ns_out_$i &
	#alchemical_analysis  -t 323 -s 50000 -j ${i}_f -q ${i}.xvg -e false >& out &
#done
done

done

