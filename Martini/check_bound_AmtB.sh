#!/bin/bash

CD=`pwd`
LIPIDS=(POPE POPG CARD)

CheckBound3 () {
rm -f occ*gro
echo -e 'r156 | r341' '\n' 'aPO4 | aGL1 | aGL2' '\n' '"POPE" | "POPG" | "CARD"' '\n' '"PO4_GL1_GL2" & "POPE_POPG_CARD"' '\n' q | gmx_avx make_ndx -f ../md.gro -o sites_occ.ndx
echo -e 'r_156_r_341' '\n' "PO4_GL1_GL2_&_POPE_POPG_CARD" | gmx_avx trjorder -f ../md_6.pbc.fit.xtc -s ../md.tpr -n sites_occ.ndx -o occ.gro -b 500000 -da 1 -com 
}

AnalyseBound2() {
grep -A10 "369VAL    SC1" occ.gro | grep "SC1\|PO4\|GL1\|GL2" > occb.gro	
awk '/SC/{nr[NR+1]}; NR in nr' occb.gro > PO4.gro
awk '{print $1}' PO4.gro | sed 's/[^0-9]*//g' | sort -u -n > lipid_array.txt
mapfile -t ARRAY < lipid_array.txt
for (( i=0; i<${#ARRAY[@]}; i++ ))
do
	j=`echo "$i + 1" | bc`
	lipid=`grep " ${ARRAY[i]}[[:alpha:]]" PO4.gro | awk '{print $1}' | tr -d 1234567890 | tail -n 1`
	rm -f ${ARRAY[i]}_occ.xvg
	echo ${ARRAY[i]} > ${ARRAY[i]}_occ.xvg 
	echo $lipid >> ${ARRAY[i]}_occ.xvg
	while read line
	do
		if [[ $line == *${ARRAY[i]}* ]]
		then
		echo $j >> ${ARRAY[i]}_occ.xvg
		else
			echo nan >> ${ARRAY[i]}_occ.xvg
		fi
	done < PO4.gro
done
}

CalculateKinetics () {
rm -f ratios*.xvg
mapfile -t ARRAY < lipid_array.txt
for (( i=0; i<${#ARRAY[@]}; i++ ))
do
	j=`echo "$i + 1" | bc`
	lipid=`grep " ${ARRAY[i]}[[:alpha:]]" PO4.gro | awk '{print $1}' | tr -d 1234567890 | tail -n 1`
	# Remove all nan with j before, then after. Loop 8 times to remove upto 24 nans in a row = 4.8 ns
	cp ${ARRAY[i]}_occ.xvg ${ARRAY[i]}_occ_in.xvg
	for l in {0..7}
	do
		perl -0777 -pe "s/${j}.*\n\K.*nan//g" ${ARRAY[i]}_occ_in.xvg > ${ARRAY[i]}_occ_del_a.xvg
		perl -0777 -pe "s/\n\K.*nan.*\n(?=.*$j)//g" ${ARRAY[i]}_occ_del_a.xvg > ${ARRAY[i]}_occ_del_b.xvg
		sed '/^$/d' ${ARRAY[i]}_occ_del_b.xvg -i
		perl -0777 -pe "s/${j}.*\n\K.*nan.*\n(?=.*${j})//g" ${ARRAY[i]}_occ_del_b.xvg > ${ARRAY[i]}_occ_in.xvg
	done
	mv ${ARRAY[i]}_occ_in.xvg ${ARRAY[i]}_occ_kinetics.xvg
	rm -f ${ARRAY[i]}_occ_kinetics_del*.xvg
	# awk one liner to get unbroken sequence of binding events
	awk -v j=$j 'prev == j {if (!k) k++; if ($1 == j) { k++ } else { a[k]++; k = 0 } } { prev = $1 }END{ for (i in a) print i ":" a[i] }' ${ARRAY[i]}_occ_kinetics.xvg >> ratios_$lipid.xvg
done
}

PlotRatios () {
LIPIDS=(CARD POPG POPE)
for i in ${LIPIDS[@]}
do
	# 20 frames = 5 ns 
	echo @ ${i} > bound_${i}.xvg 
	awk -F ':' '$1 >= 10' ratios_${i}.xvg >> bound_${i}.xvg	
	cat ../../*/*/bound_${i}.xvg > ../../bound_${i}_all.xvg
done
}


cd $CD
mkdir -p Binding_analysis
cd Binding_analysis
#CheckBound3
#AnalyseBound2 
#CalculateKinetics 
PlotRatios 
#rm -f *#*

