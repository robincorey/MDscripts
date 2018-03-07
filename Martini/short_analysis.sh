#!/bin/bash

LIPIDS=(DHPC DLPC DPPC)

GetOcc () {
echo $1 $3
rm -f ${1}_occ.xvg
echo ${1} > ${1}_occ.xvg
echo ${3} >> ${1}_occ.xvg
while read line
do
	if [[ $line == *$1* ]]
	then
		echo $2 >> ${1}_occ.xvg
	else
		echo nan >> ${1}_occ.xvg
	fi
done < PO4.gro
}

for i in ${LIPIDS[@]}
do
	rm -f ${LIPIDS[@]}_occ.xvg
done

awk '/SC/{nr[NR+1]}; NR in nr' occ_2.gro > PO4.gro
awk '{print $1}' PO4.gro | tr -d [[:alpha:]] | sort -u -n > lipid_array.txt

mapfile -t ARRAY < lipid_array.txt

for (( i=0; i<${#ARRAY[@]}; i++ ))
do
	j=`echo "$i + 1" | bc` 
	lipid=`grep " ${ARRAY[i]}[[:alpha:]]" PO4.gro | awk '{print $1}' | tr -d 1234567890 | tail -n 1`
	GetOcc ${ARRAY[i]} $j $lipid
done

#paste -d ',' *_occ.xvg > lipid_occupancies.xvg
