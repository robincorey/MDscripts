#!/bin/bash

CD=`pwd`
LIPIDS=(DHPC DLPC DPPC)

CheckBound () {
echo "resname $2 and (name PO4) and within 1 of resid 313 and within 1 of resid 441" > site_1.dat
echo "resname $2 and (name PO4) and within 1 of resid 631 and within 1 of resid 123" > site_2.dat
echo "resname $2 and (name PO4) and within 1 of resid 488 and within 1 of resid 291" > site_3.dat
echo "resname $2 and (name PO4) and within 1 of resid 170 and within 1 of resid 609" > site_4.dat
gmx_sse "select" -f md_$1.xtc -s md_$1.tpr -sf site_1.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site1.xvg
gmx_sse "select" -f md_$1.xtc -s md_$1.tpr -sf site_2.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site2.xvg
gmx_sse "select" -f md_$1.xtc -s md_$1.tpr -sf site_3.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site3.xvg
gmx_sse "select" -f md_$1.xtc -s md_$1.tpr -sf site_4.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site4.xvg
}

CheckBound2 () {
echo -e 'r313 | r430' '\n' 'r631 | r112' '\n' 'r488 | r291' '\n' 'r170 | r609' '\n' 'aPO4 & "'${2}'"' '\n'  q | gmx_sse make_ndx -f md_${1}.gro -o sites.ndx
echo -e 'r_313_r_430' '\n' "PO4_&_$2" | gmx_sse mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site1 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
echo -e 'r_631_r_112' '\n' "PO4_&_$2" | gmx_sse mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site2 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
echo -e 'r_488_r_291' '\n' "PO4_&_$2" | gmx_sse mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site3 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
echo -e 'r_170_r_609' '\n' "PO4_&_$2" | gmx_sse mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site4 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
}

CheckBound3 () {
rm -f occ*gro
echo -e 'r313 | r430' '\n' 'r631 | r112' '\n' 'r488 | r291' '\n' 'r170 | r609' '\n' 'aPO4 | aGL1 | aGL2' '\n' '"DHPC" | "DLPC" | "DPPC"' '\n' '"PO4_GL1_GL2" & "DHPC_DLPC_DPPC"' '\n' q | gmx_sse make_ndx -f ../md_${1}.gro -o sites_occ.ndx
echo -e 'r_313_r_430' '\n' "PO4_GL1_GL2_&_DHPC_DLPC_DPPC" | gmx_sse trjorder -f ../md_${i}.xtc -s ../md_${i}.tpr -n sites_occ.ndx -o occ_1.gro -b 500000 -da 1 -com 
echo -e 'r_631_r_112' '\n' "PO4_GL1_GL2_&_DHPC_DLPC_DPPC" | gmx_sse trjorder -f ../md_${i}.xtc -s ../md_${i}.tpr -n sites_occ.ndx -o occ_2.gro -b 500000 -da 1 -com
echo -e 'r_488_r_291' '\n' "PO4_GL1_GL2_&_DHPC_DLPC_DPPC" | gmx_sse trjorder -f ../md_${i}.xtc -s ../md_${i}.tpr -n sites_occ.ndx -o occ_3.gro -b 500000 -da 1 -com
echo -e 'r_170_r_609' '\n' "PO4_GL1_GL2_&_DHPC_DLPC_DPPC" | gmx_sse trjorder -f ../md_${i}.xtc -s ../md_${i}.tpr -n sites_occ.ndx -o occ_4.gro -b 500000 -da 1 -com
}

AnalyseBound() {
for k in {1..4}
do
	notbound=`grep -v '\#' num_${1}_${2}_site${k}.xvg | grep -v '\@' | grep -c "  0$"`
	total=`grep -v '\#' num_${1}_${2}_site${k}.xvg | grep -v '\@' | grep -c [0-9]`
	bound=`echo "scale=4; $notbound / $total" | bc`
	echo ${1}_${2} $bound >> ../lipid_sites.xvg
done
}

AnalyseBound2() {
for k in {1..4}
do
	grep -A10 "636GLN    SC1" occ_${k}.gro | grep "SC1\|PO4\|GL1\|GL2" > occ_${k}b.gro	
	awk '/SC/{nr[NR+1]}; NR in nr' occ_${k}b.gro > ${k}PO4.gro
	awk '{print $1}' ${k}PO4.gro | sed 's/[^0-9]*//g' | sort -u -n > lipid_array_${k}.txt
	mapfile -t ARRAY < lipid_array_${k}.txt
	for (( i=0; i<${#ARRAY[@]}; i++ ))
	do
        	j=`echo "$i + 1" | bc`
        	lipid=`grep " ${ARRAY[i]}[[:alpha:]]" ${k}PO4.gro | awk '{print $1}' | tr -d 1234567890 | tail -n 1`
		rm -f ${ARRAY[i]}_occ_${k}.xvg
		echo ${ARRAY[i]} > ${ARRAY[i]}_occ_${k}.xvg 
		echo $lipid >> ${ARRAY[i]}_occ_${k}.xvg
		while read line
		do
			if [[ $line == *${ARRAY[i]}* ]]
			then
				echo $j >> ${ARRAY[i]}_occ_${k}.xvg
			else
				echo nan >> ${ARRAY[i]}_occ_${k}.xvg
			fi
		done < ${k}PO4.gro
	done
done		
}

CalculateKinetics () {
for k in {1..4}
do
        rm -f ratios_${k}*.xvg
        mapfile -t ARRAY < lipid_array_${k}.txt
        for (( i=0; i<${#ARRAY[@]}; i++ ))
	do
		j=`echo "$i + 1" | bc`
                lipid=`grep " ${ARRAY[i]}[[:alpha:]]" ${k}PO4.gro | awk '{print $1}' | tr -d 1234567890 | tail -n 1`
                # Remove all nan with j before, then after. Loop 8 times to remove upto 24 nans in a row = 4.8 ns
		cp ${ARRAY[i]}_occ_${k}.xvg ${ARRAY[i]}_occ_${k}_in.xvg
		for l in {0..7}
		do
			perl -0777 -pe "s/${j}.*\n\K.*nan//g" ${ARRAY[i]}_occ_${k}_in.xvg > ${ARRAY[i]}_occ_${k}_del_a.xvg
			perl -0777 -pe "s/\n\K.*nan.*\n(?=.*$j)//g" ${ARRAY[i]}_occ_${k}_del_a.xvg > ${ARRAY[i]}_occ_${k}_del_b.xvg
			sed '/^$/d' ${ARRAY[i]}_occ_${k}_del_b.xvg -i
			perl -0777 -pe "s/${j}.*\n\K.*nan.*\n(?=.*${j})//g" ${ARRAY[i]}_occ_${k}_del_b.xvg > ${ARRAY[i]}_occ_${k}_in.xvg
                done
		mv ${ARRAY[i]}_occ_${k}_in.xvg ${ARRAY[i]}_occ_${k}_kinetics.xvg
		rm -f ${ARRAY[i]}_occ_${k}_kinetics_del*.xvg
                # awk one liner to get unbroken sequence of binding events
                awk -v j=$j 'prev == j {if (!k) k++; if ($1 == j) { k++ } else { a[k]++; k = 0 } } { prev = $1 }END{ for (i in a) print i ":" a[i] }' ${ARRAY[i]}_occ_${k}_kinetics.xvg >> ratios_${k}_$lipid.xvg
	done
done
}

PlotRatios () {
LIPIDS=(DHPC DLPC DPPC)
for k in {1..4}
do
	for i in ${LIPIDS[@]}
	do
		# 100 frames = 20 ns 
		echo @ ${k}_${i} > bound_${k}_${i}.xvg 
		awk -F ':' '$1 >= 100' ratios_${k}_${i}.xvg >> bound_${k}_${i}.xvg	
		cat ../../*/*/bound_${k}_${i}.xvg > ../../bound_${k}_${i}_all.xvg
	done
done
}


for i in {1..4}
do
	cd $CD
	if [ -d Run_${i} ]
	then
		cd Run_${i}
		mkdir -p Binding_analysis
		cd Binding_analysis
		#CheckBound3 $i
		#AnalyseBound2 $i
		#CalculateKinetics $i
		PlotRatios $i
		rm -f *#*
	fi
done

