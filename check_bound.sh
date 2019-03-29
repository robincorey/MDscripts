#!/bin/bash

CD=`pwd`
LIPIDS=(DHPC DLPC DPPC)

CheckBound () {
echo "resname $2 and (name PO4) and within 1 of resid 313 and within 1 of resid 441" > site_1.dat
echo "resname $2 and (name PO4) and within 1 of resid 631 and within 1 of resid 123" > site_2.dat
echo "resname $2 and (name PO4) and within 1 of resid 488 and within 1 of resid 291" > site_3.dat
echo "resname $2 and (name PO4) and within 1 of resid 170 and within 1 of resid 609" > site_4.dat
gmx_avx "select" -f md_$1.xtc -s md_$1.tpr -sf site_1.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site1.xvg
gmx_avx "select" -f md_$1.xtc -s md_$1.tpr -sf site_2.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site2.xvg
gmx_avx "select" -f md_$1.xtc -s md_$1.tpr -sf site_3.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site3.xvg
gmx_avx "select" -f md_$1.xtc -s md_$1.tpr -sf site_4.dat -selrpos whole_res_com -pbc yes -b 1800000 -os ${2}_size_site4.xvg
}

CheckBound2 () {
echo -e 'r313 | r430' '\n' 'r631 | r112' '\n' 'r488 | r291' '\n' 'r170 | r609' '\n' 'aPO4 & "'${2}'"' '\n'  q | gmx_avx make_ndx -f md_${1}.gro -o sites.ndx
echo -e 'r_313_r_430' '\n' "PO4_&_$2" | gmx_avx mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site1 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
echo -e 'r_631_r_112' '\n' "PO4_&_$2" | gmx_avx mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site2 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
echo -e 'r_488_r_291' '\n' "PO4_&_$2" | gmx_avx mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site3 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
echo -e 'r_170_r_609' '\n' "PO4_&_$2" | gmx_avx mindist -f md_${1}.xtc -s md_${1}.tpr -on num_${1}_${2}_site4 -n sites.ndx -pbc yes -b 500000 -d 0.6 -group yes
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

rm -f lipid_sites.xvg

for i in {1..4}
do
	cd $CD
	cd Run_${i}
	for j in ${LIPIDS[@]}
	do
		CheckBound2 $i $j
		AnalyseBound $i $j
	done
done

