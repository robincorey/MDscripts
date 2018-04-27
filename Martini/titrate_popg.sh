#!/bin/bash

CD=`pwd`

titrate () {
PG_num=$1
PE_num=`echo "100 - $PG_num" | bc` 
~stansfeld/n06/CGDB/scripts/run_exchange_lipids2.2.sh -np 1 -top $2 -f $3 -ut POPE:POPG -uc $PE_num:$PG_num -lt POPE:POPG -lc $PE_num:$PG_num
}

runsims () {
gmx_avx grompp -f ../md_1us_340.mdp -c exchange_lipid_system.gro -p system.top -o md_1us_$1 -n sys.ndx
gmx_avx mdrun -v -deffnm md_1us_$1 
}

analyse () {
echo -e r244 '\n' r150 '\n' q | gmx_avx make_ndx -f md_1us_$1.tpr -o $1.ndx
echo -e r_244 '\n' POPG | gmx_avx mindist -f md_1us_$1.xtc -s md_1us_$1.tpr -n $1.ndx -od mindist_${1}_site1.xvg
echo -e r_150 '\n' POPG | gmx_avx mindist -f md_1us_$1.xtc -s md_1us_$1.tpr -n $1.ndx -od mindist_${1}_site2.xvg
sed '/#/d' mindist_${1}_site1.xvg -i
sed '/@/d' mindist_${1}_site1.xvg -i
sed '/#/d' mindist_${1}_site2.xvg -i
sed '/@/d' mindist_${1}_site2.xvg -i
bind1=`awk '($2<0.8){ ++count } END{ print count }' mindist_${1}_site1.xvg`
bind2=`awk '($2<0.8){ ++count } END{ print count }' mindist_${1}_site2.xvg`
tot=`awk '($2>0){ ++count } END{ print count }' mindist_${1}_site1.xvg`
occ1=`echo "scale=3;($bind1 / $tot) * 100" | bc`
occ2=`echo "scale=3;($bind2 / $tot) * 100" | bc`
echo $1 $occ1 > occ_1.xvg
echo $1 $occ2 > occ_2.xvg
}

final_analyse () {
bind1=`awk '($2<0.5){ ++count } END{ print count }' mindist_${1}_site1.xvg`
bind2=`awk '($2<0.5){ ++count } END{ print count }' mindist_${1}_site2.xvg`
tot=`awk '($2>0){ ++count } END{ print count }' mindist_${1}_site1.xvg`
occ1=`echo "scale=3;($bind1 / $tot) * 100" | bc`
occ2=`echo "scale=3;($bind2 / $tot) * 100" | bc`
echo $1 $occ1 > occ_1.xvg
echo $1 $occ2 > occ_2.xvg
}

analyse_control () {
echo -e r96 '\n' q | gmx_avx make_ndx -f md_1us_$1.tpr -o ${1}_negcontrol.ndx
echo -e r_96 '\n' POPG | gmx_avx mindist -f md_1us_$1.xtc -s md_1us_$1.tpr -n ${1}_negcontrol.ndx -od mindist_${1}_negcontrol1.xvg
sed '/#/d' mindist_${1}_negcontrol1.xvg -i
sed '/@/d' mindist_${1}_negcontrol1.xvg -i
bind1=`awk '($2<0.5){ ++count } END{ print count }' mindist_${1}_negcontrol1.xvg`
tot=`awk '($2>0){ ++count } END{ print count }' mindist_${1}_negcontrol1.xvg`
occ1=`echo "scale=3;($bind1 / $tot) * 100" | bc`
echo $1 $occ1 > negcontrol_1.xvg
}

OCC=(5 10 20 30 50)

for i in ${OCC[@]}
do
	itp="/sansom/s137/bioc1535/Desktop/Projects/Technical_development/Kinetics/SecY_CL_occupancy/5AWW_POPC/protein-cg.itp"
	gro="/sansom/s137/bioc1535/Desktop/Projects/Technical_development/Kinetics/SecY_CL_occupancy/5AWW_POPC/MD/md.gro"
	cd $CD
	mkdir -p popg_$i
	cd popg_$i
	#titrate $i $itp $gro
	#runsims $i
	analyse $i
	final_analyse $i
	analyse_control $i
done

