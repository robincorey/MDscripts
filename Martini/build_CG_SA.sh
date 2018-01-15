#!/bin/bash

prot1=${1}
DIR=/sansom/s137/bioc1535/Desktop/CG_KIT
#Lipids: from van Meer, Voelker and Feigenson 2008: 20% CHOL, 80% PL (50% PC, 10% PE, 5% PI, 5% PS, 10% SM)
lipids=(ERGO DLPC DLPE)
lipidratios=(4 10 2)
run=yes
CD=`pwd`

rm -f vacuum_min.gro PE.gro PECL.gro sol.gro ions.gro vacuum_min.tpr genion.tpr min_test.tpr

prot=${prot1//.pdb}

cp ${DIR}/martinize.py .
cp ${DIR}/minimization.mdp .
cp ${DIR}/minimization-vacuum.mdp .
cp ${DIR}/equilibration_1_mem.mdp .
cp ${DIR}/equilibration_2_mem.mdp .
cp ${DIR}/dynamic_mem.mdp .
for (( i=0; i<${#lipids[@]}; i++ )); do
	cp ${DIR}/GROs/${lipids[i]}_single.gro .
	if [ "$?" != "0" ]; then
    		echo "[Error] copy failed!" 1>&2
    	exit 1
	fi
done
cp ${DIR}/ITPs/martini_v2.0_lipids_all_201506.itp .
cp ${DIR}/ITPs/martini_v2.0_lipids_rcorey.itp .
cp ${DIR}/ITPs/martini_v2.2_aminoacids.itp .
cp ${DIR}/ITPs/martini_v2.2.itp .
cp ${DIR}/ITPs/martini_v2.0_ions.itp .
cp ${DIR}/GROs/water.gro .

RebuildGro () {
sed '/^$/d' ${1} -i
head -n 2 ${1} > rebuild_${1}
sed '1s/^.*$/Membrane protein/g' rebuild_${1} -i
tail -n +2 ${1} | grep "BB\|SC" >> rebuild_${1}
for (( i=0; i<${#lipids[@]}; i++ )); do
        tail -n +2 ${1} | grep "${lipids[i]}" >> rebuild_${1}
done
tail -n +2 ${1} | grep "W" >> rebuild_${1}
tail -n +2 ${1} | grep "NA" >> rebuild_${1}
tail -n +2 ${1} | grep "CL-" >> rebuild_${1}
tail -n 1 ${1} >> rebuild_${1}
}

MakeNdx () {
rm -f memsol.ndx index.ndx
echo q | gmx_sse make_ndx -f ${1} -o index >& ndx
cat ndx | grep [0-9] | grep atoms | awk '{print $2}' | sort | uniq -d > ndx2
if [[ -s ndx2 ]] ; then
        while read i; do
                sed '0,/'$i'/! s/'$i'/Dud/' index.ndx -i
        done < ndx2
fi
groups=`grep -c '\[' index.ndx`
memno=$(expr $groups + 0)
solno=$(expr $groups + 1)
lipid_ndxs=`for (( i=0; i<${#lipids[@]}; i++ )); do echo -n '"'${lipids[i]}'"' '|'; done`
echo -e ${lipid_ndxs::-1} '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | gmx_sse make_ndx -f ${1} -o memsol.ndx -n index >& ndx3
}

ReportError () {
if [ ! -f ${1} ]; then
        echo "
        #################################
        Error: no ${1}
        #################################"
        exit 0
else
        echo "${1} done"
fi
}

UpdateToplogy () {
for (( i=0; i<${#lipids[@]}; i++ )); do
        sed "/${lipids[i]}/d" ${prot}.top -i
        lipidnum=`tail -n +2 ${1} | grep ${lipids[i]} | awk '{print $1}' | sort -u | wc -l`
        echo "${lipids[i]} $lipidnum" >> ${prot}.top
done
if grep -q W ${1} ;
        then sed "/^W/d" ${prot}.top -i
        solnum=`tail -n +2 ${1}| grep -c "W "`
        echo "W $solnum" >> ${prot}.top
fi
if grep -q NA ${1} ;
        then sed "/^NA/d" ${prot}.top -i
        NAnum=`tail -n +2 ${1} | grep -c "NA"`
        echo "NA $NAnum" >> ${prot}.top
fi
if grep -q CL- ${1} ;
        then sed "/^CL-/d" ${prot}.top -i
        CLnum=`tail -n +2 ${1} | grep -c "CL-"`
        echo "CL- $CLnum" >> ${prot}.top
fi
}

ChangeIons () {
sed 's/NA+ /ION /g' ${1} -i
sed 's/ NA+/  NA/g' ${1} -i
sed 's/CL- /ION /g' ${1} -i
sed 's/ CL/CL-/g' ${1} -i
}

rm -f vacuum_min.gro sol.gro ions.gro vacuum_min.tpr genion.tpr min_test.tpr index.ndx ${prot}-CG.pdb ${prot}.top ${prot}-posres.itp minimization.tpr equilibration_1.tpr equilibration_2.tpr insanetop.txt rebuild_*.gro

gmx_sse pdb2gmx -f ${prot}.pdb -ignh -q ${prot}_clean.pdb -ff oplsaa -water spc >& dbout

ReportError "${prot}_clean.pdb"

/sansom/s137/bioc1535/Downloads/memembed-1.15/bin/memembed -o ${prot}_mem.pdb  ${prot}_clean.pdb >& mem

python martinize.py -f ${prot}_mem.pdb -o ${prot}.top -p ${prot}-posres.itp -x ${prot}-CG.pdb -dssp /sansom/s137/bioc1535/dssp/mkdssp -p backbone -ff elnedyn22 >& one

ReportError ${prot}-CG.pdb

sed -i -e "1,2d" ${prot}.top
sed '1i #include "martini_v2.0_ions.itp"' ${prot}.top -i
sed '1i #include "martini_v2.0_lipids_all_201506.itp"' ${prot}.top -i
sed '1i #include "martini_v2.2_aminoacids.itp"' ${prot}.top -i
sed '1i #include "martini_v2.2.itp"' ${prot}.top -i
echo "
" >> ${prot}.top

gmx_sse editconf -f ${prot}-CG.pdb -o ${prot}-CG.box.gro -bt triclinic -d 2 >&  edc

x=`tail -n 1 ${prot}-CG.box.gro | awk {'print $1'}`
y=`tail -n 1 ${prot}-CG.box.gro | awk {'print $2'}`
z=`tail -n 1 ${prot}-CG.box.gro | awk {'print $3'}`
area=`echo "scale=4; $x * $y * 2" | bc`
lipidnumber=`echo "scale=4; $area / 0.7" | bc` # where 0.7 is ~APL

tot=0
for (( i=0; i<${#lipidratios[@]}; i++ )); do
  let totallipid+=${lipidratios[i]}
done
echo totallipid is $totallipid

cp ${prot}-CG.box.gro 0.gro

for (( i=0; i<${#lipids[@]}; i++ )); do
	number=`echo "scale=4; ${lipidratios[i]} / ${totallipid} * $lipidnumber" | bc | cut -f1 -d'.'`
	j=`echo "${i} + 1" | bc`
	gmx_sse insert-molecules -f ${i}.gro -ci ${lipids[i]}_single.gro -nmol $number -try 500 -o ${j}.gro >& box_${i}
	cp ${j}.gro box.gro
done


ReportError box.gro
UpdateToplogy box.gro
RebuildGro box.gro

z2=`echo "$z + 2" | bc`
gmx_sse solvate -cp rebuild_box.gro -cs water.gro -o sol.gro -p -box $x $y $z2 >& sol 

ReportError sol.gro
UpdateToplogy sol.gro
MakeNdx sol.gro
RebuildGro sol.gro

gmx_sse grompp -f minimization.mdp -p ${prot}.top -c rebuild_sol.gro -o genion.tpr -n memsol >& g2
echo W | gmx_sse genion -s genion.tpr -pname NA+ -nname CL- -neutral -o ions.gro >& ions # renumbers

ReportError ions.gro
UpdateToplogy ions.gro
MakeNdx ions.gro
RebuildGro ions.gro

ChangeIons rebuild_ions.gro

gmx_sse grompp -f minimization.mdp -c rebuild_ions.gro -p ${prot}.top -o minimization.tpr -n memsol.ndx >& gp4
gmx_sse mdrun -deffnm minimization -v >& md3

gmx_sse grompp -f equilibration_1_mem.mdp -c minimization.gro -p ${prot}.top -o equilibration_1.tpr -n memsol.ndx >& gp5
gmx_sse mdrun -deffnm equilibration_1 -v >& md4  # Takes about 25 mins on 1 CPU 

gmx_sse grompp -f equilibration_2_mem.mdp -c equilibration_1.gro -p ${prot}.top -o equilibration_2.tpr -n memsol.ndx -maxwarn 2 >& gp6
gmx_sse mdrun -deffnm equilibration_2 -v >& md5  # Takes about 150 minutes on 1 CPU

gmx_sse grompp -c equilibration_2.gro -f dynamic_mem.mdp -p ${prot}.top -o dynamic -n memsol.ndx -maxwarn 2 >& gp7
gmx_sse mdrun -deffnm dynamic -v >& md6

