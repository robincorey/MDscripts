#!/bin/bash

prot1=${1}
DIR=/sansom/s137/bioc1535/Desktop/CG_KIT
#Lipids: from van Meer, Voelker and Feigenson 2008: 20% CHOL, 80% PL (50% PC, 10% PE, 5% PI, 5% PS, 10% SM)
lipids=(CHOL POPC POPE) 
lipidratios=(4 10 2)
x=15
y=15
z=9
run=yes

## Get files

cp ${DIR}/martinize.py .
cp ${DIR}/minimization-vacuum.mdp .
cp ${DIR}/minimization.mdp .
cp ${DIR}/equilibration_1_mem.mdp .
cp ${DIR}/equilibration_2_mem.mdp .
cp ${DIR}/equilibration_3_mem.mdp .
cp ${DIR}/*itp .
cp ${DIR}/water.gro .
cp ${DIR}/insane_rc.py .


## Define some functions

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

MakeNdx () {
rm -f memsol.ndx index.ndx
echo q | make_ndx -f ${1} -o index >& ndx
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
echo -e ${lipid_ndxs::-1} '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | make_ndx -f ${1} -o memsol.ndx -n index >& ndx3
echo -e ${lipid_ndxs::-1} '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q
}

ChangeIons () {
sed 's/NA+ /ION /g' ${1} -i
sed 's/ NA+/  NA/g' ${1} -i
sed 's/CL- /ION /g' ${1} -i
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

## Start process

prot=${prot1//.pdb}

rm -f vacuum_min.gro sol.gro ions.gro vacuum_min.tpr genion.tpr min_test.tpr index.ndx ${prot}-CG.pdb ${prot}.top ${prot}-posres.itp minimization.tpr equilibration_1.tpr equilibration_2.tpr insanetop.txt rebuild_*.gro

pdb2gmx -f ${prot}.pdb -ignh -q ${prot}_clean.pdb -ff oplsaa -water spc >& pdbout 

ReportError "${prot}_clean.pdb"

memembed -o ${prot}_mem.pdb  ${prot}_clean.pdb >& mem 

python martinize.py -f ${prot}_mem.pdb -o ${prot}.top -p ${prot}-posres.itp -x ${prot}-CG.pdb -dssp /sansom/s137/bioc1535/dssp/mkdssp -p backbone -ff elnedyn22 >& one

ReportError ${prot}-CG.pdb

sed -i -e "1,2d" ${prot}.top
sed '1i #include "martini_v2.0_ions.itp"' ${prot}.top -i
sed '1i #include "martini_v2.0_lipids_rcorey.itp"' ${prot}.top -i
sed '1i #include "martini_v2.2_aminoacids.itp"' ${prot}.top -i
sed '1i #include "martini_v2.2.itp"' ${prot}.top -i

declare_lipids=`for (( i=0; i<${#lipids[@]}; i++ )) ; do echo -n "-l ${lipids[i]}:${lipidratios[i]} " ; done`
echo $declare_lipids
python insane_rc.py -f ${prot}-CG.pdb -o ${prot}-mem.gro -x ${x} -y ${y} -z ${z} ${declare_lipids} -sol W -salt 0.15  > insaneout

echo "
" >> ${prot}.top
sed 's/NA+ /NA /g' insanetop.txt >> ${prot}.top

ChangeIons ${prot}-mem.gro
RebuildGro ${prot}-mem.gro
MakeNdx rebuild_${prot}-mem.gro
UpdateToplogy rebuild_${prot}-mem.gro

grompp -f minimization.mdp -p ${prot}.top -c rebuild_${prot}-mem.gro -o genion.tpr -n memsol >& g2

echo W | genion -s genion.tpr -pname NA+ -nname CL- -neutral -o ions.gro >& ions

ReportError ions.gro
ChangeIons ions.gro
RebuildGro ions.gro
UpdateToplogy rebuild_ions.gro
MakeNdx rebuild_ions.gro

grompp -f minimization.mdp -p ${prot}.top -c rebuild_ions.gro -o min_test.tpr -n memsol >& gp3

${GMX51} mdrun -c -deffnm min_test >& md2

ReportError min_test.gro
RebuildGro min_test.gro

grompp -f minimization.mdp -c rebuild_min_test.gro -p ${prot}.top -o minimization.tpr -n memsol.ndx >& gp4

ReportError minimization.tpr

${GMX51} mdrun -deffnm minimization -v >& md3

ReportError minimization.gro

${GMX51} grompp -f equilibration_1_mem.mdp -c minimization.gro -p ${prot}.top -o equilibration_1.tpr -n memsol.ndx >& gp5
${GMX51} mdrun -deffnm equilibration_1 -v >& md4  

ReportError equilibration_1.gro

${GMX51} grompp -f equilibration_2_mem.mdp -c equilibration_1.gro -p ${prot}.top -o equilibration_2.tpr -n memsol.ndx -maxwarn 2 >& gp6
${GMX51} mdrun -deffnm equilibration_2 -v >& md5  

ReportError equilibration_2.gro

${GMX51} grompp -c equilibration_2.gro -f equilibration_3_mem.mdp -p ${prot}.top -o equilibration_3 -n memsol.ndx -maxwarn 2 >& gp7

if [ $run = yes ]; then
	${GMX51} mdrun -deffnm equilibration_3 -v >& md6
else
	echo now need to run ${GMX51} mdrun -deffnm equilibration_3 -v >& md6
fi



