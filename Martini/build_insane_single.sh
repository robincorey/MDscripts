#!/bin/bash

# 128 lipid molecules in total.
# build water box

prot1=${1}
lipid1=DSPE
lipid2=CARD 
lipidratio1=9
lipidratio2=1
x=15
y=15
z=19
dm=-8
run=yes

## Get files

cp /home/birac/Desktop/CoarseGrainSims/KIT/martinize.py .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization-vacuum.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_1_mem.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_2_mem.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_3_mem.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/*itp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/water.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/insane_rc.py .


## Define some functions

function UpdateToplogy {
sed "/${lipid1}/d" ${prot}.top -i
sed "/${lipid2}/d" ${prot}.top -i
lipidnum1=`tail -n +2 ${1} | grep ${lipid1} | awk '{print $1}' | sort -u | wc -l`
lipidnum2=`tail -n +2 ${1} | grep ${lipid2} | awk '{print $1}' | sort -u | wc -l`
echo "${lipid1} $lipidnum1" >> ${prot}.top
echo "${lipid2} $lipidnum2" >> ${prot}.top
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

function MakeNdx {
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
echo -e  '"'${lipid1}'"' '|' '"'${lipid2}'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | make_ndx -f ${1} -o memsol.ndx -n index >& ndx3
}

function ChangeIons {
sed 's/NA+ /ION /g' ${1} -i
sed 's/ NA+/  NA/g' ${1} -i
sed 's/CL- /ION /g' ${1} -i
}

function ReportError {
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

function RebuildGro {
sed '/^$/d' ${1} -i
head -n 2 ${1} > rebuild_${1}
tail -n +2 ${1} | grep "BB\|SC" >> rebuild_${1}
tail -n +2 ${1} | grep "${lipid1}" >> rebuild_${1}
tail -n +2 ${1} | grep "${lipid2}" >> rebuild_${1}
tail -n +2 ${1} | grep "W" >> rebuild_${1}
tail -n +2 ${1} | grep "NA" >> rebuild_${1}
tail -n +2 ${1} | grep "CL-" >> rebuild_${1}
tail -n 1 ${1} >> rebuild_${1}
}

## Start process

prot=${prot1//.pdb}

rm -f vacuum_min.gro sol.gro ions.gro vacuum_min.tpr genion.tpr min_test.tpr index.ndx ${prot}-CG.pdb ${prot}.top ${prot}-posres.itp minimization.tpr equilibration_1.tpr equilibration_2.tpr insanetop.txt rebuild_*.gro

ReportError "${prot}.pdb"

python martinize.py -f ${prot}.pdb -o ${prot}.top -p ${prot}-posres.itp -x ${prot}-CG.pdb -dssp mkdssp -p backbone -ff elnedyn22 >& one

ReportError ${prot}-CG.pdb

sed -i -e "1,2d" ${prot}.top
sed '1i #include "martini_v2.0_ions.itp"' ${prot}.top -i
sed '1i #include "martini_v2.0_lipids_cardiolipin.itp"' ${prot}.top -i
sed '1i #include "martini_v2.2_aminoacids.itp"' ${prot}.top -i
sed '1i #include "martini_v2.2.itp"' ${prot}.top -i

python insane_rc.py -f ${prot}-CG.pdb -o ${prot}-mem.gro -x ${x} -y ${y} -z ${z} -l ${lipid1}:$lipidratio1 -l ${lipid2}:$lipidratio1 -sol W -salt 0.15 -op 500 -dm ${dm} > insaneout

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


