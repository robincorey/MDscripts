#!/bin/bash

ARRAY=(ATP ADP)
CD=`pwd`
SETUP=/home/birac/Desktop/Preprotein_folding/sol

function setup {
pdb2gmx -f ${NAME}.pdb -o ${NAME}.gro -p ${NAME}.top -i ${NAME}_posre.itp -water spc -ff oplsaa >& pdb2gmxout
genbox -cp ${NAME}.gro -o ${NAME}_sol.gro -p ${NAME}.top -cs >& solvout 
${GMX51} grompp -f ${SETUP}/forgenion.mdp -c ${NAME}_sol.gro -p ${NAME}.top -o forgenion >& grompp1
echo SOL | ${GMX51} genion -s forgenion -p ${NAME}.top -o ions.gro  -pname NA -nname CL -neutral >& genion1

if [[ ! -f ions.gro ]]; then
        echo genion failed
        exit 0
fi

echo -e "Waters | IONS" '\n' q | make_ndx -f ions.gro -o groups.ndx >& mnd1

newna=`grep -c NA ions.gro`
newcl=`grep -c CL ions.gro`
newsol=`grep -c OW ions.gro`

#sed "s/SOL /SOL ${newsol} ;/g" newtop.top -i
#sed "s/NA /NA ${newna} ;/g" newtop.top -i
#sed "s/CL /CL ${newcl} ;/g" newtop.top -i

${GMX51} grompp -f ${SETUP}/em.mdp -c ions.gro -p ${NAME}.top -o ${NAME}_em -n groups.ndx >& grompp2
${GMX51} mdrun -v -deffnm ${NAME}_em >& md2
${GMX51} grompp -f ${SETUP}/pr_GPU_5.1_fold.mdp -c ${NAME}_em -p ${NAME}.top -o ${NAME}_pr -n groups.ndx >& grompp3
${GMX51} mdrun -v -deffnm ${NAME}_pr >& md3
}

for (( i=0; i<${#ARRAY[@]}; i++ )); do
	cd $CD
	mkdir -p ${ARRAY[i]}
	cd ${ARRAY[i]}	
	if [ ${ARRAY[i]} == 'ATP' ]; then
		TIMES=(70 80 110)
	elif [ ${ARRAY[i]} == 'ADP' ]; then
		TIMES=(50 65 80)
	else
		echo NUC loop failed
		exit 0
	fi
	CD2=`pwd`
	for (( j=0; j<${#TIMES[@]}; j++ )); do
		echo ${ARRAY[i]} ${TIMES[j]} started
		cd ${CD2}
		GRO=/home/birac/Desktop/Preprotein_folding/${ARRAY[i]}/Tandem_${ARRAY[i]}_${TIMES[j]}/Tandem_${ARRAY[i]}_${TIMES[j]}_em.gro
		mkdir -p ${TIMES[j]}
		cd ${TIMES[j]}
		cp ${GRO} ${ARRAY[i]}.${TIMES[j]}.gro
		sed '/1260GLU      N2/,$!d' ${ARRAY[i]}.${TIMES[j]}.gro > ${ARRAY[i]}.${TIMES[j]}.prot.gro
		sed '/1278ASN      N2/,$d' ${ARRAY[i]}.${TIMES[j]}.prot.gro -i
		length=`wc -l ${ARRAY[i]}.${TIMES[j]}.prot.gro | awk '{print $1}'`
		sed "1 i\ ${length}" ${ARRAY[i]}.${TIMES[j]}.prot.gro -i
		sed "1 i\ Protein" ${ARRAY[i]}.${TIMES[j]}.prot.gro -i
		echo "5 5 5" >> ${ARRAY[i]}.${TIMES[j]}.prot.gro
		editconf -f ${ARRAY[i]}.${TIMES[j]}.prot.gro -o ${ARRAY[i]}.${TIMES[j]}.prot.box.gro -box 6 6 6 -bt triclinic >& edc1
		editconf -f ${ARRAY[i]}.${TIMES[j]}.prot.box.gro -o ${ARRAY[i]}.${TIMES[j]}.pdb >& edc2
		NAME=${ARRAY[i]}.${TIMES[j]}	
		setup
	done
done

