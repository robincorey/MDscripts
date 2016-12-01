#!/bin/bash 

ARRAY=(SecY*E.pdb)
ARRAYLEN=${#ARRAY[@]}
ARRAY2=`echo ${ARRAYLEN} - 1 | bc`
CURRENT_DIR=`pwd`
GMX5=/usr/local/gromacs_5/bin

## make sure gro file has edited details

for i in `seq 0 ${ARRAY2}`; do
	cd ${CURRENT_DIR}
	cd ${ARRAY[i]//.pdb}
	cp /home/birac/Desktop/CoarseGrainSims/3DIN/Elastic/Run2/Run3/*mdp .
	echo -e rCARD "|" rDPPC '\n' rW "|" rION '\n' name 17 Mem '\n' name 18 SOL '\n' q | make_ndx -f ${ARRAY[i]//.pdb}.gro -o memsol.ndx
	$GMX5/grompp -f test.mdp -c ${ARRAY[i]//.pdb}.gro -p ${ARRAY[i]//.pdb}.top -o ${ARRAY[i]//.pdb}_cg_test.tpr -n memsol
	$GMX5/mdrun -v -deffnm ${ARRAY[i]//.pdb}_cg_test # test
	if [[ ! -f ${ARRAY[i]//.pdb}_cg_test.gro ]]; then 
		echo ${ARRAY[i]//.pdb} has nae worked
		exit 0
	fi
	$GMX5/grompp -f martini_v2.x_new.mdp -c ${ARRAY[i]//.pdb}.gro -p ${ARRAY[i]//.pdb}.top -o ${ARRAY[i]//.pdb}_cg.tpr -n memsol
	$GMX5/editconf -f ${ARRAY[i]//.pdb}_cg.tpr -o ${ARRAY[i]//.pdb}.tpr.gro
done
