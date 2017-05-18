#!/bin/bash
CURRENT_DIR=`pwd`
GMX5=/usr/local/gromacs_5/bin
ARRAY=(Monomer Dimer)
ARRAY2=(CL NoCL)

sleep 11000 

for i in `seq 0 1`; do
	for j in `seq 0 1`; do
		cd $CURRENT_DIR
		cd ${ARRAY[i]}/${ARRAY2[j]}
		$GMX5/mdrun -v -deffnm equilibration_2
	done
done
