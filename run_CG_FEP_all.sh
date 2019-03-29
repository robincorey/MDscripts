#!/bin/bash

# whaddya want?

minimization=yes
equilibration=yes
production=no		# for testing only
repeats=yes		# making tpr repeats

## options

lipidonly=no 		# yes, no
charge=on 		# on, off, none

## This stuff unlikey to change

LAMBDA=40      		# but only runs 21
CURRENT_DIR=`pwd`
SETUP_DIR=/sansom/s137/bioc1535/Desktop/Projects/Technical_development/FEP/CG_FEP/FEP_setup
SYSTEM=$1
STEP=$(echo "scale=4; 1 / ${LAMBDA}" | bc) 

cd ${CURRENT_DIR}
mkdir -p ${SYSTEM}_FEP
cd ${SYSTEM}_FEP
if [[ ! -f ../${SYSTEM}.gro ]]; then echo No gro file; exit 0 ; fi
if [[ ! -f ../${SYSTEM}.top ]]; then echo No top file; exit 0 ; fi
cp ../${SYSTEM}.gro .
cp ../${SYSTEM}.top .

minimization () {
for i in `seq 0 2 ${LAMBDA}`
do
	cd ${CURRENT_DIR}
	cd ${SYSTEM}_FEP	
        mkdir -p EM_${i}
        cd EM_${i}	
	if [ ! -f EM_${i}.gro ]; then
		if [ $lipidonly = yes ]
		then
			sed "s/##INIT##/${i}/g" ${SETUP_DIR}/em_FEP_lipids_charge_${charge}.mdp > steep_FEP_${i}.mdp
		elif [ $lipidonly = no ] 
		then
                        sed "s/##INIT##/${i}/g" ${SETUP_DIR}/em_FEP_charge_${charge}.mdp > steep_FEP_${i}.mdp
		else 
			echo lipid?
			exit 0
		fi
		if [ ! -s steep_FEP_${i}.mdp ]; then echo mdp failed; exit 0; fi
		cp ${CURRENT_DIR}/*itp ../.
		gmx_avx grompp -f steep_FEP_${i}.mdp -c ../${SYSTEM} -p ../${SYSTEM} -o EM_${i} -maxwarn 2 -n ../../sys
		gmx_avx mdrun -v -deffnm EM_${i}  
	fi
done
}

equilibration () {
for i in `seq 0 2 ${LAMBDA}`
do
	cd ${CURRENT_DIR}
	cd ${SYSTEM}_FEP
        if [[ ! ${i} == 0 ]] && [[ ! ${i} == ${LAMBDA} ]]; then
                prestep=$(echo "scale=4; ( ${STEP} * ${i} ) " | bc)
                step=$(echo "0${prestep}") 
        elif [[ ${i} == ${LAMBDA} ]]; then
                step=1
        elif [[ ${i} == 0 ]]; then
                step=0
        fi
        mkdir -p EQ_${i}
        cd EQ_${i} 	
	cp ../*itp .
        if [ ! -f NPT_${i}.gro ]; then
                if [ $lipidonly = yes ]
                then
                        sed "s/##INIT##/${i}/g" ${SETUP_DIR}/eq_FEP_lipids_charge_${charge}.mdp > NPT_FEP_${i}.mdp
                elif [ $lipidonly = no ]
                then
                        sed "s/##INIT##/${i}/g" ${SETUP_DIR}/eq_FEP_charge_${charge}.mdp > NPT_FEP_${i}.mdp
                else
                        echo lipid?
                        exit 0
                fi
                if [ ! -s NPT_FEP_${i}.mdp ]; then echo mdp failed; exit 0; fi
		gmx_avx grompp -f NPT_FEP_${i}.mdp -c ../EM_${i}/EM_${i}.gro -p ../${SYSTEM}.top -o NPT_${i} -maxwarn 2 -n ../../sys
		gmx_avx mdrun -v -deffnm NPT_${i} > NPT_${i}_out
		if [ ! -f NPT_${i}.gro ]; then
			echo equilibration crashing
			exit 0
		fi
	fi
done
}

production () {	
for i in `seq 0 2 ${LAMBDA}`
do
	cd ${CURRENT_DIR}
        cd ${SYSTEM}_FEP
        mkdir -p MD_${i}
        cd MD_${i}
	cp ../*itp .
	if [ ! -f MD_${i}.gro ]; then
                if [ $lipidonly = yes ]
                then
                        sed "s/##INIT##/${i}/g" ${SETUP_DIR}/md_FEP_lipids_charge_${charge}.mdp > MD_FEP_${i}.mdp
                elif [ $lipidonly = no ]
                then
                        sed "s/##INIT##/${i}/g" ${SETUP_DIR}/md_FEP_charge_${charge}.mdp > MD_FEP_${i}.mdp
                else
                        echo lipid?
                        exit 0
                fi
                if [ ! -s MD_FEP_${i}.mdp ]; then echo mdp failed; exit 0; fi		
		gmx_avx grompp -f MD_FEP_${i}.mdp -c ../EQ_${i}/NPT_${i}.gro -p ../${SYSTEM} -o MD_${i} -maxwarn 2 -n ../../sys
		gmx_avx mdrun -v -deffnm MD_${i} -dhdl dhdl_${i}.xvg -nsteps 500000 > MD_${i}_out
	fi
done
}

repeats () {
for i in `seq 0 2 ${LAMBDA}`
do
        cd ${CURRENT_DIR}
        cd ${SYSTEM}_FEP
        mkdir -p MD_${i}_repeats
        cd MD_${i}_repeats
	if [ $lipidonly = yes ]
	then
		sed "s/##INIT##/${i}/g" ${SETUP_DIR}/md_FEP_lipids_charge_${charge}.mdp > MD_FEP_${i}.mdp
	elif [ $lipidonly = no ]
	then
		sed "s/##INIT##/${i}/g" ${SETUP_DIR}/md_FEP_charge_${charge}.mdp > MD_FEP_${i}.mdp
	else
		echo lipid?
		exit 0
        fi
        for j in {1..5}
        do
                gmx_avx grompp -f MD_FEP_${i}.mdp -c ../EQ_${i}/NPT_${i}.gro -p ../${SYSTEM} -o MD_${i}_${j} -maxwarn 2 -n ../../sys
#               gmx_avx mdrun -v -deffnm MD_${i}_$j -dhdl dhdl_${i}_$j.xvg -nsteps 12500000
        done
done
}


mkdir -p Backups
mv *#* Backups/

if [[ ${minimization} = "yes" ]]; then minimization ; fi
if [[ ${equilibration} = "yes" ]]; then equilibration ; fi
if [[ ${production} = "yes" ]]; then production ; fi
if [[ ${repeats} = "yes" ]]; then repeats ; fi

