# makin' mutated itps

ARRAY=(SecY*E.pdb)
ARRAYLEN=${#ARRAY[@]}
ARRAY2=`echo ${ARRAYLEN} - 1 | bc`
CURRENT_DIR=`pwd`
GMX5=/usr/local/gromacs_5/bin

for i in `seq 0 ${ARRAY2}`; do
	cd ${CURRENT_DIR}
	cd ${ARRAY[i]//.pdb}
	$GMX5/grompp -f martini_v2.x_new.mdp -c ${ARRAY[i]//.pdb}.gro -p ${ARRAY[i]//.pdb}.top -o ${ARRAY[i]//.pdb}_cg.tpr -n memsol
	$GMX5/mdrun -v -deffnm ${ARRAY[i]//.pdb}_cg
done
