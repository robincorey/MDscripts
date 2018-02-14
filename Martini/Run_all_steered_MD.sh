#!/bin/bash

ARRAY=(Site_1_A Site_2_A Site_1_B Site_2_B)
LIPID=(918 1005 907 1175)
RES1=(313 631 488 170)
RES2=(441 123 291 609)

CD=`pwd`

for (( i=0; i<${#ARRAY[@]}; i++ ))
do
	cd $CD
	cd ${ARRAY[i]}
	cp ../${ARRAY[i]}.gro .
	echo -e "aPO4 | aGL1 | aGL2" '\n' "r${LIPID[i]}" '\n' '"PO4_GL1_GL2" & "'r_${LIPID[i]}'"' '\n' "r${RES1[i]} | r${RES2[i]}" '\n' aBB '\n' '"BB" & "'r_${RES1[i]}_r${RES2[i]}'"' '\n' q | gmx_avx make_ndx -f ${ARRAY[i]}.gro -o memprot_pull -n ../memsol.ndx
	sed "s/PO4_GL1_GL2_&_r_${LIPID[i]}/Pull_lipid/g" memprot_pull.ndx -i
	sed "s/BB_&_r_${RES1[i]}_r_${RES2[i]}/Pull_site/g" memprot_pull.ndx -i
	gmx_avx grompp -f ../cg-pull-dir_push.mdp -c ${ARRAY[i]}.gro -p ../VRG4*top -o ${ARRAY[i]}_push -n memprot_pull.ndx -maxwarn 2
	gmx_avx grompp -f ../cg-pull-dir_pull.mdp -c ${ARRAY[i]}.gro -p ../VRG4*top -o ${ARRAY[i]}_pull -n memprot_pull.ndx -maxwarn 2
	gmx_avx mdrun -v -deffnm ${ARRAY[i]}_push
	gmx_avx mdrun -v -deffnm ${ARRAY[i]}_pull
done
