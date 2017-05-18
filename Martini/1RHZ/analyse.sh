#!/bin/bash
CURRENT_DIR=`pwd`
GMX5=/usr/local/gromacs_5/bin
ARRAY=(monomer dimer)
ARRAY2=(CL NoCL)

for i in `seq 0 1`; do
	for j in `seq 0 1`; do
		cd $CURRENT_DIR
		cd ${ARRAY[i]}/${ARRAY2[j]}
		        <<'END'
			echo -e aBB'\n'q | $GMX5/make_ndx -f ${ARRAY[i]}_${ARRAY2[j]}.tpr -o BB
        		echo -e BB'\n'BB | $GMX5/g_rms -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -f ${ARRAY[i]}_${ARRAY2[j]}.xtc -o ${ARRAY[i]}_${ARRAY2[j]}.rmsd.xvg -n BB >& ca_rms &
        		echo -e Protein'\n'Protein | $GMX5/g_rms -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -f ${ARRAY[i]}_${ARRAY2[j]}.xtc -o ${ARRAY[i]}_${ARRAY2[j]}.protein.rmsd.xvg >& ca_rms &
			echo Protein | $GMX5/g_gyrate -f ${ARRAY[i]}_${ARRAY2[j]}.xtc  -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -o rgyr_protein.xvg
        		echo BB | $GMX5/g_gyrate -f ${ARRAY[i]}_${ARRAY2[j]}.xtc  -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -o rgyr_BB.xvg -n BB
        		echo -e aPO1 "|" aPO2 '\n' q | $GMX5/make_ndx -f  ${ARRAY[i]}_${ARRAY2[j]}.tpr -o ${ARRAY[i]}_${ARRAY2[j]}.phos.ndx
        		cp ../../selection_*.head.dat .
			$GMX5/g_select -sf selection_0.5.head.dat -f ${ARRAY[i]}_${ARRAY2[j]}.xtc  -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -selrpos atom -os ${ARRAY[i]}_${ARRAY2[j]}.0.5_head.xvg -n  ${ARRAY[i]}_${ARRAY2[j]}.phos.ndx >& gselhead &
        		$GMX5/g_select -sf selection_1.head.dat -f ${ARRAY[i]}_${ARRAY2[j]}.xtc  -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -selrpos atom -os ${ARRAY[i]}_${ARRAY2[j]}.1_head.xvg -n  ${ARRAY[i]}_${ARRAY2[j]}.phos.ndx >& gsel_1_head &
        		$GMX5/g_select -sf selection_2.head.dat -f ${ARRAY[i]}_${ARRAY2[j]}.xtc  -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -selrpos atom -os ${ARRAY[i]}_${ARRAY2[j]}.2_head.xvg -n  ${ARRAY[i]}_${ARRAY2[j]}.phos.ndx >& gsel_2_head &
END
		#echo System | $GMX5/trjconv -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -f ${ARRAY[i]}_${ARRAY2[j]}.xtc -pbc mol -o ${ARRAY[i]}_${ARRAY2[j]}.pbc.xtc
		echo System | $GMX5/trjconv -s ${ARRAY[i]}_${ARRAY2[j]}.tpr -f ${ARRAY[i]}_${ARRAY2[j]}.pbc.xtc -sep -skip 100 -o ${ARRAY[i]}_${ARRAY2[j]}.pdb
	done
done
