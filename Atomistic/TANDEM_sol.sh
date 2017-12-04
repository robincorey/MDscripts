#!/bin/bash

ARRAY=(ATP ADP)
CD=`pwd`

for (( i=0; i<${#ARRAY[@]}; i++ )); do
        cd $CD
        mkdir -p Tandem_${ARRAY[i]}
        cd Tandem_${ARRAY[i]}
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
                mkdir -p ${TIMES[j]}
                cd ${TIMES[j]}
                NAME=${ARRAY[i]}.${TIMES[j]}
 		GRO=../../../${ARRAY[i]}/Tandem_${ARRAY[i]}_${TIMES[j]}/Tandem_${ARRAY[i]}_${TIMES[j]}_em.gro
		grep "^ 1260\|^ 1261\|^ 1262\|^ 1263\|^ 1264\|^ 1265\|^ 1266\|^ 1267\|^ 1268\|^ 1269\|^ 1270\|^ 1271\|^ 1272\|^ 1273\|^ 1274\|^ 1275\|^ 1276\|^ 1277" ${GRO} > Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro
		sed '/POP/d' -i  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro 
		sed '/SOL/d' -i  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro
		sed '/NA/d' -i  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro
		sed '/CL/d' -i  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro
		sed '/MG/d' -i  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro
		sed '/ATP/d' -i  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro
		tail -n 1 ${GRO} >>  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro
		atomnumbercount1=`grep -c [0-9] Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro | awk '{print $1 }'`
		atomnumbercount2=`expr $atomnumbercount1 - 1`
		sed '1i '$atomnumbercount2''  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro -i
		sed '1i Membrane protein' Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro -i
		$GMX51 editconf -f Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.gro -o Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em.pdb >& edc
		sh /home/birac/Desktop/MD_KIT/fold_prot.sh  Tandem_${ARRAY[i]}_${TIMES[j]}_sol_em >& setup 
	done
done

