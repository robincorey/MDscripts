#!/bin/bash
# Script to measure dist between key residues and CL head groups

ARRAY=(*final.xtc)  # processed xtc, in folder
CURRENTDIR=`pwd`
TPR=unbiased.tpr
GRO=minimization.gro
RES=$1

if [ -z "${RES}" ]; then
	echo "add a variable"
	exit 0
fi

##############################
## Need to split by chain
echo -e "aPO1 | aPO2 | aGL0" '\n' q | make_ndx -f $GRO -o index.ndx
## First, define chains
chain1=`grep " 1MET     BB" $GRO | awk '{print $3}'`
chain2=`grep " 1ALA     BB" $GRO | awk '{print $3}'`; chain1b=`echo $chain2 - 1 | bc`
chain3=`grep " 1GLU     BB" $GRO | awk '{print $3}'`; chain2b=`echo $chain3 - 1 | bc`
chain4=`grep " 1HIS     BB" $GRO | awk '{print $3}'`; chain3b=`echo $chain4 - 1 | bc`
chain4b=`grep " 817DPPC   NC3" $GRO | awk '{print $3}'`
chain4c=`echo "$chain4b - 1" | bc`
echo -e a1-$chain0b'\n'a$chain1-$chain1b'\n'a$chain2-$chain2b'\n'a$chain3-$chain3b'\n'a$chain4-$chain4c'\n'q | make_ndx -f $GRO -o chains.ndx -n index.ndx
indices=`grep -c "\[" chains.ndx`
## Groups for each chain
PO4=`echo "$indices - 5" | bc`
A=`echo "$indices - 4" | bc`
B=`echo "$indices - 3" | bc`
C=`echo "$indices - 2" | bc`
D=`echo "$indices - 1" | bc`
echo $A | editconf -f $GRO -o A.pdb -n chains.ndx
echo $B |editconf -f $GRO -o Y.pdb -n chains.ndx
echo $C |editconf -f $GRO -o E.pdb -n chains.ndx
echo $D |editconf -f $GRO -o G.pdb -n chains.ndx

ARRAYA=(`grep "BB  ${RES}" A.pdb | awk '{print $5}'`)  # GLN as reasonably few, similar size
ARRAYY=(`grep "BB  ${RES}" Y.pdb | awk '{print $5}'`)  # to lys, polar, and likely to be
ARRAYE=(`grep "BB  ${RES}" E.pdb | awk '{print $5}'`)  # located at membrane surface (unlike
ARRAYG=(`grep "BB  ${RES}" G.pdb | awk '{print $5}'`)  # hydrophobic res)

ARRAYALEN=`echo ${#ARRAYA[@]} - 1 | bc`
ARRAYYLEN=`echo ${#ARRAYY[@]} - 1 | bc`
ARRAYELEN=`echo ${#ARRAYE[@]} - 1 | bc`
ARRAYGLEN=`echo ${#ARRAYG[@]} - 1 | bc`

# Do distance chain by chain
## SecA
for j in `seq 0 ${ARRAYLEN}`; do

rm timespent_${ARRAY[j]}_${RES}.xvg
echo timespent > timespent_${ARRAY[j]}_${RES}.xvg

## :%s/${ARRAYG\[i\]}\./${ARRAYG\[i\]}_${ARRAY\[j\]}\./g
for i in `seq 0 ${ARRAYALEN}`; do
	rm tempA.ndx
	echo -e "$A & r${ARRAYA[i]}"'\n'q | make_ndx -f $GRO -o tempA.ndx -n chains.ndx
	if [[ ! -f tempA.ndx ]]; then
		echo "#################
		make_ndx failed!
		##############"
		exit 0
	fi
	echo -e $indices'\n'PO1_PO2_GL0 | g_mindist -f ${ARRAY[j]}.pbc.xtc -s $TPR -n tempA -od SecA_${ARRAYA[i]}_${ARRAY[j]}_${RES}.xvg -nice -5 -tu us
	echo SecA_${ARRAYA[i]} >> timespent_${ARRAY[j]}_${RES}.xvg 
	num=`egrep -v "\#|\@" SecA_${ARRAYA[i]}_${ARRAY[j]}_${RES}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
	tot=`egrep -v "\#|\@" SecA_${ARRAYA[i]}_${ARRAY[j]}_${RES}.xvg | wc -l`
	occ=`echo "scale=3 ; $num / $tot * 100" | bc`
	echo $occ >> timespent_${ARRAY[j]}_${RES}.xvg
done

for i in `seq 0 ${ARRAYYLEN}`; do
        rm tempY.ndx
	echo -e "$B & r${ARRAYY[i]}"'\n'q | make_ndx -f $GRO -o tempY.ndx -n chains.ndx
        if [[ ! -f tempY.ndx ]]; then
                echo "#################
                make_ndx failed!
                ##############"
                exit 0
        fi
	echo -e $indices'\n'PO1_PO2_GL0 | g_mindist -f ${ARRAY[j]}.pbc.xtc -s $TPR -n tempY -od SecY_${ARRAYY[i]}_${ARRAY[j]}_${RES}.xvg -nice -5 -tu us
        echo SecY_${ARRAYY[i]} >> timespent_${ARRAY[j]}_${RES}.xvg
        num=`egrep -v "\#|\@" SecY_${ARRAYY[i]}_${ARRAY[j]}_${RES}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
        tot=`egrep -v "\#|\@" SecY_${ARRAYY[i]}_${ARRAY[j]}_${RES}.xvg | wc -l`
        occ=`echo "scale=3 ; $num / $tot * 100" | bc`
        echo $occ >> timespent_${ARRAY[j]}_${RES}.xvg
done

for i in `seq 0 ${ARRAYELEN}`; do
        rm tempE.ndx
	echo -e "$C & r${ARRAYE[i]}"'\n'q | make_ndx -f $GRO -o tempE.ndx -n chains.ndx
        if [[ ! -f tempE.ndx ]]; then
                echo "#################
                make_ndx failed!
                ##############"
                exit 0
        fi
	echo -e $indices'\n'PO1_PO2_GL0 | g_mindist -f ${ARRAY[j]}.pbc.xtc -s $TPR -n tempE -od SecE_${ARRAYE[i]}_${ARRAY[j]}_${RES}.xvg -nice -5 -tu us
        echo SecE_${ARRAYE[i]} >> timespent_${ARRAY[j]}_${RES}.xvg
        num=`egrep -v "\#|\@" SecE_${ARRAYE[i]}_${ARRAY[j]}_${RES}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
        tot=`egrep -v "\#|\@" SecE_${ARRAYE[i]}_${ARRAY[j]}_${RES}.xvg | wc -l`
        occ=`echo "scale=3 ; $num / $tot * 100" | bc`
	echo $occ >> timespent_${ARRAY[j]}_${RES}.xvg
done

for i in `seq 0 ${ARRAYGLEN}`; do
        rm tempG.ndx
	echo -e "$D & r${ARRAYG[i]}"'\n'q | make_ndx -f $GRO -o tempG.ndx -n chains.ndx
        if [[ ! -f tempG.ndx ]]; then
                echo "#################
                make_ndx failed!
                ##############"
                exit 0
        fi
	echo -e $indices'\n'PO1_PO2 | g_mindist -f ${ARRAY[j]}.pbc.xtc -s $TPR -n tempG -od SecG_${ARRAYG[i]}_${ARRAY[j]}_${RES}.xvg -nice -5 -tu us
        echo SecG_${ARRAYG[i]} >> timespent_${ARRAY[j]}_${RES}.xvg
        num=`egrep -v "\#|\@" SecG_${ARRAYG[i]}_${ARRAY[j]}_${RES}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
        tot=`egrep -v "\#|\@" SecG_${ARRAYG[i]}_${ARRAY[j]}_${RES}.xvg | wc -l`
        occ=`echo "scale=3 ; $num / $tot * 100" | bc`
	echo $occ >> timespent_${ARRAY[j]}_${RES}.xvg

grep Sec timespent_${ARRAY[j]}_${RES}.xvg > Sec_${RES}.xvg

done

done

if [[ ! -d Backups/ ]]; then
	mkdir Backups/
fi

mv *#* Backups/

paste -d ',' timespent_*_${RES}.xvg > timespentall_${RES}.xvg
egrep -v Sec timespentall_${RES}.xvg > values_${RES}.xvg
paste -d ',' Sec_${RES}.xvg values_${RES}.xvg > timespentall_${RES}.xvg
