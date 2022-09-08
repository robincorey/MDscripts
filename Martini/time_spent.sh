#!/bin/bash
# Script to measure dist between key residues and CL head groups
# Can be used for every sim

# this array needs to match the name of the inout files
ARRAY=(192.small 192_reseed.1.small 192_reseed.2.small 192_reseed.3.small 192_reseed.4.small unbiased_BC3)
ARRAYLEN=`echo ${#ARRAY[@]} - 1 | bc`
CURRENTDIR=`pwd`

# these files need to exist
TPR=unbiased.tpr
GRO=minimization.gro
RES=$1

if [ -z "${RES}" ]; then
	echo "add a fcuking variable you tard"
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
indices=`grep -c "\[" chains.ndx`
## Groups for each chain
PO4=`echo "$indices - 5" | bc`
A=`echo "$indices - 4" | bc`
B=`echo "$indices - 3" | bc`
C=`echo "$indices - 2" | bc`
D=`echo "$indices - 1" | bc`

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

rm linger_${ARRAY[j]}_${RES}.xvg
echo linger > linger_${ARRAY[j]}_${RES}.xvg

for i in `seq 0 ${ARRAYALEN}`; do
	if [[ ! -f tempA.ndx ]]; then
		echo "#################
		make_ndx failed!
		##############"
		exit 0
	fi
	echo SecA_${ARRAYA[i]} >> linger_${ARRAY[j]}_${RES}.xvg 
	while IFS= read -r line
        do
                dist=`echo $line | awk '{print $2}'`
                if [[ $dist =~ "-01" ]] && [ $dist == 7* ] || [ $dist == 6* ] || [ $dist == 5* ] || [ $dist == 4* ] || [ $dist == 3* ] || [ $dist == 2* ] || [ $dist == 1* ]; then
                        echo $line
                        #rep "$line" $system >> del.gro
                        #sed "/$delme/d" new.gro -i
                fi
        done < SecA_${ARRAYA[i]}_${ARRAY[j]}.xvg
done

exit 0
for i in `seq 0 ${ARRAYYLEN}`; do
        if [[ ! -f tempY.ndx ]]; then
                echo "#################
                make_ndx failed!
                ##############"
                exit 0
        fi
        echo SecY_${ARRAYY[i]} >> linger_${ARRAY[j]}_${RES}.xvg
        num=`egrep -v "\#|\@" SecY_${ARRAYY[i]}_${ARRAY[j]}_${RES}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
        echo $num >> linger_${ARRAY[j]}_${RES}.xvg
done

for i in `seq 0 ${ARRAYELEN}`; do
        if [[ ! -f tempE.ndx ]]; then
                echo "#################
                make_ndx failed!
                ##############"
                exit 0
        fi
        echo SecE_${ARRAYE[i]} >> linger_${ARRAY[j]}_${RES}.xvg
        num=`egrep -v "\#|\@" SecE_${ARRAYE[i]}_${ARRAY[j]}_${RES}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
	echo $num >> linger_${ARRAY[j]}_${RES}.xvg
done

for i in `seq 0 ${ARRAYGLEN}`; do
        if [[ ! -f tempG.ndx ]]; then
                echo "#################
                make_ndx failed!
                ##############"
                exit 0
        fi
        echo SecG_${ARRAYG[i]} >> linger_${ARRAY[j]}_${RES}.xvg
        num=`egrep -v "\#|\@" SecG_${ARRAYG[i]}_${ARRAY[j]}_${RES}.xvg | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
	echo $num >> linger_${ARRAY[j]}_${RES}.xvg

grep Sec linger_${ARRAY[j]}_${RES}.xvg > Sec_${RES}.xvg

done

done

if [[ ! -d Backups/ ]]; then
	mkdir Backups/
fi

mv *#* Backups/

paste -d ',' linger_*_${RES}.xvg > lingerall_${RES}.xvg
egrep -v Sec lingerall_${RES}.xvg > values_${RES}.xvg
paste -d ',' Sec_${RES}.xvg values_${RES}.xvg > lingerall_${RES}.xvg
