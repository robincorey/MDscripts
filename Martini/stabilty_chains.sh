#!/bin/bash
# Script to measure dist between key residues and CL head groups
# Can be used for every sim

ARRAY=(192.small 192_reseed.1.small 192_reseed.2.small 192_reseed.3.small 192_reseed.4.small unbiased_BC3)
#ARRAY=(unbiased_BC3)
ARRAYLEN=`echo ${#ARRAY[@]} - 1 | bc`
#ARRAYLEN=0
CURRENTDIR=`pwd`
TPR=unbiased.tpr
GRO=minimization.gro
GMX5=/usr/local/gromacs_5/bin

##############################
echo -e "aBB" '\n' q | make_ndx -f $GRO -o BB.ndx
## First, define chains
chain1=`grep " 1MET     BB" $GRO | awk '{print $3}'`
chain2=`grep " 1ALA     BB" $GRO | awk '{print $3}'`; chain1b=`echo $chain2 - 1 | bc`
chain3=`grep " 1GLU     BB" $GRO | awk '{print $3}'`; chain2b=`echo $chain3 - 1 | bc`
chain4=`grep " 1HIS     BB" $GRO | awk '{print $3}'`; chain3b=`echo $chain4 - 1 | bc`
chain4b=`grep " 817DPPC   NC3" $GRO | awk '{print $3}'`
chain4c=`echo "$chain4b - 1" | bc`
echo -e a$chain1-$chain1b "&" '"BB"' '\n'a$chain2-$chain2b "&" '"BB"' '\n'a$chain3-$chain3b "&" '"BB"' '\n'a$chain4-$chain4c "&" '"BB"' '\n'q | make_ndx -f $GRO -o BBchains.ndx -n BB.ndx
indices=`grep -c "\[" BBchains.ndx`
## Groups for each chain
A=`echo "$indices - 4" | bc`
B=`echo "$indices - 3" | bc`
C=`echo "$indices - 2" | bc`
D=`echo "$indices - 1" | bc`
#echo $A | editconf -f $GRO -o A.pdb -n BBchains.ndx
#echo $B |editconf -f $GRO -o Y.pdb -n BBchains.ndx
#echo $C |editconf -f $GRO -o E.pdb -n BBchains.ndx
#echo $D |editconf -f $GRO -o G.pdb -n BBchains.ndx


for i in `seq 0 ${ARRAYLEN}`; do
	echo -e $A'\n'$A | $GMX5/g_rms -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_A.rmsd.xvg -n BBchains >& ca_rms_$A &
        echo -e $B'\n'$B | $GMX5/g_rms -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_Y.rmsd.xvg -n BBchains >& ca_rms_$B &
	echo -e $C'\n'$C | $GMX5/g_rms -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_E.rmsd.xvg -n BBchains >& ca_rms_$C &
	echo -e $D'\n'$D | $GMX5/g_rms -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_G.rmsd.xvg -n BBchains >& ca_rms_$D &
	echo -e $A | $GMX5/g_rmsf -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_A.rmsf.xvg -n BBchains >& ca_rmsf_$A &
	echo -e $B | $GMX5/g_rmsf -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_Y.rmsf.xvg -n BBchains >& ca_rmsf_$B &
	echo -e $C | $GMX5/g_rmsf -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_E.rmsf.xvg -n BBchains >& ca_rmsf_$C &
	echo -e $D | $GMX5/g_rmsf -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_G.rmsf.xvg -n BBchains >& ca_rmsf_$D &
	echo -e BB'\n'BB | $GMX5/g_rms -s ${TPR} -f ${ARRAY[i]}.pbc.xtc -o ${ARRAY[i]}_BB.rmsd.xvg -n BBchains
done
