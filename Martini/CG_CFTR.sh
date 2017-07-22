#!/bin/bash

prot1=5uak_rotate.pdb
CD=`pwd`

rm -f vacuum_min.gro *_box.gro sol.gro ions.gro vacuum_min.tpr genion.tpr min_test.tpr

prot=${prot1//.pdb}

mkdir -p CH_test3
cd CH_test3

cp /home/birac/Desktop/CoarseGrainSims/KIT/martinize.py .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization-vacuum.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_1_CFTR.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_2_CFTR.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_3_CFTR.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/*itp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/waterbox.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/DPPC_single.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/CHOL_single.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/DLiPC_single.gro .

#<<'END'
rm -f ${prot}-CG.pdb

python martinize.py -f ../${prot}.pdb -o ${prot}.top -p ${prot}-posres.itp -x ${prot}-CG.pdb -dssp mkdssp -p backbone -ff elnedyn22 >& one

if [ ! -f ${prot}-CG.pdb ]; then
        echo "
        #####################################
        Error in martinize - PDB not created!
        #####################################"
        exit 0
else
        echo martinize done
fi

sed -i -e "1,2d" $prot.top
sed '1i #include "martini_v2.0_ions.itp"' $prot.top -i
sed '1i #include "martini_v2.0_lipids_CFTR.itp"' $prot.top -i
sed '1i #include "martini_v2.2_aminoacids.itp"' $prot.top -i
sed '1i #include "martini_v2.2.itp"' $prot.top -i

editconf -f ${prot}-CG.pdb -o ${prot}-CG.box.pdb -box 25 25 18   >&  edc1
genbox -cp ${prot}-CG.box.pdb -ci DPPC_single.gro -nmol 840 -try 500 -o DPPC_box.gro -vdwd 0.21 -p ${prot}.top >& box1

if [ ! -f  DPPC_box.gro ]; then 
	echo box1 fails
	exit 0
fi

genbox -cp DPPC_box.gro -ci DLiPC_single.gro -nmol 560 -try 500 -o DPPC_DLiPC_box.gro -vdwd 0.21 -p ${prot}.top >& box2

genbox -cp DPPC_DLiPC_box.gro -ci CHOL_single.gro -nmol 600 -try 500 -o DPPC_DLiPC_CHOL_box.gro -vdwd 0.21 -p ${prot}.top >& box3 # with virtual sites!: M.N. Melo, H.I. Ingolfsson, S.J. Marrink. Parameters for Martini sterols and hopanoids based on a virtual-site description.JCP, 143:243152, 2015. doi:10.1063/1.4937783

if [ ! -f  DPPC_DLiPC_CHOL_box.gro ]; then
        echo box2 fails
        exit 0
fi

echo "box building done!"

PCnum=`grep -c "DPPC   PO4" DPPC_DLiPC_CHOL_box.gro`
DLnum=`grep -c "DIPC   PO4" DPPC_DLiPC_CHOL_box.gro`
CHnum=`grep -c "CHO    ROH" DPPC_DLiPC_CHOL_box.gro`

echo "
DPPC $PCnum" >> ${prot}.top
echo "DIPC $DLnum" >> ${prot}.top
echo "CHOL $CHnum" >> ${prot}.top


genbox -o sol.gro -cs waterbox.gro -cp DPPC_DLiPC_CHOL_box.gro -vdwd 0.21 -p ${prot}.top >& sol #-maxsol 768
sed '/W/d' ${prot}.top -i
SOLnum=`grep -c "     W" sol.gro`
echo "W $SOLnum" >> ${prot}.top

echo q | make_ndx -f sol.gro -o index >& ndx
cat ndx | grep [0-9] | grep atoms | awk '{print $2}' | sort | uniq -d >& ndx2

if [[ -s ndx2 ]] ; then
        while read i; do
                sed '0,/'$i'/! s/'$i'/Dud/' index.ndx -i
        done < ndx2
fi

groups=`grep -c '\[' index.ndx`
memno=$(expr $groups + 0)
solno=$(expr $groups + 1)
 
echo -e  '"'DPPC'"' '|' '"'DIPC'"' '|' '"'CHOL'"' '\n' '"'W'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | make_ndx -f sol.gro -o memsol1 -n index >& ndx3

rm -f genion.tpr
rm -f ions.gro

grompp -f minimization.mdp -p ${prot}.top -c sol.gro -o genion.tpr -n memsol1 >& g2

echo W | genion -s genion.tpr -pname NA+ -nname CL- -neutral -o ions.gro -p ${prot}.top >& ions

if [ ! -f ions.gro ]; then
        echo "
        ##################
        Error in genion!!
        ##################"
        exit 0
else
        echo genion done
fi

sed 's/NA+/NA /g' ions.gro -i
END
echo q | make_ndx -f ions.gro -o index >& ndx
cat ndx | grep [0-9] | grep atoms | awk '{print $2}' | sort | uniq -d > ndx2

if [[ -s ndx2 ]] ; then
        while read i; do
                sed '0,/'$i'/! s/'$i'/Dud/' index.ndx -i
        done < ndx2
fi

groups=`grep -c '\[' index.ndx`
memno=$(expr $groups + 0)
solno=$(expr $groups + 1)
echo -e  '"'DPPC'"'  '|' '"'DIPC'"' '|' '"'CHOL'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q | make_ndx -f ions.gro -o memsol -n index >& ndx4

grompp -f minimization.mdp -p ${prot}.top -c ions.gro -o min_test.tpr -n memsol -maxwarn 2 >& gp3

$GMX51 mdrun -c -deffnm min_test -rdd 3 >& md2

if [[ ! -f min_test.gro ]]; then
	echo min test failed
	exit 0
else
	echo min test worked!
fi

rm -f minimization.tpr equilibration_1.tpr equilibration_2.tpr 

# Do a short minin with posres - need define = -DPOSRES in mdp

grompp -f minimization.mdp -c ions.gro -p ${prot}.top -o minimization.tpr -n memsol.ndx -maxwarn 2  >& gp4

if [ ! -f minimization.tpr ]; then
        echo "
        #################################
        Error in building minimization !!
        #################################"
        exit 0
fi


$GMX51 mdrun -deffnm minimization -v  >& md3

if [ ! -f minimization.gro ]; then
        echo "
        ########################
        Error in minimization !!
        ########################"
        exit 0
else
	echo minimization worked!
fi
# Then crack out an equilibration run
$GMX51 grompp -f equilibration_1_CFTR.mdp -c minimization.gro -p ${prot}.top -o equilibration_1.tpr -n memsol.ndx -maxwarn 2 >& gp5
$GMX51 mdrun -deffnm equilibration_1 -v  >& md4  # Takes about 25 mins on 1 CPU 
$GMX51 grompp -f equilibration_2_CFTR.mdp -c equilibration_1.gro -p ${prot}.top -o equilibration_2.tpr -n memsol.ndx -maxwarn 2 >& gp6
$GMX51 mdrun -deffnm equilibration_2 -v  >& md5  # Takes about 150 minutes on 1 CPU
$GMX51 grompp -c equilibration_2.gro -f equilibration_3_CFTR.mdp -p ${prot}.top -o equilibration_3 -n memsol.ndx  -maxwarn 2 >& gp7

#$GMX51 mdrun -deffnm equilibration_3 -v >& md5
