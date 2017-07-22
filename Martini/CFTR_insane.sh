#!/bin/bash

prot1=5uak_rotate.pdb
CD=`pwd`

rm -f vacuum_min.gro *_box.gro sol.gro ions.gro vacuum_min.tpr genion.tpr min_test.tpr

prot=${prot1//.pdb}

mkdir -p CH_insane5
cd CH_insane5

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
cp /home/birac/Desktop/CoarseGrainSims/KIT/DUPC_single.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/insane_CFTR.py .
#<<'END'
rm -f ${prot}-CG.pdb -f $prot-mem.gro *tpr

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

python insane_CFTR.py -f $prot-CG.pdb -o $prot-mem.gro  -x 25 -y 25 -z 18 -l DPPC:42 -l DUPC:28 -l CHOL:30 -sol W -salt 0.15 -op 500  > insaneout

if [ ! -f $prot-mem.gro ]; then
        echo box2 fails
        exit 0
fi

echo "box building done!"

echo "
" >> ${prot}.top
cat ${prot}.top insanetop.txt > new.top

<<'END'
PCnum=`grep -c "DPPC   PO4" $prot-mem.gro`
DLnum=`grep -c "DUPC   PO4" $prot-mem.gro`
CHnum=`grep -c "CHOL   ROH" $prot-mem.gro`

echo "
DPPC $PCnum" >> ${prot}.top
echo "DUPC $DLnum" >> ${prot}.top
echo "CHOL $CHnum" >> ${prot}.top


genbox -o sol.gro -cs waterbox.gro -cp $prot-mem.gro -vdwd 0.21 -p ${prot}.top >& sol #-maxsol 768
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
 
echo -e  '"'DPPC'"' '|' '"'DUPC'"' '|' '"'CHOL'"' '\n' '"'W'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | make_ndx -f sol.gro -o memsol1 -n index >& ndx3

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

cp $prot-mem.gro ions.gro
mv new.top $prot.top 
sed 's/NA+/NA /g' ions.gro -i
sed 's/NA+/NA /g' $prot.top -i

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
echo -e  '"'DPPC'"'  '|' '"'DUPC'"' '|' '"'CHOL'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q | make_ndx -f ions.gro -o memsol -n index >& ndx4

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
