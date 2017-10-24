#!/bin/bash

# 128 lipid molecules in total.
# build water box

prot1=$1
lipid1=DSPE # these don't currently do anything 
lipid2=CL   #
CD=`pwd`

rm -f vacuum_min.gro PE.gro PECL.gro sol.gro ions.gro vacuum_min.tpr genion.tpr min_test.tpr

prot=${prot1//.pdb}

cp /home/birac/Desktop/CoarseGrainSims/KIT/martinize.py .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization-vacuum.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_1.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_2.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_3.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/*itp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/DSPE_single.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/CL_single.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/water.gro .

rm -f ${prot}-CG.pdb

python martinize.py -f ${prot}.pdb -o ${prot}.top -p ${prot}-posres.itp -x ${prot}-CG.pdb -dssp mkdssp -p backbone -ff elnedyn22 >& one

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
sed '1i #include "martini_v2.0_lipids_cardiolipin.itp"' $prot.top -i
sed '1i #include "martini_v2.2_aminoacids.itp"' $prot.top -i
sed '1i #include "martini_v2.2.itp"' $prot.top -i

editconf -f ${prot}-CG.pdb -o ${prot}-CG.box.pdb -box 15 15 15 >&  edc

rm -f vacuum_min.tpr

grompp -f minimization-vacuum.mdp -p ${prot}.top -c ${prot}-CG.box.pdb -o vacuum_min >& g1
mdrun -v -deffnm vacuum_min >& md1

if [ ! -f vacuum_min.tpr ]; then
        echo "
        ##################
        Error in vacuum MD
        ##################"
        exit 0
else
        echo vacuum min done
fi

# Add 460 DSPE

genbox -cp vacuum_min.gro -ci DSPE_single.gro -nmol 460 -try 500 -o PE.gro -vdwd 0.21 >& box1

# add 52 CL
genbox -cp PE.gro -ci CL_single.gro -nmol 52 -try 500 -o PECL.gro -vdwd 0.21 >& box2

sed '/DSPE/d' ${prot}.top -i 
sed '/CARD/d' ${prot}.top -i
sed '/W/d' ${prot}.top -i

echo "
" >> ${prot}.top

PEnum=`grep -c "DSPE   PO4" PECL.gro`
CLnum=`grep -c "CARD   GL0" PECL.gro`

echo "DSPE $PEnum" >> ${prot}.top
echo "CARD $CLnum" >> ${prot}.top

genbox -o sol.gro -cs water.gro -cp PECL.gro -vdwd 0.21 >& sol #-maxsol 768

sed '/W/d' ${prot}.top -i
SOLnum=`grep -c "W" sol.gro`
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
 
echo -e  '"'DSPE'"' '|' '"'CARD'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | make_ndx -f sol.gro -o memsol1 -n index >& ndx3

rm -f genion.tpr
rm -f ions.gro

grompp -f minimization.mdp -p ${prot}.top -c sol.gro -o genion.tpr -n memsol1 >& g2

echo W | genion -s genion.tpr -pname NA+ -nname CL- -neutral -o ions.gro >& ions

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

sed '/W/d' ${prot}.top -i
SOLnum=`grep -c "W" ions.gro`
echo "W $SOLnum" >> ${prot}.top

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

echo -e  '"'DSPE'"' '|' '"'CARD'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q | make_ndx -f ions.gro -o memsol2 -n index >& ndx4

sed '/NA/d' ${prot}.top -i
NAnum=`grep -c "NA" ions.gro`
echo "NA $NAnum" >> ${prot}.top

sed '/CL-/d' ${prot}.top -i
CLnum=`grep -c "CL-" ions.gro`
echo "CL- $CLnum" >> ${prot}.top

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

echo -e  '"'DSPE'"' '|' '"'CARD'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | make_ndx -f ions.gro -o memsol2 -n index >& ndx5

grompp -f minimization.mdp -p ${prot}.top -c ions.gro -o min_test.tpr -n memsol2 >& gp3

$GMX51 mdrun -c -deffnm min_test >& md2

if [[ ! -f min_test.gro ]]; then
	echo min test failed
	exit 0
else
	echo min test worked!
fi

rm -f minimization.tpr equilibration_1.tpr equilibration_2.tpr 

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

echo -e  '"'DSPE'"' '|' '"'CARD'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | make_ndx -f ions.gro -o memsol -n index >&ndx6

head -n 2 ions.gro > all.gro
grep "BB\|SC" ions.gro >> all.gro 
grep "DSPE  " ions.gro >> all.gro
grep "CARD  " ions.gro >> all.gro
grep "W  " ions.gro >> all.gro
grep "NA  " ions.gro >> all.gro
grep "CL  " ions.gro >> all.gro
tail -n 1 ions.gro >> all.gro

# Do a short minin with posres - need define = -DPOSRES in mdp

grompp -f minimization.mdp -c all.gro -p ${prot}.top -o minimization.tpr -n memsol.ndx >& gp4

if [ ! -f minimization.tpr ]; then
        echo "
        #################################
        Error in building minimization !!
        #################################"
        exit 0
fi


$GMX51 mdrun -deffnm minimization -v >& md3

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

$GMX51 grompp -f equilibration_1.mdp -c minimization.gro -p ${prot}.top -o equilibration_1.tpr -n memsol.ndx >& gp5
$GMX51 mdrun -deffnm equilibration_1 -v >& md4  # Takes about 25 mins on 1 CPU 

$GMX51 grompp -f equilibration_2.mdp -c equilibration_1.gro -p ${prot}.top -o equilibration_2.tpr -n memsol.ndx >& gp6
$GMX51 mdrun -deffnm equilibration_2 -v >& md5  # Takes about 150 minutes on 1 CPU

$GMX51 grompp -c equilibration_2.gro -f equilibration_3.mdp -p ${prot}.top -o equilibration_3 >& gp7

echo "

Check orientation with checkbox.sh!

"
