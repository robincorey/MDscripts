#!/bin/bash
# Finally, a setup script for martini, wahey!!
#
# Tailored for 3DIN with card+insane
#

#import required files

rm *mdp insane.py insane_rc.py martinize.py

cp /home/birac/Desktop/CoarseGrainSims/KIT/martinize.py .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization-vacuum.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/minimization.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_1.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/equilibration_2.mdp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/cgmd_1000.mdp .
#cp /home/birac/Desktop/CoarseGrainSims/KIT/waterbox.gro .
cp /home/birac/Desktop/CoarseGrainSims/KIT/*itp .
cp /home/birac/Desktop/CoarseGrainSims/KIT/insane_rc.py .

#
## Define variables
#

prot=$1

#
## Firstly, make the CG protein
## Use -p to specify posres
#

mv $prot-CG.pdb $prot-CG_backup.pdb 

python martinize.py -f $prot.pdb -o $prot.top -p $prot-posres.itp -x $prot-CG.pdb -dssp mkdssp -p backbone -ff martini22

#
## Test to see if martinize worked
#

if [ ! -f $prot-CG.pdb ]; then
	echo "
	#####################################
	Error in martinize - PDB not created!
	#####################################"
	exit 0
fi

#
## Amend topology file accordingly
#

sed -i -e "1,2d" $prot.top
sed '1i #include "martini_v2.0_ions.itp"' $prot.top -i
sed '1i #include "martini_v2.0_lipids_cardiolipin.itp"' $prot.top -i 
sed '1i #include "martini_v2.2_aminoacids.itp"' $prot.top -i
sed '1i #include "martini_v2.2.itp"' $prot.top -i

#
## Build into a membrane with insane.py
#

mv $prot-mem.gro $prot-mem_backup.gro 

### Insane membrane construction using input from above

lr1=$(expr $lipidratio1 / 1)
lr2=$(expr $lipidratio2 / 1)

mv insaneout insaneout_backup

if [ -z "$lipid2" ];
then
	`python insane_rc.py -f $prot-CG.pdb -o $prot-mem.gro  -x 17 -y 17 -z 15 -l $lipid1:$lr1 -sol W -salt 0.15 -dm 3 -op 500` > insaneout
else
	`python insane_rc.py -f $prot-CG.pdb -o $prot-mem.gro  -x 17 -y 17 -z 15 -l $lipid1:$lr1 -l $lipid2:$lr2 -sol W -salt 0.15 -dm 3 -op 500` > insaneout
fi	


if [ ! -f $prot-mem.gro  ]; then
        echo "
        ######################
        Error in INSANE.py!!!!
        ######################"
        exit 0
fi

<< 'END'
# Update topology file

sed -e '/ W   /d' $prot.top  -i
sed -e '/NA   /d' $prot.top  -i
sed -e '/CL   /d' $prot.top  -i


solnum=`grep -c "W " $prot-mem.gro`
nanum=`grep -c "NA+" $prot-mem.gro`
clnum=`grep -c "CL-" $prot-mem.gro`

echo -e '\n' >> $prot.top
echo 'W   ' `echo $solnum` >> $prot.top
echo 'NA+  ' `echo $nanum` >> $prot.top
echo 'CL-  ' `echo $clnum` >> $prot.top

END

# To break up the colour!

## Section to generate new topology file

echo -e '\n' >> $prot.top
cat insanetop.txt >> $prot.top # Insanetop created by modified insane.py

<< 'END'
mv minimization-vacuum.gro minimization-vacuum_backup.gro

grompp -f minimization-vacuum.mdp -p $prot.top -c $prot-mem.gro -o minimization-vacuum.tpr
mdrun -deffnm minimization-vacuum -v


if [ ! -f minimization-vacuum.gro ]; then
        echo "
        ###############################
        Error in vacuum minimization !!
        ###############################"
        exit 0
fi

END

# Make an appropriate index file to combine all lipids and solvents

rm memsol.ndx
rm index.ndx
rm ndx
rm ndx2

echo q | make_ndx -f $prot-mem.gro -o index > ndx

cat ndx | grep [0-9] | grep atoms | awk '{print $2}' | sort | uniq -d > ndx2

if [[ -s ndx2 ]] ; then
        while read i; do
                sed '0,/'$i'/! s/'$i'/Dud/' index.ndx -i
        done < ndx2
fi

#echo q | make_ndx -f $prot-mem.gro -o index.ndx

groups=`grep -c '\[' index.ndx`

mem=$(expr $groups + 0)
sol=$(expr $groups + 1)

echo -e '"'$lipid1'"' '|' '"'$lipid2'"' '\n' '"'W'"' '|' '"'ION'"' '\n' name $mem Mem '\n' name $sol Sol '\n'q  | make_ndx -f $prot-mem.gro -o memsol.ndx -n index.ndx

# Do a short minin with posres - need define = -DPOSRES in mdp

grompp -f minimization.mdp -c $prot-mem.gro -p $prot.top -o minimization.tpr -n memsol.ndx

if [ ! -f minimization.tpr ]; then
        echo "
        #################################
        Error in building minimization !!
        #################################"
        exit 0
fi


mdrun -deffnm minimization -v

if [ ! -f minimization.gro ]; then
        echo "
        ########################
        Error in minimization !!
        ########################"
        exit 0
fi

# Then crack out an equilibration run

grompp -f equilibration_1.mdp -c minimization.gro -p $prot.top -o equilibration_1.tpr -n memsol.ndx
mdrun -deffnm equilibration_1 -v   # Takes about 25 mins on 1 CPU

grompp -f equilibration_2.mdp -c equilibration_1.gro -p $prot.top -o equilibration_2.tpr -n memsol.ndx
mdrun -deffnm equilibration_2 -v   # Takes about 150 minutes on 1 CPU

