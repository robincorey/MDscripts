#!/bin/bash
# To run protein folding simulations

############################
prot=input_RQH          #### >>> Needs to match input PDB
forcefield=opls     	#### >>> amber, charmm, opls 
############################

# * OPLS-AA
# ** Amberff99
# *** Charmm27

###### Change log 
###### MODIFIED 24/11/2014 for new gmembed - including NUC in gmembed option###############
###### MODIFIED 04/06/2015 for Gromos ff and new lipids ##################################
###### ADDED OPLS on 16/06/2015 ############
###### ADDED CHARMM 36 22/06/2015 ############
###### Made into prot dolfing script on 21/09/2015


###############################################################################

#
## Shouldn't have to do anything after this bit
#

cp /home/birac/Desktop/MD_KIT/em.mdp .
cp /home/birac/Desktop/MD_KIT/md_fold.mdp .
cp /home/birac/Desktop/MD_KIT/eq_fold1.mdp .
cp /home/birac/Desktop/MD_KIT/eq_fold2.mdp .
cp /home/birac/Desktop/MD_KIT/eq_npt.mdp .


## Next line checks if forcefield is correctly entered
## If $forcefield is not opls and is not gromos then exit

if [ ! "$forcefield" == "opls" ] && [ ! "$forcefield" == "amber" ] && [ ! "$forcefield" == "charmm" ]
        then echo "
####################################
Not a forcefield - check your input!
####################################"; exit 0
fi

if [ ! -f ${prot}.pdb ]
	then echo" 
##########################
No PDB - check your input!
##########################"; exit 0
fi

firstres=`grep ATOM ${prot}.pdb | head -n 1 | awk '{print $5}'`

if [ ! ${firstres} -eq 1 ]
        then echo "
#####################################
PDB numbering off - check your input!
#####################################"; exit 0
fi

## To stop shit building up

mkdir Overwrites
mv *#* Overwrites

###############################################################################

#
## Main section. Numbered for ease
#

#
## 1;; Generate topology file from pdb
#

mv $prot.gro $prot_backup.gro

if [[ "$forcefield" == "opls" ]]
	then ff=oplsaa
elif [[ "$forcefield" == "charmm" ]]
	then ff=charmm27
elif [[ "$forcefield" == "amber" ]]
        then ff=amber99
fi

# add protonation options if required

pdb2gmx -f $prot.pdb -o $prot.gro -p $prot.top -ff $ff -water tip3p -ignh    	

if [ ! -f $prot.gro ]; then
        echo "
        ################
        Error in pdb2gmx
        ################"
        exit 0
fi

#
## 2;; Building simulation box and preparing for g_membed
#

# Do a very short minimization of protein to prevent anything terrible happening
# Mostly useful for mutations etc.

#grompp -f em_short.mdp -c $prot.gro -p $prot.top -o em_pre
#mdrun -v -deffnm em_pre

#editconf -f em_pre.gro -o $prot.gro.pdb
editconf -f ${prot}.gro -o ${prot}.gro.pdb
sed -e 's/TER/   /g' -e 's/ENDMDL/      /g' ${prot}.gro.pdb -i

#
## Build simulation box. 
#

editconf -f ${prot}.gro.pdb -o ${prot}.pdb.gro -d 1 -bt cubic

if [ ! -f ${prot}.pdb.gro ]; then
        echo "
        ###############################
        Error in making simulation box!
        ###############################"
        exit 0
fi

size=`tail -n 1 ${prot}.pdb.gro`
grep -v "^$" ${prot}.pdb.gro > prot.gro


#
## 3;; Solvating protein
#

mv solmem.gro solmem_backup.gro

genbox -cp prot.gro -o solprot.gro -p $prot.top -cs

if [ ! -f solprot.gro ]; then
        echo "
        #########################
        Error in solvating system
        #########################"
        exit 0
fi

editconf -f solprot.gro -o solprot.pdb


mv ions_in.tpr ions_in_backup.tpr

grompp -f em.mdp -c solprot.gro -p $prot.top -o ions_in

if [ ! -f ions_in.tpr ]; then
        echo "
        ########################
        Error in building genion
        ########################"
        exit 0
fi

echo SOL | genion -s ions_in.tpr -p $prot.top -o ions.gro -neutral -conc 0.150

editconf -f ions.gro -o ions.pdb

#
## 4;; Make a new gro file with everything in order, yo
#

grep SOL ions.pdb> SOLx.pdb
grep NA ions.pdb> NAx.pdb
grep CL ions.pdb> CLx.pdb
grep -e ALA -e CYS -e ASP -e GLU -e PHE -e GLY -e HIS -e ILE -e LYS -e LEU -e MET -e ASN -e PRO -e GLN -e ARG -e SER -e THR -e VAL -e TRP -e TYR ions.pdb | grep ATOM > protein.pdb

editconf -f protein.pdb -o protein.gro
editconf -f SOLx.pdb -o SOLx.gro
editconf -f NAx.pdb -o NAx.gro
editconf -f CLx.pdb -o CLx.gro

mv all.pdb all_backup.pdb

cat protein.pdb SOLx.pdb NAx.pdb  CLx.pdb | grep -e ALA -e CYS -e ASP -e GLU -e PHE -e GLY -e HIS -e ILE -e LYS -e LEU -e MET -e ASN -e PRO -e GLN -e ARG -e SER -e THR -e VAL -e TRP -e TYR  -e SOL -e NA -e CL >> all.pdb

editconf -f all.pdb -o all.gro

if [[ ! -f all.gro ]]; then
	echo "
        #######################
        Error in making all.gro
        #######################"
	exit 0
fi

sed -i '$ d' all.gro
tail -n 1 solprot.gro > bottom
cat bottom >> all.gro

#
## 5;; Update topology for new components 
#

## Updating solvent section. Na and Cl only

sed -i '/SOL/d' $prot.top
uniqatomssol=`grep SOL all.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep SOL all.gro | awk '{print $1}' | wc -l`
echo SOL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> $prot.top

sed -i '/NA/d' $prot.top
uniqatomssol=`grep NA all.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep NA all.gro | awk '{print $1}' | wc -l`
echo NA `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> $prot.top

sed -i '/CL/d' $prot.top
uniqatomssol=`grep CL all.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep CL all.gro | awk '{print $1}' | wc -l`
echo CL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> $prot.top

#
## 9;; Minimizing and solvating system, yo
#

mv em1.tpr em1_backup.tpr

grompp -v -f em.mdp -o em1 -e em1 -c all.gro  -p $prot.top 

if [ ! -f em1.tpr ]; then
        echo "
        #######################
        Error in building min 1
        #######################"
        exit 0
fi

mv after_em1.gro after_em1_backup.gro

mdrun -v -s em1 -o em1 -e e1 -c after_em1.gro 

editconf -f after_em1.gro -o after_em1.pdb

if [ ! -f after_em1.gro ]; then
        echo "
        #######################
        Error in minimization 1
        #######################"
        exit 0
fi


#
## 10;; Prepare for equilibration
#

# short timestep run to ease in 

grompp -f eq_fold1.mdp -c after_em1.gro -p $prot.top -o eq1.tpr

mdrun -v -deffnm eq1

# nvt - longer nvt run

grompp -f eq_fold2.mdp -c eq1.gro -p $prot.top -o eq2.tpr

mdrun -v -deffnm eq2

# npt - add pressure

grompp -f eq_npt.mdp -c eq2.gro -p $prot.top -o eq3.tpr

mdrun -v -deffnm eq3 -c ${prot}_${forcefield}.gro

editconf -f ${prot}_${forcefield}.gro -o ${prot}_${forcefield}.pdb

grompp -f md_fold.mdp -c ${prot}_${forcefield}.gro -p $prot.top -o ${forcefield}${prot}.tpr

# Analyse eq data

#echo -e Potential'\n'0 | g_energy -f eq1.edr -o potential1.xvg
#echo -e Potential'\n'0 | g_energy -f eq2.edr -o potential2.xvg
#echo -e Potential'\n'0 | g_energy -f eq3.edr -o potential3.xvg

#echo -e Temperature'\n'0 | g_energy -f eq1.edr -o temp1.xvg
#echo -e Temperature'\n'0 | g_energy -f eq2.edr -o temp2.xvg
#echo -e Temperature'\n'0 | g_energy -f eq3.edr -o temp3.xvg

mkdir Backup
mv *backup* Backup/

