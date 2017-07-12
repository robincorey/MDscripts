#!/bin/bash
# Script to run a membrane protein from pdb state to membrane embedded state
# Need to make sure pdb file is numbered correctly and TER is used to break chains
# Modified for SecA inclusion - specific for this system but should work for any prot
# If no nucleotide might cause issues - need to check

#
## Need to define prot and nuc options (ATP/ADP)
## Check variables are correctly assigned!!
#

############################
prot=5aww_MTS		#### >>> Needs to match input PDB
nuc=			#### >>> ADP or ATP (if none leave blank)
forcefield=opls      	#### >>> Either opls* or gromos** at the moment
gromoslipidonly=popc	#### >>> POPE, POPC and POPG all working
z=10			#### >>> Simulation box z axis. 15 as standard
sol=spc			#### >>> spc, tip3p etc
ters=no 		#### >>> Use for capping peptides
propka=no		#### >>> do propka?
############################

# * OPLS-AA
# ** Gromos53a6

if [[ "$1" == "-h" ]]; then echo "
 2;; Building simulation box and preparing for g_membed"; sleep 0.2; echo "
 3;; Solvating membrane (Opls only)"; sleep 0.2; echo "
 4;; Make a new gro file with everything in order, yo"; sleep 0.2; echo "
 5;; Update topology for new components"; sleep 0.2; echo "
 6;; Make an appropriate index file to combine all lipids and solvents"; sleep 0.2; echo "
 7;; The DREADED g_membed stage - this is where most setups fail as membrane just explodes grrr."; sleep 0.2; echo "
 8;; Reupdate topogy, yo"; sleep 0.2; echo "
 9;; Minimizing and solvating system, yo"; sleep 0.2; echo "
 10;; Prepare for equilibration"; sleep 0.2; echo "
"
exit 0
fi

###### Change log 
###### MODIFIED 24/11/2014 for new gmembed - including NUC in gmembed option###############
###### MODIFIED 04/06/2015 for Gromos ff and new lipids ##################################
###### ADDED OPLS on 16/06/2015 ############
###### ADDED CHARMM 36 22/06/2015 ############
###### Added peptide capping stuff 08/03/2016 #############
###### Added propka stuff 28/02/2017 #############

if [[ "$forcefield" == "gromos" ]]
        then lipid=$gromoslipidonly	## Can change if using ff with different lipids
fi

###############################################################################

#
## Shouldn't have to do anything after this bit
#


# get correct membrane and other setup files

if [[ "$forcefield" == "opls" ]]
	then cp /home/birac/Desktop/MD_KIT/popc_opls.itp . 
	cp /home/birac/Desktop/MD_KIT/512-popc-eq.pdb .
	cp /home/birac/Desktop/MD_KIT/*mdp .
	cp /home/birac/Desktop/MD_KIT/membed.dat .
fi

if [[ "$forcefield" == "gromos" ]]
	then cp /home/birac/Desktop/Gromos53A6/MD_KIT/*itp . 
	cp /home/birac/Desktop/Gromos53A6/MD_KIT/512-$lipid-gromos53a6.pdb .
	cp /home/birac/Desktop/Gromos53A6/MD_KIT/*mdp .
	cp /home/birac/Desktop/Gromos53A6/MD_KIT/membed.dat .
fi


## Next line checks if forcefield is correctly entered
## If $forcefield is not opls and is not gromos then exit

if [ ! "$forcefield" == "opls" ] && [ ! "$forcefield" == "gromos" ]
        then echo "
####################################
Not a forcefield - check your input!
####################################"; exit 0
fi

## To stop shit building up

mkdir Overwrites
mv *#* Overwrites

###############################################################################

#
## Main section. Numbered for ease
#

#
## 2;; Building simulation box and preparing for g_membed
#

editconf -f ${prot}.gro -o ${prot}.gro.pdb
sed -e 's/TER/   /g' -e 's/ENDMDL/      /g' ${prot}.gro.pdb -i

mv protmem.gro.pdb protmem.gro_backup.pdb

if [[ "$forcefield" == "opls" ]]
        then cat ${prot}.gro.pdb 512-popc-eq.pdb > protmem.gro.pdb
elif [[ "$forcefield" == "gromos" ]]
        then cat ${prot}.gro.pdb 512-$lipid-gromos53a6.pdb > protmem.gro.pdb
fi

if [ ! -f protmem.gro.pdb ]; then
        echo "
        ###############################
        Error in building prot-membrane
        ###############################"
        exit 0
fi

mv ${prot}-mem.pdb.gro ${prot}-mem_backup.pdb.gro

#
## Build simulation box. X+Y based on input pdb. Z based on input
#

if [ "$forcefield" == "gromos" ] && [ "$lipid" == "pope" ]
        then editconf -f protmem.gro.pdb -o ${prot}-mem.pdb.gro -c -box 12.2528 12.3878 $z -bt triclinic
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
	then editconf -f protmem.gro.pdb -o ${prot}-mem.pdb.gro -c -box 13.2756 13.2546 $z -bt triclinic
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popc" ]
        then editconf -f protmem.gro.pdb -o ${prot}-mem.pdb.gro -c -box 13.3356 13.3040 $z -bt triclinic
elif [[ "$forcefield" == "opls" ]]
	then editconf -f protmem.gro.pdb -o ${prot}-mem.pdb.gro -c -box 13.1708 13.1708 $z -bt triclinic
fi

if [ ! -f ${prot}-mem.pdb.gro ]; then
        echo "
        ###############################
        Error in making simulation box!
        ###############################"
        exit 0
fi

size=`tail -n 1 ${prot}-mem.pdb.gro`
grep -v "^$" ${prot}-mem.pdb.gro > protmem.gro

#
## Update topology
#

if [ "$forcefield" == "gromos" ] && [ "$lipid" == "pope" ]
	then sed '/Include water topology/i \
# include "'pope_GROMOS-CKP.itp'"\ ' ${prot}.top -i
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
        then sed '/Include water topology/i \
# include "'dpopgGromacs4_53a6.itp'"\ ' ${prot}.top -i ; sed '/# include "'dpopgGromacs4_53a6.itp'"/i \
# include "'lpopgGromacs4_53a6.itp'"\
' ${prot}.top -i
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popc" ]
	then sed '/Include water topology/i \
# include "'popcGromacs4_53a6.itp'"\ ' ${prot}.top -i
elif [[ "$forcefield" == "opls" ]]
	then sed '/Include water topology/i \
# include "'popc_opls.itp'"\
' ${prot}.top -i
fi

#
## 3;; Solvating membrane (Opls only)
#

mv solmem.gro solmem_backup.gro

if [[ "$forcefield" == "opls" ]]
	then genbox -cp protmem.gro -o solmem.gro -p ${prot}.top -cs
fi

if [ "$forcefield" == "opls" ] && [ ! -f solmem.gro ]; then
        echo "
        ################################
        Error in solvating opls membrane
        ################################"
        exit 0
fi

mv premembed.gro premembed_backup.gro

if [[ "$forcefield" == "opls" ]]
	then cp solmem.gro premembed.gro
fi

#
## 4;; Make a new gro file with everything in order, yo
#

if [[ "$lipid" == "popg" ]]
	then sed -i 's/NA+/NA /g' premembed.gro
	sed -i 's/CL-/CL /g' premembed.gro
fi

mv all.gro all_backup.gro
mv premembed.pdb premembed_backup.pdb
if [[ ! "$forcefield" == "opls" ]]
        then cp protmem.gro premembed.gro
fi

editconf -f premembed.gro -o premembed.pdb

if [ ! -f premembed.pdb ]; then
        echo "
        ################################
        Error in solvating opls membrane
        ################################"
        exit 0
fi

if [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
	then grep -e "DPO\|LPO" premembed.pdb> Memx.pdb
else grep POP premembed.pdb> Memx.pdb
fi

grep SOL premembed.pdb> SOLx.pdb
grep NA premembed.pdb> NAx.pdb
grep CL premembed.pdb> CLx.pdb
if [[ ! -z "${nuc}" ]]; then
grep ${nuc} premembed.pdb> Nucx.pdb; fi
grep MG premembed.pdb> MGx.pdb
grep -e ACE -e NH2 -e NAC -e ALA -e CYS -e ASP -e GLU -e PHE -e GLY -e HIS -e ILE -e LYS -e LEU -e MET -e ASN -e PRO -e GLN -e ARG -e SER -e THR -e VAL -e TRP -e TYR premembed.pdb | grep ATOM > protein.pdb

editconf -f protein.pdb -o protein.gro
editconf -f Memx.pdb -o Memx.gro
editconf -f SOLx.pdb -o SOLx.gro
editconf -f NAx.pdb -o NAx.gro
editconf -f CLx.pdb -o CLx.gro
editconf -f Nucx.pdb -o Nucx.gro
editconf -f MGx.pdb -o MGx.gro

mv all.pdb all_backup.pdb

### This section might seem superfluous but is here to make sure gro file is all in block order (i.e. easier to topologise)

cat MGx.pdb Nucx.pdb protein.pdb Memx.pdb SOLx.pdb NAx.pdb  CLx.pdb | grep -e MG -e ATP -e ADP -e ACE -e ALA -e CYS -e ASP -e GLU -e PHE -e GLY -e HIS -e ILE -e LYS -e LEU -e MET -e ASN -e PRO -e GLN -e ARG -e SER -e THR -e VAL -e TRP -e TYR -e NAC -e NH2 -e POP -e DPO -e LPO -e DDM -e SOL -e NA -e CL >> all.pdb

editconf -f all.pdb -o all.gro

if [[ ! -f all.gro ]]; then
	echo "
        #######################
        Error in making all.gro
        #######################"
	exit 0
fi

sed -i '$ d' all.gro
tail -n 1 premembed.gro > bottom
cat bottom >> all.gro

#
## 5;; Update topology for new components 
#

if [[ "$forcefield" == "opls" ]]
	then sed -i '/POP/d' ${prot}.top 
	lipidcount=`grep -c P8 all.gro`
	echo POP $lipidcount >> ${prot}.top
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popc" ]
	then sed -i '/POPC/d' ${prot}.top
	lipidcount=`grep -c P8 all.gro`
	echo POPC $lipidcount >> ${prot}.top
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
        then sed -i '/DPO/d' ${prot}.top 
	sed -i '/LPO/d' ${prot}.top
	lpocount=`grep -c "LPO     P9" all.gro`
	dpocount=`grep -c "DPO     P9" all.gro`
	echo DPO $dpocount >> ${prot}.top 
	echo LPO $lpocount >> ${prot}.top
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "pope" ]
        then sed -i '/POPE/d' ${prot}.top 
	lipidcount=`grep -c P8 all.gro` 
	echo POPE $lipidcount >> ${prot}.top
fi

### Above section needs to be checked for each forcefield and lipid to make sure still works...

## Updating solvent section. Na and Cl only

sed -i '/SOL/d' ${prot}.top
uniqatomssol=`grep SOL all.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep SOL all.gro | awk '{print $1}' | wc -l`
echo SOL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top

# no ions yet?

<<'END'
sed -i '/NA/d' ${prot}.top
uniqatomssol=`grep NA all.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep NA all.gro | awk '{print $1}' | wc -l`
echo NA `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top

sed -i '/CL/d' ${prot}.top
uniqatomssol=`grep CL all.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep CL all.gro | awk '{print $1}' | wc -l`
echo CL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top
END
#
## 6;; Make an appropriate index file to combine all lipids and solvents
#

cp all.gro solmem.gro

rm memsol.ndx
rm index.ndx
rm ndx
rm ndx2

echo q | make_ndx -f solmem.gro -o index > ndx

cat ndx | grep [0-9] | grep atoms | awk '{print $2}' | sort | uniq -d > ndx2

if [[ -s ndx2 ]] ; then
        while read i; do
                sed '0,/'$i'/! s/'$i'/Dud/' index.ndx -i
        done < ndx2
fi

## Below bit reads current index file to guess future number of Mem and Nuc index
## Ugly, but effective (a bit like your mum)
## Only for Nuc sims

groups=`grep -c '\[' index.ndx`

function nucndx {

memno=$(expr $groups + 0)
nucno=$(expr $groups + 1)

if [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
        then echo -e '"DPO"' '|' '"LPO"' '\n' '"MG"' '|' '"'$nuc'"' '|' '"Protein"' '\n' name $memno Mem '\n' name $nucno Nuc '\n'q  | make_ndx -f solmem.gro -o memsolpre.ndx -n index.ndx
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "pope" ]
        then echo -e '"POP"' '\n' '"MG"' '|' '"'$nuc'"' '|' '"Protein"' '\n' name $memno Mem '\n' name $nucno Nuc '\n'q  | make_ndx -f solmem.gro -o memsolpre.ndx -n index.ndx
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popc" ]
        then echo -e '"POP"' '\n' '"MG"' '|' '"'$nuc'"' '|' '"Protein"' '\n' name $memno Mem '\n' name $nucno Nuc '\n'q  | make_ndx -f solmem.gro -o memsolpre.ndx -n index.ndx
elif [[ "$forcefield" == "opls" ]]
	then echo -e '"POP"' '\n' '"MG"' '|' '"'$nuc'"' '|' '"Protein"' '\n' name $memno Mem '\n' name $nucno Nuc '\n'q  | make_ndx -f solmem.gro -o memsolpre.ndx -n index.ndx
fi
}

### make ndx file either way

if [[ ! -z ${nuc} ]]; then
	nucndx
else 
	memno=$(expr $groups + 0)
	echo -e '"POP"' '\n' name $memno Mem'\n'q | make_ndx -f solmem.gro -o memsolpre.ndx -n index.ndx
fi

## Could consolidate some of the above loops but useful to keep them separated for future use

groups2=`grep -c '\[' memsolpre.ndx`
memsolno=$(expr $groups2 + 0)

if [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
	then echo -e '"Mem"' '|' '"Water_and_ions"' '\n'q | make_ndx -f solmem.gro -n memsolpre.ndx -o memsol.ndx
	sed 's/Mem_Water_and_ions/Mem_SOL/' memsol.ndx -i
else echo -e '"Mem"' '|' '"SOL"' '\n'q | make_ndx -f solmem.gro -n memsolpre.ndx -o memsol.ndx
fi
#
## 7;; The DREADED g_membed stage - this is where most setups fail as membrane just explodes grrr.
#

mv gmembed.tpr gmembed_backup.tpr

# Below section adds extra coupling group for if ions present

if [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
	then sed -i 's/Nuc Mem SOL/Nuc Mem SOL Ion/' gmembed_nuc.mdp
	sed -i 's/1 1 1/1 1 1 1/' gmembed_nuc.mdp
	sed -i 's/300 300 300/300 300 300 300/' gmembed_nuc.mdp
fi

## for awkward bonded interactions in ACE stuff

if [[ $ters == "yes" ]]
	then sed -i 's/1     5     7     9     3/;1     5     7     9     3/g' ${prot}.top
fi

#End

### Here the maxwarn is justifiable as sometimes the toplogy will have blank NA+CL lines

if [[ ! -z "${nuc}" ]]; then
grompp -v -f gmembed_nuc.mdp -o gmembed -c solmem.gro -p ${prot}.top -n memsol.ndx -maxwarn 4
fi

if [[ -z "${nuc}" ]]; then
grompp -v -f gmembed.mdp -o gmembed -c solmem.gro -p ${prot}.top -n memsol.ndx -maxwarn 4
fi

if [ ! -f gmembed.tpr ]; then
        echo "
        ###############################
        Error in making gmembed.tpr :-(
        ###############################"
        exit 0
fi

mv membedded.gro membedded_backup.gro

if [[ ! -z "${nuc}" ]]; then
echo Nuc Mem | mdrun -s gmembed.tpr -membed membed.dat -o traj.trr -c membedded.gro -e ener.edr -cpt -1 -v -stepout 100 -mn memsol
elif [[ -z "${nuc}" ]]; then
echo Protein Mem | mdrun -s gmembed.tpr -membed membed.dat -o traj.trr -c membedded.gro -e ener.edr -cpt -1 -v -stepout 100 -mn memsol
fi

if [ ! -f membedded.gro ]; then
        echo "
        ######################
        Error in doing gmembed
        ######################"
        exit 0
fi

#
## Embedding done!
#

#
## 8;; Reupdate topogy, yo
#

if [[ "$forcefield" == "opls" ]]
        then sed -i '/POP/d' ${prot}.top
        lipidcount=`grep -c P8 membedded.gro`
        echo POP $lipidcount >> ${prot}.top
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popg" ]
        then sed -i '/DPO/d' ${prot}.top
        sed -i '/LPO/d' ${prot}.top
        lpocount=`grep -c "LPOPG   P9" membedded.gro`
        dpocount=`grep -c "DPOPG   P9" membedded.gro`
        echo DPO $dpocount >> ${prot}.top
        echo LPO $lpocount >> ${prot}.top
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "pope" ]
        then sed -i '/POPE/d' ${prot}.top
        lipidcount=`grep -c P8 membedded.gro`
        echo POPE $lipidcount >> ${prot}.top
elif [ "$forcefield" == "gromos" ] && [ "$lipid" == "popc" ]
        then sed -i '/POPC/d' ${prot}.top
        lipidcount=`grep -c P8 membedded.gro`
        echo POPC $lipidcount >> ${prot}.top
fi

sed -i '/SOL/d' ${prot}.top
uniqatomssol=`grep SOL membedded.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep SOL membedded.gro | awk '{print $1}' | wc -l`
echo SOL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top

sed -i '/NA/d' ${prot}.top
uniqatomssol=`grep NA membedded.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep NA membedded.gro | awk '{print $1}' | wc -l`
echo NA `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top

sed -i '/CL/d' ${prot}.top
uniqatomssol=`grep CL membedded.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep CL membedded.gro | awk '{print $1}' | wc -l`
echo CL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top


#
## 9;; Minimizing and solvating system, yo
#

mv em1.tpr em1_backup.tpr

grompp -v -f em.mdp -o em1 -e em1 -c membedded.gro  -p ${prot}.top -maxwarn 3

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


mv em_sol.gro em_sol_backup.gro

grep -Ev 'SOL|NA|CL' after_em1.gro > presol.gro


tail -n+3 presol.gro > presol_mod.gro
atomnumbercount1=`grep -c [0-9] presol_mod.gro | awk '{print $1 }'`
atomnumbercount2=`expr $atomnumbercount1 - 1`
sed '1i '$atomnumbercount2''  presol_mod.gro -i
sed '1i Membrane protein whaaat no solvent??' presol_mod.gro -i
#sed -i '$ d' all.gro

sed -i '/SOL/d' ${prot}.top
sed -i '/NA/d' ${prot}.top
sed -i '/CL/d' ${prot}.top

genbox -cp presol_mod.gro -o em_sol.gro -p ${prot}.top -cs

if [ ! -f em_sol.gro ]; then
        echo "
        ######################
        Error in solvation, yo
        ######################"
        exit 0
fi

editconf -f em1.tpr -o ref1.pdb

mv ions_in.tpr ions_in_backup.tpr

grompp -f em.mdp -c em_sol.gro -p ${prot}.top -o ions_in

if [ ! -f ions_in.tpr ]; then
        echo "
        ########################
        Error in building genion
        ########################"
        exit 0
fi

echo SOL | genion -s ions_in.tpr -p ${prot}.top -o ions.gro -neutral -conc 0.150

sed -i '/SOL/d' ${prot}.top
sed -i '/NA/d' ${prot}.top
sed -i '/CL/d' ${prot}.top

tail -n+3 ions.gro > ions_mod.gro
atomnumbercount1=`grep -c [0-9] ions_mod.gro | awk '{print $1 }'`
atomnumbercount2=`expr $atomnumbercount1 - 1`
sed '1i '$atomnumbercount2''  ions_mod.gro -i
sed '1i Membrane protein wiv ions and everything' ions_mod.gro -i
#sed -i '$ d' all.gro
#tail -n 1 protmem.gro > bottom
#cat bottom >> all.gro

sed -i '/SOL/d' ${prot}.top
uniqatomssol=`grep SOL ions_mod.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep SOL ions_mod.gro | awk '{print $1}' | wc -l`
echo SOL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top

sed -i '/NA/d' ${prot}.top
uniqatomssol=`grep NA ions_mod.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep NA ions_mod.gro | awk '{print $1}' | wc -l`
echo NA `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top

sed -i '/CL/d' ${prot}.top
uniqatomssol=`grep CL ions_mod.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep CL ions_mod.gro | awk '{print $1}' | wc -l`
echo CL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> ${prot}.top

mv em2.tpr em2_backup.tpr

grompp -v -f em -o em2 -c ions_mod.gro -p ${prot}.top

if [ ! -f em2.tpr ]; then
        echo "
        #######################
        Error in building min 2
        #######################"
        exit 0
fi

mdrun -v -s em2 -o em2 -e em2 -c after_em2

editconf -f after_em2.gro -o after_em2.pdb

editconf -f em2.tpr -o ref2.pdb

#
## 10;; Prepare for equilibration
#

#grompp -f pr_POPC_berger_test.mdp -c after_em2.gro -p ${prot}.top -o pr_test.tpr
#grompp -f pr_POPC_berger.mdp -c after_em2.gro -p ${prot}.top -o pr_full.tpr
# pr_POPC_berger.mdp is fucked! No pressure!

grompp -f pr_GPU_POPC_test.mdp -c after_em2.gro -p ${prot}.top -o pr_test.tpr
grompp -f pr_GPU_POPC.mdp -c after_em2.gro -p ${prot}.top -o pr_full.tpr

mv pr_test.gro pr_test_backup.gro

mdrun -v -deffnm pr_test

if [ ! -f pr_test.gro ]; then
        echo "
        ####################################
        Your equilibration has exploded!!!!!
        ####################################"
        exit 0
elif [ -f pr_test.gro ]; then
        echo "
        ##########################################
        Equilibration seems fine. Run a 1 ns equil
        with: mdrun -v -deffnm pr_full >& out &
        ##########################################"
        exit 0
fi

mkdir Backup
mv *backup* Backup/

echo "

.....phase 3: profit!


"

## clean up

rm -f arg bottom chains his ndx pdboutput protonation asp chain1 glu lys ndx2 preprotonation_* prop1 512-popc-eq.pdb ${prot}_backup.gro ${prot}.gro ${prot}.gro.pdb ${prot}-mem_backup.pdb.gro ${prot}-mem.pdb.gro ${prot}_postpdb2gmx.gro all_backup.gro all_backup.pdb all.gro all.pdb CLx.pdb conf.gro cons_pull.mdp e1.edr em1_backup.tpr em1.mdp em1.tpr em2.edr em2.tpr em_gmx5.mdp em.mdp em_short.mdp emshort.mdp em_sol.gro eq_fold1.mdp eq_fold2.mdp eq_npt.mdp forpropka_*.pdb forpropka_*.propka_input gmembed_backup.tpr gmembed.mdp gmembed_nuc.mdp gmembed_NUC.mdp gmembed_nuc_sol.mdp gmembed_pop.mdp gmembed.tpr index.ndx ions.gro ions_in.tpr ions_mod.gro md_fold.mdp md_GPU_gmx5_pull.mdp md_GPU_gmx5_soluble.mdp md_GPU.mdp md_GPU_POPE.mdp md.log md_posre.mdp md_sol.mdp membed.dat membedded_backup.gro memsol.ndx memsolpre.ndx Memx.gro Memx.pdb MGx.pdb NAx.pdb NPT.mdp NVT.mdp premembed_backup.gro premembed_backup.pdb premembed.gro premembed.pdb presol.gro presol_mod.gro pr_GPU_POPC_gmx5.mdp pr_GPU_POPC.mdp pr_GPU_POPC_test_gmx5.mdp pr_GPU_POPC_test.mdp protein.gro protein.pdb protmem.gro protmem.gro_backup.pdb protmem.gro.pdb pr_test.cpt pr_test.edr pr_test.gro pr_test.log pr_test.tpr ref1.pdb ref2.pdb solmem_backup.gro solmem.gro SOLx.gro SOLx.pdb state.cpt state_prev.cpt topol.top

