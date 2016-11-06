#
# Script to run a membrane protein from pdb state to membrane embedded state
# Need to make sure pdb file is numbered correctly and any chains are labelled
# Modified for SecA inclusion - specific for this system
# 

###### MODIFIED 24/11/2014 for new gmembed - including NUC in gmembed option###############

#
## Need to define prot and nuc options (ATP/ADP)
#

prot=AYEG
nuc=ATP
mem=POP

# get equilibrated POPC membrane and other setup files

cp /home/robinc/Desktop/MD_KIT/512-popc-eq.pdb .
cp /home/robinc/Desktop/MD_KIT/SOL.gro .   # Added to get around g_membed issues in 4.6.5
cp /home/robinc/Desktop/MD_KIT/*mdp .
cp /home/robinc/Desktop/MD_KIT/*gro .
cp /home/robinc/Desktop/MD_KIT/*itp .
cp /home/robinc/Desktop/MD_KIT/popc_opls.itp .
cp /home/robinc/Desktop/MD_KIT/membed.dat .

#
## Generate topology file from pdb
#

pdb2gmx -f $prot.pdb -o $prot.gro -p $prot.top -ff oplsaa -water spc -ignh    	# add protonation options if required
editconf -f $prot.gro -o $prot.gro.pdb
sed -e 's/TER/   /g' -e 's/ENDMDL/      /g' $prot.gro.pdb -i
cat $prot.gro.pdb 512-popc-eq.pdb > protmem.gro.pdb

editconf -f protmem.gro.pdb -o $prot-mem.pdb.gro -c -box 13.1708 13.1708 15 -bt triclinic
size=`tail -n 1 $prot-mem.pdb.gro`
grep -v "^$" $prot-mem.pdb.gro > protmem.gro

uniqatoms=`grep $mem $prot-mem.pdb.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatoms=`grep $mem $prot-mem.pdb.gro | awk '{print $1}' | wc -l`
echo $mem `echo $totalatoms $uniqatoms | awk '{print $1/$2}'` >> $prot.top

sed '/Include water topology/i \
# include "'popc_opls.itp'"\
' $prot.top -i 

# 
## Topology created
#

#
## Setting up gmembed - MODIFIED FOR MUC INCLUSION
#

# ADD SOL FOR GMEMBED

genbox -cp $prot-mem.pdb.gro -o solmem.gro -p $prot.top -cs
sed -i -e "1,2d" solmem.gro
cat $nuc.gro solmem.gro > nucsolmem.gro

sed -i '/; Compound        #mols/a \Other_chain_B 1\' $prot.top
sed -e '/; Include chain topologies/a \#include "'"$nuc.itp"'"' $prot.top > tmptp
cat tmptp > $prot.top


#tail -n 1 $prot-mem.pdb.gro > bottom
#cat bottom >> solmem.gro

atomnumbercount1=`grep -c [0-9] nucsolmem.gro | awk '{print $1 }'`
atomnumbercount2=`expr $atomnumbercount1 - 1`
sed '1i '$atomnumbercount2'' nucsolmem.gro -i
sed '1i Membrane protein with Nucleotide' nucsolmem.gro -i

sed -i '/SOL/d' $prot.top
uniqatomssol=`grep SOL nucsolmem.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatomssol=`grep SOL nucsolmem.gro | awk '{print $1}' | wc -l`
echo SOL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> $prot.top

cp /home/robinc/Desktop/MD_KIT/*ndx .
grompp -v -f gmembed_NUC.mdp -o gmembed -c nucsolmem.gro -p $prot.top -n gmem_$nuc.ndx
echo NUC $mem | mdrun -s gmembed.tpr -membed membed.dat -o traj.trr -c membeddedsol.gro -e ener.edr -nt 1 -cpt -1 -v -stepout 100 -mn gmem_$nuc
##############################
grep -Ev 'SOL' membeddedsol.pdb > membedded.gro
sed -i '/SOL/d' $prot.top

sed -i -e "1,2d" membedded.gro
atomnumbercount1=`grep -c [0-9] membedded.gro | awk '{print $1 }'`
atomnumbercount2=`expr $atomnumbercount1 - 1`
sed '1i '$atomnumbercount2''  membedded.gro -i
sed '1i Membrane protein' membedded.gro -i

#
## gmembed completed
#

#
## Minimizing and solvating system
#

grompp -v -f em.mdp -o em1 -e em1 -c membedded.gro  -p $prot.top
mdrun -v -s em1 -o em1 -e e1 -c after_em1.gro 
editconf -f em1.tpr -o ref1.pdb

genbox -cp after_em1.gro -o b4em -p $prot.top -cs
grompp -v -f em -o em -c b4em.gro -p $prot.top 
echo SOL | genion -s em.tpr -p $prot.top -o ions.gro -neutral -conc 0.150
grompp -v -f em -o em2 -c ions.gro -p $prot.top
mdrun -v -s em2 -o em2 -e em2 -c after_em2
editconf -f em2.tpr -o ref2.pdb

#
## Preparing final .gro file
#

editconf -f after_em2.gro -o after_em2.pdb
grep POP after_em2.pdb> POPCx.pdb
grep SOL after_em2.pdb> SOLx.pdb
grep NA after_em2.pdb> NAx.pdb
grep CL after_em2.pdb> CLx.pdb    
grep $nuc after_em2.pdb> Nucx.pdb   
grep MG after_em2.pdb> MGx.pdb   

editconf -f after_em2.gro -o after_em2.pdb

grep -e ALA -e CYS -e ASP -e GLU -e PHE -e GLY -e HIS -e ILE -e LYS -e LEU -e MET -e ASN -e PRO -e GLN -e ARG -e SER -e THR -e VAL -e TRP -e TYR  after_em2.pdb| grep ATOM > protein.pdb

dowser protein.pdb

editconf -f protein.pdb -o protein.gro 
editconf -f POPCx.pdb -o POPCx.gro
editconf -f SOLx.pdb -o SOLx.gro
editconf -f NAx.pdb -o NAx.gro
editconf -f CLx.pdb -o CLx.gro
editconf -f Nucx.pdb -o Nucx.gro
editconf -f MGx.pdb -o MGx.gro
sed 's/HOH/SOL/g' dowserwat.gro -i

rm -f all.gro

cat MGx.gro Nucx.gro protein.gro POPCx.gro SOLx.gro dowserwat.gro NAx.gro  CLx.gro | grep -e MG -e ATP -e ADP  -e ALA -e CYS -e ASP -e GLU -e PHE -e GLY -e HIS -e ILE -e LYS -e LEU -e MET -e ASN -e PRO -e GLN -e ARG -e SER -e THR -e VAL -e TRP -e TYR -e POP -e SOL -e POG -e NA -e CL -e DDM >> all.gro

waternumber=`grep -c OW all.gro`
sed 's/SOL.*/SOL '$waternumber'/g' $prot.top -i

atomnumbercount=`grep -c [0-9] all.gro | awk '{print $1 }'`
sed '1i '$atomnumbercount'' all.gro -i
sed '1i Membrane protein' all.gro -i
tail -n 1 after_em2.gro > bottom2
cat bottom2 >> all.gro

grompp -v -f em -o em3 -c all.gro -p $prot.top
mdrun -v -s em3 -o em3 -e em3 -c after_em3



