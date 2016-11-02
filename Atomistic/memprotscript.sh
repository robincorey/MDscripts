#
# Script to run a membrane protein from pdb state to membrane embedded state
# Need to make sure pdb file is numbered correctly and any chains are labelled
# Modified for SecA inclusion - specific for this system
# 


#
## Need to define prot and nuc options (ATP/ADP)
#

prot=ATPpos
nuc=TPp
mem=POP

# get equilibrated POPC membrane and other setup files

cp /home/robinc/Desktop/MD_KIT/512-popc-eq.pdb .
cp /home/robinc/Desktop/MD_KIT/SOL.gro .   # Added to get around g_membed issues in 4.6.5
cp /home/robinc/Desktop/MD_KIT/*mdp .
cp /home/robinc/Desktop/MD_KIT/*gro .
cp /home/robinc/Desktop/MD_KIT/membed.dat .
cp /home/robinc/Desktop/MD_KIT/*itp .
cp /home/robinc/Desktop/MD_KIT/popc_opls.itp .

#
## Generate topology file from pdb
#

pdb2gmx -f $prot.pdb -o $prot.gro -p $prot.top -ff oplsaa -water spc -ignh    	# add protonation options if required
editconf -f $prot.gro -o $prot.gro.pdb
sed -e 's/TER/   /g' -e 's/ENDMDL/      /g' $prot.gro.pdb -i
cat $prot.gro.pdb 512-popc-eq.pdb > protmem.gro.pdb

grep -Ev 'REMARK|MODEL|TER|ENDMDL|TITLE|CRYST1' protmem.gro.pdb > protmem_out.gro.pdb
grep -v "^$" protmem_out.gro.pdb > protmem_out_2.gro.pdb

# 
## Topology created
#

editconf -f protmem_out_2.gro.pdb -o $prot-mem.pdb.gro -c -box 13.1708 13.1708 15 -bt triclinic
#grep -e POP $prot-mem.pdb.gro | awk '{print $1}' |  sed 's/[0-9]//g' | uniq > species

#for list in `cat species`
#do
#searchterm=`echo $list | sed 's/-/ /g' | awk '{print $1}'`
#uniqatoms=`grep $searchterm $prot-mem.pdb.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
#totalatoms=`grep $searchterm $prot-mem.pdb.gro | awk '{print $1}' | wc -l`
#echo $searchterm  `echo $totalatoms $uniqatoms | awk '{print $1/$2}'` >> $prot.top
#done

uniqatoms=`grep $mem $prot-mem.pdb.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatoms=`grep $mem $prot-mem.pdb.gro | awk '{print $1}' | wc -l`
echo $mem  `echo $totalatoms $uniqatoms | awk '{print $1/$2}'` >> $prot.top


# adding itp file to top file

sed '/Include water topology/i \
#include "'popc_opls.itp'"\
' $prot.top -i 

#
## Setting up gmembed - assumes gmembed.mdp is already completed. Use section from John's script if not
#

######## ADD SOL FOR GMEMBED

genbox -cp $prot-mem.pdb.gro -o solmem.gro -p $prot.top -cs


#sed '$d' < $prot-mem.pdb.gro > prot.gro
#cat prot.gro SOL2.gro > solmem.gro
#tail -n 1 $prot-mem.pdb.gro > bottom
#cat bottom >> solmem.gro

#sed -i -e "1,2d" solmem.gro
#atomnumbercount1=`grep -c [0-9] solmem.gro | awk '{print $1 }'`
#atomnumbercount2=`expr $atomnumbercount1 - 1`
#sed '1i '$atomnumbercount2'' solmem.gro -i
#sed '1i Membrane protein' solmem.gro -i

#sed -i '/SOL/d' $prot.top

#uniqatomssol=`grep SOL solmem.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
#totalatomssol=`grep SOL solmem.gro | awk '{print $1}' | wc -l`
#echo SOL `echo $totalatomssol $uniqatomssol | awk '{print $1/$2}'` >> $prot.top

grompp -v -f gmembed.mdp -o gmembed -c solmem.gro -p $prot.top
echo 1 POP | mdrun -s gmembed.tpr -membed membed.dat -o traj.trr -c membeddedsol.gro -e ener.edr -nt 1 -cpt -1 -v -stepout 100

grep -Ev 'SOL' membeddedsol.gro > membedded.gro

sed -i '/SOL/d' $prot.top

#sed -i -e "1,2d" membedded.gro
#cat $nuc.gro membedded.gro > $nuc-mem.gro
#tail -n 1 $prot-mem.pdb.gro > bottom
#atomnumbercount1=`grep -c [0-9] $nuc-mem.gro | awk '{print $1 }'`
#atomnumbercount2=`expr $atomnumbercount1 - 1`
#sed '1i '$atomnumbercount2''  $nuc-mem.gro -i
#sed '1i Membrane protein' $nuc-mem.gro -i
#cat bottom >> $nuc-mem.gro
#######  EXTRA SOL REMOVED


#
## gmembed completed
#



#
## Adding nucleotide to the protein - specify at beginning of script
#

sed -i -e "1,2d" membedded.gro
cat $nuc.gro membedded.gro > $nuc-mem.gro
tail -n 1 $prot-mem.pdb.gro > bottom
atomnumbercount1=`grep -c [0-9] $nuc-mem.gro | awk '{print $1 }'`
atomnumbercount2=`expr $atomnumbercount1 - 1`
sed '1i '$atomnumbercount2''  $nuc-mem.gro -i
sed '1i Membrane protein' $nuc-mem.gro -i
cat bottom >> $nuc-mem.gro

sed -i '/; Compound        #mols/a \Other_chain_B 1\' $prot.top
sed -e '/; Include chain topologies/a \#include "'"$nuc.itp"'"' $prot.top > tmptp
cat tmptp > $prot.top

sed -i '/POP/d' $prot.top

uniqatoms=`grep $mem $nuc-mem.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
totalatoms=`grep $mem $nuc-mem.gro | awk '{print $1}' | wc -l`
echo $mem  `echo $totalatoms $uniqatoms | awk '{print $1/$2}'` >> $prot.top


#for list in `cat species`
#do 
#searchterm=`echo $list | sed 's/-/ /g' | awk '{print $1}'`
#uniqatoms=`grep $searchterm $nuc-mem.gro | awk '{print $1}' | uniq -c | tail -n 1 | awk '{print $1}'`
#totalatoms=`grep $searchterm $nuc-mem.gro | awk '{print $1}' | wc -l`
#echo $searchterm `echo $totalatoms $uniqatoms | awk '{print $1/$2}'` >> $prot.top 
#done

#
## Minimizing and solvating system
#

grompp -v -f em.mdp -o em1 -e em1 -c $nuc-mem.gro  -p $prot.top
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


