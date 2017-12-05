#!/home/birac/anaconda2/bin/python
import sys
import os
from Bio.PDB import *
import fileinput

GRO = 'pr_full.gro'
PDB_chain = 'pr_test_chainA.pdb'
RES = [ 'E643', 'D646', 'E789', 'E793' ]

def prep_pdb(input_pdb):
    os.system("convpdb.pl %s > pre_mutate.pdb" % input_pdb)
    for line in fileinput.input("pre_mutate.pdb", inplace=True):
        if "HSE" in line:
            newhse = line.replace("HSE", "HIS")
            print newhse,
            continue
        if "HSD" in line:
            newhsd = line.replace("HSD", "HIS")
            print newhse,
            continue
        print(line),
    for line in fileinput.input("pre_mutate.pdb", inplace=True):
        if 'O2' in line:
            continue
        if 'O1' in line:
            newter = line.replace("O1", "O ")
            print(newter),
            continue
        if 'OT1' in line:
            newter = line.replace("OT1", "O  ")
            print(newter),
            continue
        if 'OT2' in line:
            continue
        print(line),

def make_mutation(input_pdb, number, residue, mutation, prepdb):
    io = PDBIO()
    ppb = PPBuilder()
    parser = PDBParser()
    structure = parser.get_structure('chain', input_pdb)
    with open('%s_%s_%s.seq'  % (residue, number, mutation), 'w') as f:
        for pp in ppb.build_peptides(structure):
            seq = str(pp.get_sequence())
            previous = int(number) - 1
            nextup = int(number) + 0
            mutated_sequence = seq[:previous] + mutation + seq[nextup:]
            f.write(mutated_sequence)
    os.system("Scwrl4 -i %s -s %s_%s_%s.seq -o %s_%s_%s.pdb > %s_%s_%s.txt" % (prepdb, residue, number, mutation, residue, number, mutation, residue, number, mutation))


def remake_itp_files(residue, number, mutation):
    name = "%s_%s_%s" % (residue, number, mutation)
    os.system("$GMX51 pdb2gmx -f %s.pdb -p %s -i %s_posre -o %s.gro -ignh -ff oplsaa -water spc" % (name, name, name, name))
    delete_these = ['forcefield']
    with open('%s.top' % name) as oldfile, open('%s.itp' % name, 'w') as newfile:
        for line in oldfile:
            if not any(delete_this in line for delete_this in delete_these):
                newfile.write(line)
    fin = fileinput.input('%s.itp' % name, inplace=1)
    with open('%s_2.itp' % name, 'w') as newfile2:
        for line in fin:
            newfile2.write(line)
            if '#endif' in line:
               break
    fin.close()

def update_topology(residue, number, mutation, input_topology, input_itp):
    name = "%s_%s_%s" % (residue, number, mutation)
    with open(input_topology) as oldfile, open('%s_remade.top' % name, 'w') as newfile2:
        for line in oldfile:
            if input_itp in line:
                newline = line.replace(input_itp, '%s_2.itp' % name)
                newfile2.write(newline)
            elif 'Protein_chain_A' in line:
                newline = line.replace('Protein_chain_A' , 'Protein')
                newfile2.write(newline)
            else:
                newfile2.write(line)

def minimize_system(residue, number, mutation):
    name = "%s_%s_%s" % (residue, number, mutation)
    os.system("$GMX51 grompp -f /home/birac/Desktop/MD_KIT/em_gmx5.mdp -c %s.gro -p %s.top -o %s_em.tpr" % (name, name, name))
    os.system("$GMX51 mdrun -v -deffnm %s_em" % name )

#prep_pdb(PDB_chain)

for res in RES:
    residue_type = ''.join([i for i in res if not i.isdigit()])
    residue_number = ''.join([i for i in res if i.isdigit()])
    if residue_type == "E":
        mutate_residue = "Q"
    elif residue_type == "D":
        mutate_residue = "N"
    mutate_residue2 = "A"
    print("changing %s%s to %s / %s" % (residue_type, residue_number, mutate_residue, mutate_residue2))
    #make_mutation('pre_mutate.pdb', residue_number, residue_type, mutate_residue, 'pre_mutate.pdb')
    #make_mutation(PDB_chain, residue_number, residue_type, mutate_residue2, 'pre_mutate.pdb')
    #remake_itp_files(residue_type, residue_number, mutate_residue2)
    #update_topology(residue_type, residue_number, mutate_residue2, 'AYEG_WT.top', 'AYEG_WT_Protein_chain_A.itp')
    minimize_system(residue_type, residue_number, mutate_residue2)
