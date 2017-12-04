# coding: utf-8
io = PDBIO()
from Bio.PDB import *
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
    
io = PDBIO()
structure = parser.get_structure('WT', 'WT.pdb')
parser = PDBParser()
structure = parser.get_structure('WT', 'WT.pdb')
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
    
seq = polypeptide.get_sequence()
seq = pp.get_sequence()
print seq
structure = parser.get_structure('WT', 'WT_C.pdb')
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
    
get_ipython().magic(u'save fasta_stuff 1-15')
