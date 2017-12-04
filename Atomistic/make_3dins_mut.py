#!/home/birac/anaconda2/bin/python
from modeller import *
from modeller.automodel import *
import re
from Bio.PDB import *
import textwrap 

# Define these
array = ['WT', 'E643Q', 'D646N', 'BOTH']
chains = ['C']
residues_to_model = ['42','61']

def get_chains(chain_id):
  parser = PDBParser()
  structure = parser.get_structure(protein, '%s.pdb' % protein)
  class ChainSelect(Select):
    def accept_chain(self, chain):
      if chain.get_id() == '%s' % chain_id:
          return 1
      else:
          return 0
  io = PDBIO()
  io.set_structure(structure)
  io.save('%s_%s.pdb' % (protein, chain_id), ChainSelect() )

def make_sequence_align(chain_id):
  io = PDBIO()
  parser = PDBParser()
  ppb = PPBuilder()
  structure = parser.get_structure(protein, '%s_%s.pdb' % (protein, chain_id))
  f = open('%s_%s.ali' % (protein, chain_id), 'w')
  f.write('>P1;%s\n' % protein)
  f.write('sequence:%s::::::::\n' % protein)
  for pp in ppb.build_peptides(structure):
    seq = str(pp.get_sequence()) 
    f.write('%s' % seq)
  f.write('*')
  env = environ()
  aln = alignment(env)
  mdl = model(env, file='3din_%s.pdb' % chain_id) #, model_segment=('FIRST:A','LAST:A'))
  aln.append_model(mdl, align_codes='3din', atom_files='3din_%s.pdb' % chain_id)
  aln.append(file='%s_%c.ali' % (protein, chain_id), align_codes='%s' % protein)
  aln.align2d()
  aln.write(file='%s-3din.ali' % protein, alignment_format='PIR')

def loop(prot, start, end):
  log.verbose()
  env = environ()
  class MyLoop(loopmodel):
    def select_loop_atoms(self):
        startint = int(start)
        endint = int(end)
        return selection(self.residue_range(startint, endint))
  m = MyLoop(env, inimodel='%s.pdb' % prot, sequence=prot)
  m.loop.starting_model= 1
  m.loop.ending_model  = 10
  m.loop.md_level = refine.very_fast
  m.make()
   
class MyModel(automodel):
    def select_atoms(self):
        return selection(self.residue_range('246', '249'))
'''
a = MyModel(env, alnfile = 'mafft.ali',
            knowns = 'Tom_model', sequence = 'Tom_loops')
a.starting_model= 1
a.ending_model  = 1
'''
 
for protein in array:
  for chain_id in chains:
    get_chains(chain_id)
    make_sequence_align(chain_id)
    #loop('%s_%s' % (protein, chain_id), residues_to_model[0], resdiues_to_model[1])
