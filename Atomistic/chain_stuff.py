# coding: utf-8
from Bio.PDB import *
parser = PDBparser()
parser = PDBParser()
structure = parser.get_structure('WT', 'WT.pdb')
io = PDBIO()
io.set_structure(structure)
chain = structure[0]['C']
io.save('test.pdb', chain)
model = structure[0]
chain = model['C']
io.save('test.pdb', chain)
print chain
io.save('test.pdb')
io.chain = model['C']
io.save('test.pdb')
io.save('test.pdb', chain2)
chain2 = model['C']
io.save('test.pdb', chain2)
class ChainSelect(Select):
    def accept_chain(self, chain):
        if chain.get_name() =='C':
            return 1
        else:
            return 0
        
io.save('test.pdb', ChainSelect() )
select.accept_chain(chain)
class ChainSelect(Select):
    def accept_chain(self, chain):
        if chain.get_id() =='C':
            return 1
        else:
            return 0
  
io.save('test.pdb', ChainSelect() )
get_ipython().magic(u'save chain_stuff')
get_ipython().magic(u'save chain_stuff 1-25')
