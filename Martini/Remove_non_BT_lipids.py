#!/sansom/s137/bioc1535/anaconda2/bin/

import MDAnalysis
from MDAnalysis.analysis.leaflet import LeafletFinder
import MDAnalysis.analysis.leaflet

u = MDAnalysis.Universe("md.pdb","md-center.xtc")

for ts in u.trajectory[4999:5000]:
	L = MDAnalysis.analysis.leaflet.LeafletFinder(u, '(name C4*) and resname POP*',cutoff=20,pbc=True)

bilayer = L.groups(0).residues.atoms 
extras = L.groups(1).residues.atoms 

o = u.select_atoms("protein or resname CYS*") + bilayer + u.select_atoms("resname W or resname CL or resname NA or resname ION")

o.write("md-cut.pdb")

extras.write("extra.pdb")

with MDAnalysis.Writer("md-cut.xtc", o.n_atoms) as W:
    for ts in u.trajectory:
`        W.write(o)
