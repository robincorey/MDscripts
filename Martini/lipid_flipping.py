#!/sansom/s137/bioc1535/anaconda2/bin/python
#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path


#=========================================================================================
# create parser
#=========================================================================================
version_nb="0.0.7"
parser = argparse.ArgumentParser(prog='ff_contacts', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
[ USAGE ]
	
Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		1	: process every t-frames
Lipids identification  
-----------------------------------------------------
--flipflops		: input file with flipflopping lipids, see note 2
--beads			: leaflet identification technique, see note 3(a)
--leaflets	optimise: leaflet identification technique, see note 3(b)
--use_gro		: use gro file instead of xtc, see note 3(b)
Protein clusters identification and contacts
-----------------------------------------------------
--proteins		: protein selection file, (optional, see note 6)
--groups		: cluster groups definition file, see note 4(c)
--algorithm	min	: 'cog','min' or 'density', see 'DESCRIPTION'
--nx_cutoff 	6	: networkX cutoff distance for protein-protein contact (Angstrom)
--db_radius 	20	: DBSCAN search radius (Angstrom)
--db_neighbours	3	: DBSCAN minimum number of neighbours within a circle of radius --db_radius	
 Contacts profile options
-----------------------------------------------------
--profile 		: calculate local bilayer normal and store position of each contacts
--normal	svd	: local normal to bilayer ('z', 'cog' or 'svd'), see note 5
--normal_d	50	: distance of points to take into account for local normal, see note 5
--pl_cutoff 	6	: cutoff distance for protein-lipid contact (Angstrom)
--slices_dist	40 	: max distance from center of bilayer (Angstrom)
--slices_thick	0.5 	: thickness of the profile slices (Angstrom)
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[1], type=int, help=argparse.SUPPRESS)

#lipids identification options
parser.add_argument('--beads', nargs=1, dest='beadsfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)
parser.add_argument('--use_gro', dest='use_gro', action='store_true', help=argparse.SUPPRESS)

#protein options
parser.add_argument('--proteins', nargs=1, dest='selection_file_prot', default=['auto'], help=argparse.SUPPRESS)
parser.add_argument('--groups', nargs=1, dest='cluster_groups_file', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--algorithm', dest='m_algorithm', choices=['min','cog','density'], default='min', help=argparse.SUPPRESS)
parser.add_argument('--nx_cutoff', nargs=1, dest='nx_cutoff', default=[6], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_radius', nargs=1, dest='dbscan_dist', default=[20], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_neighbours', nargs=1, dest='dbscan_nb', default=[3], type=int, help=argparse.SUPPRESS)

#profile options
parser.add_argument('--profile', dest='profile', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--normal', dest='normal', choices=['z','cog','svd'], default='svd', help=argparse.SUPPRESS)
parser.add_argument('--normal_d', nargs=1, dest='normal_d', default=[50], type=float, help=argparse.SUPPRESS)
parser.add_argument('--pl_cutoff', nargs=1, dest='cutoff_pl', default=[6], type=float, help=argparse.SUPPRESS)
parser.add_argument('--slices_dist', nargs=1, dest='slices_dist', default=[40], type=float, help=argparse.SUPPRESS)
parser.add_argument('--slices_thick', nargs=1, dest='slices_thick', default=[0.5], type=float, help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

#parse user inputs
#-----------------
args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start = args.t_start[0]
args.t_end = args.t_end[0]
args.frames_dt = args.frames_dt[0]
#lipids identification options
args.beadsfilename = args.beadsfilename[0]
args.cutoff_leaflet = args.cutoff_leaflet[0]
args.selection_file_ff = args.selection_file_ff[0]
#protein options
args.selection_file_prot = args.selection_file_prot[0]
args.cluster_groups_file = args.cluster_groups_file[0]
args.nx_cutoff = args.nx_cutoff[0]
args.dbscan_dist = args.dbscan_dist[0]
args.dbscan_nb = args.dbscan_nb[0]
#profile options
args.normal_d = args.normal_d[0]
args.cutoff_pl = args.cutoff_pl[0]
args.slices_dist = args.slices_dist[0]
args.slices_thick = args.slices_thick[0]

#process options
#---------------
global lipids_ff_nb
lipids_ff_nb = 0
	
#leaflet identification
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should be a number or 'large', see note 2"
		sys.exit(1)

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy as sp
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================

if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.selection_file_ff != "no" and not os.path.isfile(args.selection_file_ff):
	print "Error: file " + str(args.selection_file_ff) + " not found."
	sys.exit(1)
if args.selection_file_prot != "auto" and not os.path.isfile(args.selection_file_prot):
	print "Error: file " + str(args.selection_file_prot) + " not found."
	sys.exit(1)
if args.cluster_groups_file != "no" and not os.path.isfile(args.cluster_groups_file):
	print "Error: file " + str(args.cluster_groups_file) + " not found."
	sys.exit(1)
if args.beadsfilename != "no" and not os.path.isfile(args.beadsfilename):
	print "Error: file " + str(args.beadsfilename) + " not found."
	sys.exit(1)
if args.t_end != -1 and args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)
if args.normal != 'z' and args.normal_d <= 0:
	print "Error: --normal_d (" + str(args.normal__d) + " AA) should be greater 0. Or choose 'z' for --normal, see note 5."
	sys.exit(1)

if args.xtcfilename == "no":
	if args.use_gro:
		print "Error: you can only use the --use_gro file if you specify an xtc file."
		sys.exit(1)
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	args.output_folder = "ff_ctcts_" + args.xtcfilename[:-4]

if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	os.mkdir(args.output_folder)
	os.mkdir(args.output_folder + '/by_size')
	os.mkdir(args.output_folder + '/by_type')
	if args.cluster_groups_file != "no":
		os.mkdir(args.output_folder + '/by_group')
	if args.profile:
		os.mkdir(args.output_folder + '/TM_profile')
		os.mkdir(args.output_folder + '/TM_profile/xvg')
		os.mkdir(args.output_folder + '/TM_profile/png')
		if args.cluster_groups_file != "no":
			os.mkdir(args.output_folder + '/TM_profile/groups')
			os.mkdir(args.output_folder + '/TM_profile/groups/xvg')
			os.mkdir(args.output_folder + '/TM_profile/groups/png')
	
	#create log
	#----------
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/ff_contacts.log'
	output_log=open(filename_log, 'w')		
	output_log.write("[ff_contacts v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python ff_contacts.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()

	#copy input files
	#----------------
	if args.selection_file_ff != "no":
		shutil.copy2(args.selection_file_ff,args.output_folder + "/")	
	if args.selection_file_prot != "no" and args.selection_file_prot != "auto":
		shutil.copy2(args.selection_file_prot,args.output_folder + "/")
	if args.cluster_groups_file != "no":
		shutil.copy2(args.cluster_groups_file,args.output_folder + "/")
	if args.beadsfilename != "no":
		shutil.copy2(args.beadsfilename,args.output_folder + "/")

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def set_lipids_beads():

	global leaflet_sele_string

	#set default beads
	leaflet_beads = {}
	leaflet_beads['martini'] = "name PO4 or name PO3 or name B1A"
	leaflet_sele_string = leaflet_beads['martini']

	#use users input
	if args.beadsfilename != "no":
		with open(args.beadsfilename) as f:
			lines = f.readlines()
		if len(lines) > 1:
			print "Error: the file " + str(args.beadsfilename) + " should conly ontain 1 line (" + str(len(lines)) + " found), see note 2(a)."
			sys.exit(1)
		else:
			if lines[0][-1] == "\n":
				lines[0] = lines[0][:-1]
			leaflet_sele_string = lines[0]

	return
def load_MDA_universe():
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global frames_to_write
	global nb_frames_to_process
	global f_start
	global f_end
	global radial_bins
	global radial_bin_max
	global radial_radius_max
	f_start = 0

	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		frames_to_write = [True]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		if args.use_gro:
			global U_gro
			U_gro = Universe(args.grofilename)
		U = Universe(args.grofilename, args.xtcfilename)
		U_timestep = U.trajectory.dt
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes		

		#sanity check
		if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
			print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
			sys.exit(1)
		if U.trajectory.numframes < args.frames_dt:
			print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

		#rewind traj (very important to make sure that later the 1st frame of the xtc will be used for leaflet identification)
		U.trajectory.rewind()

		#create list of index of frames to process
		if args.t_end != -1:
			f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_end < 0:
				print "Error: the starting time specified is before the beginning of the xtc."
				sys.exit(1)
		else:
			f_end = nb_frames_xtc - 1		
		if args.t_start != -1:
			f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_start > f_end:
				print "Error: the starting time specified is after the end of the xtc."
				sys.exit(1)
		if (f_end - f_start)%args.frames_dt == 0:
			tmp_offset = 0
		else:
			tmp_offset = 1		
		frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
		nb_frames_to_process = len(frames_to_process)

	#check for the presence of proteins
	test_prot = U.selectAtoms("protein")
	if test_prot.numberOfAtoms() == 0:
			print "Error: no protein found in the system."
			sys.exit(1)
			
	#check for the presence of lipids
	if args.cutoff_leaflet != "no":
		test_beads = U.selectAtoms(leaflet_sele_string)
		if test_beads.numberOfAtoms() == 0:
			print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
			print "-> no lipid particles selected, check the --beads option."
			sys.exit(1)

	return
def identify_ff():
	print "\nReading selection file for flipflopping lipids..."
	
	#declare variables
	global lipids_ff_nb
	global lipids_ff_info
	global lipids_ff_resnames
	global lipids_ff_leaflet
	global lipids_ff_u2l_index
	global lipids_ff_l2u_index
	global lipids_sele_ff
	global lipids_sele_ff_bead
	global lipids_sele_ff_bonds
	global lipids_sele_ff_VMD_string
	global leaflet_sele_string
	lipids_ff_nb = 0
	lipids_ff_info = {}
	lipids_ff_resnames = []
	lipids_ff_leaflet = []
	lipids_ff_u2l_index = []
	lipids_ff_l2u_index = []
	lipids_sele_ff = {}
	lipids_sele_ff_bead = {}
	lipids_sele_ff_bonds = {}
	lipids_sele_ff_VMD_string={}
		
	with open(args.selection_file_ff) as f:
		lines = f.readlines()
	lipids_ff_nb = len(lines)
	print " -found " + str(lipids_ff_nb) + " flipflopping lipids"
	leaflet_sele_string = leaflet_sele_string + " and not ("
	for l_index in range(0,lipids_ff_nb):
		line = lines[l_index]
		if line[-1] == "\n":
			line = line[:-1]
		try:
			line_content = line.split(',')
			if len(line_content) != 6:
				print "Error: wrong format for line " + str(l_index+1) + " in " + str(args.selection_file_ff) + ", see note 4 in bilayer_perturbations --help."
				print " ->", line
				sys.exit(1)
			#read current lipid details
			lip_resname = line_content[0]
			lip_resnum = int(line_content[1])
			lip_leaflet = line_content[2]
			lip_bead = line_content[3]
			lip_tstart = float(line_content[4])
			lip_tend = float(line_content[5])
			if lip_tend == 0:
				lip_tend = U.trajectory.totaltime/float(1000)			#this is to handle lipids which haven't finished flip-flopping
			lipids_ff_info[l_index] = [lip_resname,lip_resnum,lip_leaflet,lip_bead,lip_tstart,lip_tend]
						
			#update: starting leaflets
			if lip_leaflet not in lipids_ff_leaflet:
				lipids_ff_leaflet.append(lip_leaflet)

			#update: index in directional lists
			if lip_leaflet == "upper":
				lipids_ff_u2l_index.append(l_index)
			elif lip_leaflet == "lower":
				lipids_ff_l2u_index.append(l_index)
			else:
				print "->unknown starting leaflet '" + str(lip_leaflet) + "'."
				sys.exit(1)
			
			#update: resnames
			if lip_resname not in lipids_ff_resnames:
				lipids_ff_resnames.append(lip_resname)
	
			#update: leaflet selection string
			if l_index==0:
				leaflet_sele_string+="(resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"
			else:
				leaflet_sele_string+=" or (resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"

			#create selections
			lipids_sele_ff[l_index] = U.selectAtoms("resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + " and name " + str(lipids_ff_info[l_index][3]))
			lipids_sele_ff_bead[l_index] = lipids_sele_ff[l_index].selectAtoms("name " + str(lip_bead))
			lipids_sele_ff_VMD_string[l_index]="resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1])
			if lipids_sele_ff[l_index].numberOfAtoms() == 0:
				print "Error:"
				print line
				print "-> no such lipid found."
				sys.exit(1)	
		except:
			print "Error: invalid flipflopping lipid selection string on line " + str(l_index+1) + ": '" + line + "'"
			sys.exit(1)
	leaflet_sele_string+=")"		

	return
def identify_proteins():	
	print "\nIdentifying proteins..."
	
	#import modules
	if args.m_algorithm == "density":
		global DBSCAN
		from sklearn.cluster import DBSCAN
	else:
		global nx
		import networkx as nx

	#declare variables
	global proteins_nb
	global proteins_sele
	global proteins_sele_string
	global proteins_sele_string_VMD
	global proteins_boundaries
	global proteins_nb_atoms
	global nb_atom_per_protein	
	global proteins_max_at_number	
	proteins_nb = 0
	proteins_sele = {}
	proteins_sele_string = {}
	proteins_sele_string_VMD = {}
	proteins_boundaries = {}
	
	#check for protein presence
	if U.selectAtoms("protein").numberOfAtoms() == 0:
		print "Error: no protein detected."
		sys.exit(1)
	
	#case: selection file provided
	if args.selection_file_prot != "auto":
		print " -reading protein selection file..."
		with open(args.selection_file_prot) as f:
			lines = f.readlines()
		proteins_nb=len(lines)
		proteins_sele["all"] = MDAnalysis.core.AtomGroup.AtomGroup([])
		for p_index in range(0,proteins_nb):
			line = lines[p_index]
			if line[-1] == "\n":
				line = line[:-1]
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			try:
				print " p[" + str(p_index) + "]=U.selectAtoms(" + line + ")"
				proteins_sele[p_index] = U.selectAtoms(line[1:-2])
				proteins_sele["all"] += proteins_sele[p_index]
				proteins_boundaries[p_index] = [proteins_sele[p_index].indices()[0] + 1, proteins_sele[p_index].indices()[proteins_sele[p_index].numberOfAtoms()]+1]
				proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
				proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			except:
				print "Error:"
				print line
				print "->invalid selection string."
				sys.exit(1)
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
	
	#case: automatic detection
	else:
		#declare local variables
		proteins_ca_nb = {}
		proteins_ca_nmax = 0
		proteins_ca_group = {}
		proteins_boundaries = {}
	
		#retrieve 1st atom info
		proteins_sele["all"] = U.selectAtoms("protein")
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
		prec_resnum = proteins_sele["all"][0].resnum
		prec_segid = proteins_sele["all"][0].segid
		prec_atnum = proteins_sele["all"][0].number+1
		prev_atnum = proteins_sele["all"][0].number+1						#atom corresponding to the beginning of the current protein
		#browse following atoms
		for a in proteins_sele["all"][1:]:
			delta_res = a.resnum-prec_resnum
			delta_atm = a.number+1-prec_atnum
			if delta_res < 0 or a.segid != prec_segid or delta_atm > 1:
				proteins_boundaries[proteins_nb] = [prev_atnum,prec_atnum]
				proteins_nb += 1
				prev_atnum = a.number + 1
			prec_resnum = a.resnum
			prec_atnum = a.number + 1
			prec_segid = a.segid		
		#add last protein section
		if prev_atnum < proteins_sele["all"][proteins_nb_atoms-1].number:
			proteins_boundaries[proteins_nb] = [prev_atnum,proteins_sele["all"][proteins_nb_atoms-1].number+1]
			proteins_nb += 1
		
		#display results
		print " -protein found:", proteins_nb
		print " -protein boundaries (atom numbers): see protein.sele file"
		#create protein selections and save into a txt file
		filename_sele=os.getcwd() + '/' + str(args.output_folder) + '/proteins.sele'
		output_stat = open(filename_sele, 'w')	
		output_stat.write("#This file was generated by the script bilayer_perturbations v" + str(version_nb) +"\n")
		output_stat.write("#The lines below correspond to MDAnalysis section string, e.g. U.selectAtoms(LINE)\n")
		output_stat.write("\n")	
		for p_index in range(0, proteins_nb):
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
			proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			proteins_sele[p_index] = U.selectAtoms(proteins_sele_string[p_index])
			output_stat.write(proteins_sele_string[p_index] + "\n")
		output_stat.close()

	#get protein details
	nb_atom_per_protein = proteins_sele[0].numberOfAtoms()
	proteins_max_at_number = np.max(proteins_sele["all"].indices())
	print ""

	return
def identify_leaflets():
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			if args.use_gro:
				cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U_gro, leaflet_sele_string)
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U_gro, leaflet_sele_string, cutoff_value[0])
			else:
				cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])			
		else:
			if args.use_gro:
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U_gro, leaflet_sele_string, args.cutoff_leaflet)
			else:
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if np.shape(L.groups())[0] < 2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			if args.use_gro:
				tmp_up = L.group(0)
				tmp_lw = L.group(1)
			else:
				leaflet_sele["upper"]["all species"] = L.group(0)
				leaflet_sele["lower"]["all species"] = L.group(1)
		else:
			if args.use_gro:
				tmp_up = L.group(1)
				tmp_lw = L.group(0)			
			else:
				leaflet_sele["upper"]["all species"] = L.group(1)
				leaflet_sele["lower"]["all species"] = L.group(0)
		
		if args.use_gro:
			tmp_up_indices = tmp_up.indices()
			tmp_lw_indices = tmp_lw.indices()
			tmp_up_nb_atoms = len(tmp_up_indices)
			tmp_lw_nb_atoms = len(tmp_lw_indices)
			leaflet_sele["upper"]["all species"] = U.selectAtoms("bynum " + str(tmp_up_indices[0] + 1))
			leaflet_sele["lower"]["all species"] = U.selectAtoms("bynum " + str(tmp_lw_indices[0] + 1))
			for index in range(1,tmp_up_nb_atoms):
				leaflet_sele["upper"]["all species"] += U.selectAtoms("bynum " + str(tmp_up_indices[index] + 1))
				progress = '\r -identifying upper leaflet from gro file... ' + str(round(index/float(tmp_up_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			print ''
			for index in range(1,tmp_lw_nb_atoms):
				leaflet_sele["lower"]["all species"] += U.selectAtoms("bynum " + str(tmp_lw_indices[index] + 1))
				progress = '\r -identifying lower leaflet from gro file... ' + str(round(index/float(tmp_lw_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			print ''
		
		leaflet_sele["both"]["all species"] = leaflet_sele["lower"]["all species"] + leaflet_sele["upper"]["all species"]
		if np.shape(L.groups())[0] == 2:
			print " -found 2 leaflets: ", leaflet_sele["upper"]["all species"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"]["all species"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, np.shape(L.groups())[0]):
				other_lipids += L.group(g).numberOfResidues()
			print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"]["all species"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"]["all species"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	#use cog and z coordinates in the GRO file supplied:
	else:
		if args.use_gro:
			tmp_all = U_gro.selectAtoms(leaflet_sele_string)
			tmp_lipids_avg_z = tmp_all.centerOfGeometry()[2]
			tmp_up = tmp_all.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
			tmp_lw = tmp_all.selectAtoms("prop z < " + str(tmp_lipids_avg_z))
			tmp_up_indices = tmp_up.indices()
			tmp_lw_indices = tmp_lw.indices()
			tmp_up_nb_atoms = len(tmp_up_indices)
			tmp_lw_nb_atoms = len(tmp_lw_indices)
			leaflet_sele["upper"]["all species"] = U.selectAtoms("bynum " + str(tmp_up_indices[0] + 1))
			leaflet_sele["lower"]["all species"] = U.selectAtoms("bynum " + str(tmp_lw_indices[0] + 1))
			for index in range(1,tmp_up_nb_atoms):
				leaflet_sele["upper"]["all species"] += U.selectAtoms("bynum " + str(tmp_up_indices[index] + 1))
				progress = '\r -identifying upper leaflet from gro file... ' + str(round(index/float(tmp_up_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			print ''
			for index in range(1,tmp_lw_nb_atoms):
				leaflet_sele["lower"]["all species"] += U.selectAtoms("bynum " + str(tmp_lw_indices[index] + 1))
				progress = '\r -identifying lower leaflet from gro file... ' + str(round(index/float(tmp_lw_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			leaflet_sele["both"]["all species"] = leaflet_sele["upper"]["all species"] + leaflet_sele["lower"]["all species"]
			print ''
		else:
			leaflet_sele["both"]["all species"] = U.selectAtoms(leaflet_sele_string)
			tmp_lipids_avg_z = leaflet_sele["both"]["all species"].centerOfGeometry()[2]
			leaflet_sele["upper"]["all species"] = leaflet_sele["both"]["all species"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
			leaflet_sele["lower"]["all species"] = leaflet_sele["both"]["all species"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))
		print " -found 2 leaflets: ", leaflet_sele["upper"]["all species"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"]["all species"].numberOfResidues(), '(lower) lipids'
			
	#store full selections
	for l in ["lower","upper","both"]:
		leaflet_sele_atoms[l]["all species"] = leaflet_sele[l]["all species"].residues.atoms
	
	return
def initialise_groups():

	global groups_labels
	global groups_nb
	global groups_boundaries
	global groups_sizes_dict

	groups_labels = {}
	groups_boundaries = {}
	groups_sizes_dict = {}

	#read file
	print "\nReading cluster groups definition file..."
	with open(args.cluster_groups_file) as f:
		lines = f.readlines()
	groups_nb = len(lines)
	for g_index in range(0,groups_nb):
		line = lines[g_index]
		if line[-1] == "\n":
			line = line[:-1]
		l_content = line.split(',')
		#check format
		if len(l_content) < 2:
			print "Error: the format of line " + str(g_index+1) + " should be 'min,max'."
			print "->", line
			sys.exit(1)
		tmp_beg = int(l_content[0])
		tmp_end = l_content[1]
		if tmp_end == "max":
			tmp_end = 100000											#put a stupidly big size to cap the open ended group
		else:
			tmp_end = int(tmp_end)
		groups_boundaries[g_index] = [tmp_beg,tmp_end]
		
	#display results
	print " -found " + str(groups_nb) + " cluster groups:"
	for g_index in range(0,groups_nb):
		if groups_boundaries[g_index][1] == 100000:
			print "   g" + str(g_index) + " = " + str(groups_boundaries[g_index][0]) + "+"
		else:
			print "   g" + str(g_index) + " = " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])

	#check for boundaries overlapping
	prev_beg = groups_boundaries[0][0]
	prev_end = groups_boundaries[0][1]
	if prev_end < prev_beg:
		print "Error: the max size is smaller than the min size for specified cluster groups " + str(groups_boundaries[0]) + "."
		sys.exit(1)
	for g_index in range(1,groups_nb):
		if groups_boundaries[g_index][1] < groups_boundaries[g_index][0]:
			print "Error: the max size is smaller than the min size for group " + str(g_index) + "(" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + ")."
			sys.exit(1)
		if groups_boundaries[g_index][0] <= prev_end:
			print "Error: specified cluster groups " + str(prev_beg) + "-" + str(prev_end) + " and " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + " overlap or are not in increasing order (boundaries are inclusive, see note 7 in bilayer_perturbations --help)."
			sys.exit(1)
		prev_beg = groups_boundaries[g_index][0]
		prev_end = groups_boundaries[g_index][1]
	
	#create equivalency table between groups and sizes
	for g_index in range(0,groups_nb):
		bb = groups_boundaries[g_index]
		tmp_beg = bb[0]
		tmp_end = bb[1]
		for tmp_size in range(tmp_beg, tmp_end+1):
			groups_sizes_dict[tmp_size] = g_index
	for tmp_size in list(set(range(1,max(groups_sizes_dict.keys()))) - set(groups_sizes_dict.keys())): 		#this handles potentially unaccounted for sizes up to the maximum specified by the user
		groups_sizes_dict[tmp_size] = groups_nb
	if max(groups_sizes_dict.keys())!= 100000:															    #this handles potentially unaccounted for sizes above the maximum specified by the user (in case it's not an open group)
		for tmp_size in range(max(groups_sizes_dict.keys())+1,100001):
			groups_sizes_dict[tmp_size] = groups_nb
				
	#create label for each group
	for g_index in range(0, groups_nb):
		if groups_boundaries[g_index][1] == 100000:
			groups_labels[g_index] = str(groups_boundaries[g_index][0]) + "+"
		elif groups_boundaries[g_index][0] == groups_boundaries[g_index][1]:
			groups_labels[g_index] = str(groups_boundaries[g_index][0])
		else:
			groups_labels[g_index] = str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])

	#add labels for the group "other", "lower" and "upper"
	groups_labels[groups_nb] = "other"
	groups_labels[-1] = "lower"
	groups_labels[groups_nb+1] = "upper"
		
	return

#=========================================================================================
# data structures
#=========================================================================================

def data_struct_time():													

	global frames_nb
	global frames_time
	frames_nb = np.zeros(nb_frames_to_process)
	frames_time = np.zeros(nb_frames_to_process)

	return
def data_ff_contacts():

	#create residues selection strings
	global residues_types_sele_string, residues_types_colours
	residues_types_colours = {}
	residues_types_sele_string = {}
	
	#TO DO: make these definitions user selectable
	type_res={}
	type_res["basic"]=['ARG','LYS']
	type_res["acidic"]=['ASP','GLU']
	type_res["polar"]=['SER','THR','ASN','GLN','HIS']
	type_res["hydrophobic"]=['VAL','ILE','LEU','MET','PHE','PRO','CYS','TYR','TRP']
	type_res["bb_only"]=['ALA','GLY']
	all_types = ["basic","polar","hydrophobic","bb_only"]

	residues_types_colours["basic"] = "b"
	residues_types_colours["acidic"] = "m"
	residues_types_colours["polar"] = "g"
	residues_types_colours["hydrophobic"] = "r"
	residues_types_colours["bb_only"] = "y"

	#create selection string for each type: residues
	residues_types_sele_string = {}
	for t in all_types:
		residues_types_sele_string[t] = "resname " + str(type_res[t][0])
		for r in type_res[t][1:]:
			residues_types_sele_string[t] += " or resname " + str(r)
	
	#get protein composition
	global protein_composition
	protein_composition = np.zeros(4)
	for t_index in range(0,len(all_types)):
		protein_composition[t_index] = proteins_sele[0].selectAtoms(residues_types_sele_string[all_types[t_index]]).numberOfAtoms()/float(nb_atom_per_protein)*100
	
	#structure for proteins cluster distribution
	global protein_TM_distribution_sizes
	protein_TM_distribution_sizes = np.zeros(proteins_nb)
	if args.cluster_groups_file != "no":
		global protein_TM_distribution_groups
		protein_TM_distribution_groups = np.zeros(groups_nb+1)
	
	#store contacts: global
	#----------------------
		
	#contacts stored in matrix where types are rows and cluster sizes are columns (#0:basic, 1:polar, 2: hydrophobic, 3:backbone)
	global lipids_ff_contacts_during_nb
	global lipids_ff_contacts_outside_nb
	global lipids_ff_contacts_during_pc
	global lipids_ff_contacts_outside_pc
	lipids_ff_contacts_during_nb = {}
	lipids_ff_contacts_outside_nb = {}
	lipids_ff_contacts_during_pc = {}
	lipids_ff_contacts_outside_pc = {}
	for l_index in range(0,lipids_ff_nb):
		lipids_ff_contacts_during_nb[l_index] = np.zeros((4,proteins_nb))
		lipids_ff_contacts_during_pc[l_index] = np.zeros((4,proteins_nb))
		lipids_ff_contacts_outside_nb[l_index] = np.zeros((4,proteins_nb))
		lipids_ff_contacts_outside_pc[l_index] = np.zeros((4,proteins_nb))
	#groups
	if args.cluster_groups_file != "no":
		global lipids_ff_contacts_during_nb_groups
		global lipids_ff_contacts_outside_nb_groups
		global lipids_ff_contacts_during_pc_groups
		global lipids_ff_contacts_outside_pc_groups
		lipids_ff_contacts_during_nb_groups = {}
		lipids_ff_contacts_outside_nb_groups = {}
		lipids_ff_contacts_during_pc_groups = {}
		lipids_ff_contacts_outside_pc_groups = {}
		for l_index in range(0,lipids_ff_nb):
			lipids_ff_contacts_during_nb_groups[l_index] = np.zeros((4,groups_nb+1))
			lipids_ff_contacts_during_pc_groups[l_index] = np.zeros((4,groups_nb+1))
			lipids_ff_contacts_outside_nb_groups[l_index] = np.zeros((4,groups_nb+1))
			lipids_ff_contacts_outside_pc_groups[l_index] = np.zeros((4,groups_nb+1))
		
	#store contacts: profile
	#-----------------------
	if args.profile:
		#relative leaflets positions
		global z_upper_avg, z_lower_avg, z_leaflets
		z_upper_avg = 0
		z_lower_avg = 0
		z_leaflets = 0

		#determine number of bins
		global bins_nb, bins_nb_max, bins_labels, bins_range
		bins_nb = int(np.floor(args.slices_dist/float(args.slices_thick))) 			#actually it's twice that as (-bins_nb,bins_nb) has to be filled
		bins_nb_max = bins_nb
		bins_labels = [str((n+0.5)*args.slices_thick) for n in range(-bins_nb,bins_nb)]
		bins_range = [(n+0.5)*args.slices_thick for n in range(-bins_nb,bins_nb)]
		
		#contacts stored in matrix, 1st dim = type, 2nd dim = z position, 3rd dim = cluster size
		global lipids_ff_contacts_during_nb_profile, lipids_ff_contacts_during_pc_profile
		global lipids_ff_contacts_outside_nb_profile, lipids_ff_contacts_outside_pc_profile
		lipids_ff_contacts_during_nb_profile = {}
		lipids_ff_contacts_during_pc_profile = {}
		lipids_ff_contacts_outside_nb_profile = {}
		lipids_ff_contacts_outside_pc_profile = {}
		for l_index in range(0,lipids_ff_nb):
			lipids_ff_contacts_during_nb_profile[l_index] = np.zeros((4, 2*bins_nb, proteins_nb))
			lipids_ff_contacts_during_pc_profile[l_index] = np.zeros((4, 2*bins_nb, proteins_nb))
			lipids_ff_contacts_outside_nb_profile[l_index] = np.zeros((4, 2*bins_nb, proteins_nb))
			lipids_ff_contacts_outside_pc_profile[l_index] = np.zeros((4, 2*bins_nb, proteins_nb))
		#same for groups
		if args.cluster_groups_file != "no":
			global lipids_ff_contacts_during_nb_groups_profile, lipids_ff_contacts_during_pc_groups_profile
			global lipids_ff_contacts_outside_nb_groups_profile, lipids_ff_contacts_outside_pc_groups_profile
			lipids_ff_contacts_during_nb_groups_profile = {}
			lipids_ff_contacts_during_pc_groups_profile = {}
			lipids_ff_contacts_outside_nb_groups_profile = {}
			lipids_ff_contacts_outside_pc_groups_profile = {}
			for l_index in range(0,lipids_ff_nb):
				lipids_ff_contacts_during_nb_groups_profile[l_index] = np.zeros((4, 2*bins_nb, groups_nb+1))
				lipids_ff_contacts_during_pc_groups_profile[l_index] = np.zeros((4, 2*bins_nb, groups_nb+1))
				lipids_ff_contacts_outside_nb_groups_profile[l_index] = np.zeros((4, 2*bins_nb, groups_nb+1))
				lipids_ff_contacts_outside_pc_groups_profile[l_index] = np.zeros((4, 2*bins_nb, groups_nb+1))

	#store statistics
	#----------------
	#calculate total nb of contacts
	global lipids_ff_contacts_u2l_during_tot_nb
	global lipids_ff_contacts_u2l_outside_tot_nb
	global lipids_ff_contacts_l2u_during_tot_nb
	global lipids_ff_contacts_l2u_outside_tot_nb
	lipids_ff_contacts_u2l_during_tot_nb = {}
	lipids_ff_contacts_u2l_outside_tot_nb = {}
	lipids_ff_contacts_l2u_during_tot_nb = {}
	lipids_ff_contacts_l2u_outside_tot_nb = {}

	#total nb and distribution of contacts with each cluster size
	global lipids_ff_contacts_u2l_during_tot_nb_by_size
	global lipids_ff_contacts_u2l_during_tot_pc_by_size
	global lipids_ff_contacts_u2l_outside_tot_nb_by_size
	global lipids_ff_contacts_u2l_outside_tot_pc_by_size
	global lipids_ff_contacts_l2u_during_tot_nb_by_size
	global lipids_ff_contacts_l2u_during_tot_pc_by_size
	global lipids_ff_contacts_l2u_outside_tot_nb_by_size
	global lipids_ff_contacts_l2u_outside_tot_pc_by_size
	lipids_ff_contacts_u2l_during_tot_nb_by_size = {}
	lipids_ff_contacts_u2l_during_tot_pc_by_size = {}
	lipids_ff_contacts_u2l_outside_tot_nb_by_size = {}
	lipids_ff_contacts_u2l_outside_tot_pc_by_size = {}
	lipids_ff_contacts_l2u_during_tot_nb_by_size = {}
	lipids_ff_contacts_l2u_during_tot_pc_by_size = {}
	lipids_ff_contacts_l2u_outside_tot_nb_by_size = {}
	lipids_ff_contacts_l2u_outside_tot_pc_by_size = {}

	#calculate total nb and distribution of contacts with each cluster group
	if args.cluster_groups_file != "no":
		global lipids_ff_contacts_u2l_during_tot_nb_by_size_group
		global lipids_ff_contacts_u2l_during_tot_pc_by_size_group
		global lipids_ff_contacts_u2l_outside_tot_nb_by_size_group
		global lipids_ff_contacts_u2l_outside_tot_pc_by_size_group
		global lipids_ff_contacts_l2u_during_tot_nb_by_size_group
		global lipids_ff_contacts_l2u_during_tot_pc_by_size_group
		global lipids_ff_contacts_l2u_outside_tot_nb_by_size_group
		global lipids_ff_contacts_l2u_outside_tot_pc_by_size_group
		lipids_ff_contacts_u2l_during_tot_nb_by_size_group = {}
		lipids_ff_contacts_u2l_during_tot_pc_by_size_group = {}
		lipids_ff_contacts_u2l_outside_tot_nb_by_size_group = {}
		lipids_ff_contacts_u2l_outside_tot_pc_by_size_group = {}
		lipids_ff_contacts_l2u_during_tot_nb_by_size_group = {}
		lipids_ff_contacts_l2u_during_tot_pc_by_size_group = {}
		lipids_ff_contacts_l2u_outside_tot_nb_by_size_group = {}
		lipids_ff_contacts_l2u_outside_tot_pc_by_size_group = {}

	#calculate total nb and distribution of contacts with each residue type
	global lipids_ff_contacts_u2l_during_tot_nb_by_type
	global lipids_ff_contacts_u2l_during_tot_pc_by_type
	global lipids_ff_contacts_u2l_outside_tot_nb_by_type
	global lipids_ff_contacts_u2l_outside_tot_pc_by_type
	global lipids_ff_contacts_l2u_during_tot_nb_by_type
	global lipids_ff_contacts_l2u_during_tot_pc_by_type
	global lipids_ff_contacts_l2u_outside_tot_nb_by_type
	global lipids_ff_contacts_l2u_outside_tot_pc_by_type
	lipids_ff_contacts_u2l_during_tot_nb_by_type = {}
	lipids_ff_contacts_u2l_during_tot_pc_by_type = {}
	lipids_ff_contacts_u2l_outside_tot_nb_by_type = {}
	lipids_ff_contacts_u2l_outside_tot_pc_by_type = {}
	lipids_ff_contacts_l2u_during_tot_nb_by_type = {}
	lipids_ff_contacts_l2u_during_tot_pc_by_type = {}
	lipids_ff_contacts_l2u_outside_tot_nb_by_type = {}
	lipids_ff_contacts_l2u_outside_tot_pc_by_type = {}
	
	return

#=========================================================================================
# core functions
#=========================================================================================

def get_z_coords(f_index):														
	
	tmp_zu = leaflet_sele["upper"]["all species"].centerOfGeometry()[2]
	tmp_zl = leaflet_sele["lower"]["all species"].centerOfGeometry()[2]
	tmp_zm = tmp_zl + (tmp_zu - tmp_zl)/float(2)
	z_upper[f_index] = tmp_zu - tmp_zm
	z_lower[f_index] = tmp_zl - tmp_zm
	for l in range(0,lipids_ff_nb):	
		z_ff[l][f_index] = lipids_sele_ff_bead[l].centerOfGeometry()[2] - tmp_zm

	return
def fit_coords_into_box(coords, box_dim):
	
	coords[:,0] -= np.floor(coords[:,0]/float(box_dim[0])) * box_dim[0]
	coords[:,1] -= np.floor(coords[:,1]/float(box_dim[1])) * box_dim[1]
	
	return coords
def get_distances(box_dim):												
	
	#method: use minimum distance between proteins
	#---------------------------------------------
	if args.m_algorithm == "min":
		#pre-process: get protein coordinates
		tmp_proteins_coords = np.zeros((proteins_nb, nb_atom_per_protein, 3))
		for p_index in range(0, proteins_nb):
			tmp_proteins_coords[p_index,:] = fit_coords_into_box(proteins_sele[p_index].coordinates(), box_dim)

		#store min distance between each proteins
		dist_matrix = 100000 * np.ones((proteins_nb,proteins_nb))
		for n in range(proteins_nb,1,-1):
			dist_matrix[proteins_nb-n,proteins_nb-n+1:proteins_nb] = map(lambda pp: np.min(MDAnalysis.analysis.distances.distance_array(np.float32(tmp_proteins_coords[proteins_nb-n,:]), np.float32(tmp_proteins_coords[pp,:]), box_dim)), range(proteins_nb-n+1,proteins_nb))
			dist_matrix[proteins_nb-n+1:proteins_nb,proteins_nb-n] = dist_matrix[proteins_nb-n,proteins_nb-n+1:proteins_nb]
											
	#method: use distance between cog
	#--------------------------------
	else:
		tmp_proteins_cogs = np.asarray(map(lambda p_index: calculate_cog(fit_coords_into_box(proteins_sele[p_index].coordinates(), box_dim), box_dim), range(0,proteins_nb)))
		dist_matrix = MDAnalysis.analysis.distances.distance_array(np.float32(tmp_proteins_cogs), np.float32(tmp_proteins_cogs), box_dim)

	return dist_matrix
def calculate_cog(tmp_coords, box_dim):										
	
	#this method allows to take pbc into account when calculcating the center of geometry 
	#see: http://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
	
	cog_coord = np.zeros(3)
	tmp_nb_atoms = np.shape(tmp_coords)[0]
	
	for n in range(0,3):
		tet = tmp_coords[:,n] * 2 * math.pi / float(box_dim[n])
		xsi = np.cos(tet)
		zet = np.sin(tet)
		tet_avg = math.atan2(-np.average(zet),-np.average(xsi)) + math.pi
		cog_coord[n] = tet_avg * box_dim[n] / float(2*math.pi)
	
	return cog_coord
def detect_clusters_connectivity(dist, box_dim):						
	
	#use networkx algorithm
	connected = (dist<args.nx_cutoff)
	network = nx.Graph(connected)
	groups = nx.connected_components(network)
	
	return groups
def detect_clusters_density(dist, box_dim):
	
	#run DBSCAN algorithm
	dbscan_output = DBSCAN(eps = args.dbscan_dist, metric = 'precomputed', min_samples = args.dbscan_nb).fit(dist)

	#build 'groups' structure i.e. a list whose element are all the clusters identified
	groups = []
	for c_lab in np.unique(dbscan_output.labels_):
		tmp_pos = np.argwhere(dbscan_output.labels_ == c_lab)
		if c_lab == -1:
			groups += map(lambda p:p[0] , tmp_pos)
		else:
			groups.append(map(lambda p:p[0] , tmp_pos))

	return groups
def identify_ff_contacts(box_dim, f_time, f_nb):

	global protein_TM_distribution_sizes
	global lipids_ff_contacts_during_nb, lipids_ff_contacts_outside_nb
	if args.cluster_groups_file != "no":
		global protein_TM_distribution_groups
	if args.profile:
		global z_upper_avg, z_lower_avg, z_leaflets
		loc_z_axis = np.array([0,0,1])
		loc_z_axis = loc_z_axis.reshape((3,1))

	#initialise array allowing to retrieve cluster size from atom indices
	tmp_atindices_2_csize = np.zeros(proteins_max_at_number+1)
	tmp_atindices_2_csele = {}
	
	#retrieve coordinates arrays (pre-processing saves time as MDAnalysis functions are quite slow and we need to make such calls a few times)
	tmp_lip_coords = {l: fit_coords_into_box(leaflet_sele[l]["all species"].coordinates(), box_dim) for l in ["lower","upper"]}
	
	#identify clusters
	#=================	
	if args.m_algorithm != "density":
		clusters = detect_clusters_connectivity(get_distances(box_dim), box_dim)
	else:
		clusters = detect_clusters_density(get_distances(box_dim), box_dim)
	
	#process each cluster
	#====================
	c_sele_all = MDAnalysis.core.AtomGroup.AtomGroup([])
	nb_clusters = len(clusters)
	c_counter = 0
	for cluster in clusters:		
		#display update
		c_counter += 1
		progress = '\r -processing frame ' + str(f_nb+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' from ' + str(f_start) + ' to ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ') and cluster ' + str(c_counter) + '/' + str(nb_clusters) + '              '
		sys.stdout.flush()
		sys.stdout.write(progress)

		#create selection for current cluster and only process it if it's TM (find closest PO4 particles for each particles of clusters, if all are in the same leaflet then it's surfacic [NB: this is done at the CLUSTER level (the same criteria at the protein level would probably fail)])
		c_sele = MDAnalysis.core.AtomGroup.AtomGroup([])
		for p_index in cluster:
			c_sele += proteins_sele[p_index]
		c_sele_all += c_sele
		tmp_c_sele_coordinates = fit_coords_into_box(c_sele.coordinates(), box_dim)
		dist_min_lower = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["lower"], box_dim), axis = 1)
		dist_min_upper = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["upper"], box_dim), axis = 1)
		dist = dist_min_upper - dist_min_lower
		#case: interfacial lower
		if np.size(dist[dist>0]) == np.size(dist):
			tmp_atindices_2_csize[c_sele.indices()] = -1
		#case: interfacial upper
		elif np.size(dist[dist>0]) == 0:
			tmp_atindices_2_csize[c_sele.indices()] = 99999				#stupidly big number
		#case: TM
		else:
			c_size = np.size(cluster)
			tmp_atindices_2_csize[c_sele.indices()] = c_size
			tmp_atindices_2_csele.update(dict.fromkeys(c_sele.indices(), tmp_c_sele_coordinates))
			protein_TM_distribution_sizes[c_size-1] += c_size
			if args.cluster_groups_file != "no":
				protein_TM_distribution_groups[groups_sizes_dict[c_size]] += c_size
				
	#process each ff lipid
	#=====================
	ff_counter = 0
	for l_index in range(0,lipids_ff_nb):
		
		#display update
		ff_counter += 1
		progress = '\r -processing frame ' + str(f_nb+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' from ' + str(f_start) + ' to ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ') and flip-flopping lipid ' + str(ff_counter) + '/' + str(lipids_ff_nb) + '   '
		sys.stdout.flush()
		sys.stdout.write(progress)

		#detect contacts
		ff_lip_and_prot_TM = lipids_sele_ff[l_index] + c_sele_all
		around_lip_prot_TM = ff_lip_and_prot_TM.selectAtoms("around " + str(args.cutoff_pl) + " (resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1]) + ")")	
		
		#get size of cluster in contact if any
		if around_lip_prot_TM.numberOfAtoms() > 0:			
			tmp_ctct_at_number = around_lip_prot_TM.atoms[0].number
			tmp_size = tmp_atindices_2_csize[tmp_ctct_at_number]
			tmp_nbct = around_lip_prot_TM.numberOfAtoms()
						
			#store contacts if cluster TM
			if tmp_size > 0 and tmp_size < 99999:
				
				#get size group 
				#--------------
				if args.cluster_groups_file != "no":
					tmp_group = groups_sizes_dict[tmp_size]

				#store overal number of each contact types
				#-----------------------------------------
				tmp_ctct_basic = around_lip_prot_TM.selectAtoms(residues_types_sele_string['basic']).numberOfAtoms()
				tmp_ctct_polar = around_lip_prot_TM.selectAtoms(residues_types_sele_string['polar']).numberOfAtoms()
				tmp_ctct_hydrophobic = around_lip_prot_TM.selectAtoms(residues_types_sele_string['hydrophobic']).numberOfAtoms()
				tmp_ctct_bb_only = around_lip_prot_TM.selectAtoms(residues_types_sele_string['bb_only']).numberOfAtoms()
				#outside ff
				if f_time < lipids_ff_info[l_index][4] or f_time > lipids_ff_info[l_index][5]:
					lipids_ff_contacts_outside_nb[l_index][0,tmp_size - 1] += tmp_ctct_basic
					lipids_ff_contacts_outside_nb[l_index][1,tmp_size - 1] += tmp_ctct_polar
					lipids_ff_contacts_outside_nb[l_index][2,tmp_size - 1] += tmp_ctct_hydrophobic
					lipids_ff_contacts_outside_nb[l_index][3,tmp_size - 1] += tmp_ctct_bb_only
					if args.cluster_groups_file != "no":
						lipids_ff_contacts_outside_nb_groups[l_index][0,tmp_group] += tmp_ctct_basic
						lipids_ff_contacts_outside_nb_groups[l_index][1,tmp_group] += tmp_ctct_polar
						lipids_ff_contacts_outside_nb_groups[l_index][2,tmp_group] += tmp_ctct_hydrophobic
						lipids_ff_contacts_outside_nb_groups[l_index][3,tmp_group] += tmp_ctct_bb_only				
				#during ff
				else:
					lipids_ff_contacts_during_nb[l_index][0,tmp_size - 1] += tmp_ctct_basic
					lipids_ff_contacts_during_nb[l_index][1,tmp_size - 1] += tmp_ctct_polar
					lipids_ff_contacts_during_nb[l_index][2,tmp_size - 1] += tmp_ctct_hydrophobic					
					lipids_ff_contacts_during_nb[l_index][3,tmp_size - 1] += tmp_ctct_bb_only
					if args.cluster_groups_file != "no":
						lipids_ff_contacts_during_nb_groups[l_index][0,tmp_group] += tmp_ctct_basic
						lipids_ff_contacts_during_nb_groups[l_index][1,tmp_group] += tmp_ctct_polar
						lipids_ff_contacts_during_nb_groups[l_index][2,tmp_group] += tmp_ctct_hydrophobic					
						lipids_ff_contacts_during_nb_groups[l_index][3,tmp_group] += tmp_ctct_bb_only

				#store them along with local depth of insertion of the lipid
				#-----------------------------------------------------------
				if args.profile:
					#get ff bead coordinate
					ff_bead = fit_coords_into_box(lipids_sele_ff[l_index].coordinates(), box_dim)
					
					#calculate local normal to bilayer
					if args.normal != 'z':
						#identify local normal to bilayer using cog of cluster in contact as reference
						cluster_cog = calculate_cog(tmp_atindices_2_csele[tmp_ctct_at_number], box_dim)
	
						#switch to flip-flop bead referential
						tmp_lip_coords_up = tmp_lip_coords["upper"] - cluster_cog
						tmp_lip_coords_lw = tmp_lip_coords["lower"] - cluster_cog
						
						#identify neighbouring particles in each leaflet
						tmp_lip_coords_up_within = tmp_lip_coords_up[tmp_lip_coords_up[:,0]**2 + tmp_lip_coords_up[:,1]**2 + tmp_lip_coords_up[:,2]**2 < args.normal_d**2]
						tmp_lip_coords_lw_within = tmp_lip_coords_lw[tmp_lip_coords_lw[:,0]**2 + tmp_lip_coords_lw[:,1]**2 + tmp_lip_coords_lw[:,2]**2 < args.normal_d**2]
						if np.shape(tmp_lip_coords_up_within)[0] == 0:
							print "\nWarning: no neighbouring particles found in the upper leaflet for current cluster (size " + str(int(tmp_size)) + "). Check the --normal_d option.\n"						
							#debug
							print "cluster cog", cluster_cog
							print "time: ", f_time
							print tmp_ctct_at_number
							print tmp_atindices_2_csele[tmp_ctct_at_number]
							toto = tmp_atindices_2_csele[tmp_ctct_at_number]
							print toto.indices()
							continue
						else:
							cog_up = np.average(tmp_lip_coords_up_within, axis = 0)
						if np.shape(tmp_lip_coords_lw_within)[0] == 0:
							print "\nWarning: no neighbouring particles found in the lower leaflet for current cluster (size " + str(int(tmp_size)) + "). Check the --normal_d option.\n"
							continue
						else:
							cog_lw = np.average(tmp_lip_coords_lw_within, axis = 0)
						
						#identify normal vector: case cog
						if args.normal == 'cog':
							norm_vec = cog_up - cog_lw
							norm_vec /= float(np.linalg.norm(norm_vec))
							norm_vec = norm_vec.reshape((3,1))
						#identify normal vector: case svd
						else:
							tmp_lip_coords_within = np.concatenate((tmp_lip_coords_up_within-cog_up,tmp_lip_coords_lw_within-cog_lw))
							svd_U, svd_D, svd_V = np.linalg.svd(tmp_lip_coords_within)
							norm_vec = svd_V[2].reshape((3,1))
							#orientate the normal vector so that it goes from inside (lower) to outside (upper) (IMPORTANT: ensures correct + sign convention)
							tmp_delta_cog = cog_up - cog_lw
							tmp_delta_cog = tmp_delta_cog.reshape((3,1))
							if np.dot(norm_vec[:,0],tmp_delta_cog[:,0]) < 0:
								norm_vec *= -1
	
						#identify rotation matrix
						norm_ax = np.cross(loc_z_axis,norm_vec,axis=0)
						norm_cos = np.dot(loc_z_axis[:,0],norm_vec[:,0])
						norm_sin = np.linalg.norm(norm_ax)
						norm_ax_skew_sym = norm_vec*loc_z_axis.T - loc_z_axis*norm_vec.T
						norm_rot = np.identity(3) - norm_ax_skew_sym + (1-norm_cos)/float(norm_sin**2)*np.dot(norm_ax_skew_sym,norm_ax_skew_sym)
					
						#rotate neighbouring bilayer in local cluster referential
						tmp_lip_coords_up_within_rotated = np.dot(norm_rot, tmp_lip_coords_up_within.T).T
						tmp_lip_coords_lw_within_rotated = np.dot(norm_rot, tmp_lip_coords_lw_within.T).T
						
						#identify z coord of local middle of bilayer after rotation
						cog_up_rotated = np.median(tmp_lip_coords_up_within_rotated, axis = 0)
						cog_lw_rotated = np.median(tmp_lip_coords_lw_within_rotated, axis = 0)
						norm_z_middle = cog_lw_rotated[2] + (cog_up_rotated[2] - cog_lw_rotated[2])/float(2)

						#store local relative positions of leaflets
						z_upper_avg += cog_up_rotated[2] - norm_z_middle
						z_lower_avg += cog_lw_rotated[2] - norm_z_middle
						z_leaflets += 1

						#rotate coordinates of ff lipid
						ff_bead -= cluster_cog
						ff_bead = np.dot(norm_rot, ff_bead.T).T

					#or simply assume use z axis
					else:
						tmp_zu = np.average(tmp_lip_coords["upper"], axis = 0)[2]
						tmp_zl = np.average(tmp_lip_coords["lower"], axis = 0)[2]
						norm_z_middle = tmp_zl + (tmp_zu - tmp_zl)/float(2)	
	
						#store local relative positions of leaflets
						z_upper_avg += tmp_zu - norm_z_middle
						z_lower_avg += tmp_zl - norm_z_middle
						z_lealflets += 1
						
					#find which bin along the normal between of local bilayer the ff bead is at
					bin_rel = np.floor((ff_bead[:,2] - norm_z_middle)/float(args.slices_thick)).astype(int)
					if abs(bin_rel) < bins_nb_max:
						#the + int(bins_nb) allows to only have positive bin indices
						bin_abs = bin_rel + int(bins_nb)
								
						#add the nb of contacts of each type for that bin in each contact type struct
						if f_time < lipids_ff_info[l_index][4] or f_time > lipids_ff_info[l_index][5]:
							lipids_ff_contacts_outside_nb_profile[l_index][0, bin_abs, tmp_size - 1] += tmp_ctct_basic
							lipids_ff_contacts_outside_nb_profile[l_index][1, bin_abs, tmp_size - 1] += tmp_ctct_polar
							lipids_ff_contacts_outside_nb_profile[l_index][2, bin_abs, tmp_size - 1] += tmp_ctct_hydrophobic
							lipids_ff_contacts_outside_nb_profile[l_index][3, bin_abs, tmp_size - 1] += tmp_ctct_bb_only
							if args.cluster_groups_file != "no":
								lipids_ff_contacts_outside_nb_groups_profile[l_index][0, bin_abs, tmp_group] += tmp_ctct_basic
								lipids_ff_contacts_outside_nb_groups_profile[l_index][1, bin_abs, tmp_group] += tmp_ctct_polar
								lipids_ff_contacts_outside_nb_groups_profile[l_index][2, bin_abs, tmp_group] += tmp_ctct_hydrophobic
								lipids_ff_contacts_outside_nb_groups_profile[l_index][3, bin_abs, tmp_group] += tmp_ctct_bb_only
						else:
							lipids_ff_contacts_during_nb_profile[l_index][0, bin_abs, tmp_size - 1] += tmp_ctct_basic
							lipids_ff_contacts_during_nb_profile[l_index][1, bin_abs, tmp_size - 1] += tmp_ctct_polar
							lipids_ff_contacts_during_nb_profile[l_index][2, bin_abs, tmp_size - 1] += tmp_ctct_hydrophobic
							lipids_ff_contacts_during_nb_profile[l_index][3, bin_abs, tmp_size - 1] += tmp_ctct_bb_only
							if args.cluster_groups_file != "no":
								lipids_ff_contacts_during_nb_groups_profile[l_index][0, bin_abs, tmp_group] += tmp_ctct_basic
								lipids_ff_contacts_during_nb_groups_profile[l_index][1, bin_abs, tmp_group] += tmp_ctct_polar
								lipids_ff_contacts_during_nb_groups_profile[l_index][2, bin_abs, tmp_group] += tmp_ctct_hydrophobic
								lipids_ff_contacts_during_nb_groups_profile[l_index][3, bin_abs, tmp_group] += tmp_ctct_bb_only
								
	return

#=========================================================================================
# statistics
#=========================================================================================

def calc_stats_ctcts():
	
	if args.profile:
		global z_lower_avg, z_upper_avg
		z_upper_avg /= float(z_leaflets)
		z_lower_avg /= float(z_leaflets)
	
	#calculate TM proteins distribution
	#==================================
	global protein_TM_distribution_sizes
	protein_TM_distribution_sizes = protein_TM_distribution_sizes / float(np.sum(protein_TM_distribution_sizes)) * 100
	protein_TM_distribution_sizes = protein_TM_distribution_sizes[0:np.amax(np.nonzero(protein_TM_distribution_sizes))+1]	
	
	#identify index of cluster sizes sampled
	global protein_sizes_sampled, protein_max_size_sampled
	protein_sizes_sampled = np.arange(0,proteins_nb)[protein_TM_distribution_sizes != 0]
	protein_max_size_sampled = len(protein_TM_distribution_sizes)
		
	#check whether the "other" group has been sampled
	if args.cluster_groups_file != "no":
		global protein_TM_distribution_groups
		global group_gmax
		protein_TM_distribution_groups = protein_TM_distribution_groups / float(np.sum(protein_TM_distribution_groups)) * 100
		if protein_TM_distribution_groups[groups_nb] != 0:
			group_gmax = groups_nb + 1
		else:
			group_gmax = groups_nb

	#create data structure to store data in order to calculate average contact profiles
	if args.profile:
		tmp_profile_u2l_during_sizes = {}
		tmp_profile_u2l_outside_sizes = {}
		tmp_profile_l2u_during_sizes = {}
		tmp_profile_l2u_outside_sizes = {}
		if args.cluster_groups_file != "no":
			tmp_profile_u2l_during_groups = {g_index: {} for g_index in range(0, group_gmax)}
			tmp_profile_u2l_outside_groups = {g_index: {} for g_index in range(0, group_gmax)}
			tmp_profile_l2u_during_groups = {g_index: {} for g_index in range(0, group_gmax)}
			tmp_profile_l2u_outside_groups = {g_index: {} for g_index in range(0, group_gmax)}
	
	# for each flip-flopping lipids: u2l
	#==============================
	for l_index in lipids_ff_u2l_index:
		#calculate total nb of contacts
		#------------------------------
		lipids_ff_contacts_u2l_during_tot_nb[l_index] = np.sum(lipids_ff_contacts_during_nb[l_index])
		lipids_ff_contacts_u2l_outside_tot_nb[l_index] = np.sum(lipids_ff_contacts_outside_nb[l_index])

		#calculate total nb and distribution of contacts with each cluster size (and group)
		#----------------------------------------------------------------------------------
		lipids_ff_contacts_u2l_during_tot_nb_by_size[l_index] = np.sum(lipids_ff_contacts_during_nb[l_index], axis = 0)
		lipids_ff_contacts_u2l_during_tot_pc_by_size[l_index] = lipids_ff_contacts_u2l_during_tot_nb_by_size[l_index] / float(lipids_ff_contacts_u2l_during_tot_nb[l_index]) * 100
		lipids_ff_contacts_u2l_outside_tot_nb_by_size[l_index] = np.sum(lipids_ff_contacts_outside_nb[l_index], axis = 0)
		lipids_ff_contacts_u2l_outside_tot_pc_by_size[l_index] = lipids_ff_contacts_u2l_outside_tot_nb_by_size[l_index] / float(lipids_ff_contacts_u2l_outside_tot_nb[l_index]) * 100
		if args.cluster_groups_file != "no":
			lipids_ff_contacts_u2l_during_tot_nb_by_size_group[l_index] = np.sum(lipids_ff_contacts_during_nb_groups[l_index], axis = 0)
			lipids_ff_contacts_u2l_during_tot_pc_by_size_group[l_index] = lipids_ff_contacts_u2l_during_tot_nb_by_size_group[l_index] / float(lipids_ff_contacts_u2l_during_tot_nb[l_index]) * 100
			lipids_ff_contacts_u2l_outside_tot_nb_by_size_group[l_index] = np.sum(lipids_ff_contacts_outside_nb_groups[l_index], axis = 0)
			lipids_ff_contacts_u2l_outside_tot_pc_by_size_group[l_index] = lipids_ff_contacts_u2l_outside_tot_nb_by_size_group[l_index] / float(lipids_ff_contacts_u2l_outside_tot_nb[l_index]) * 100
			
		#calculate total nb and distribution of contacts with each residue type
		#----------------------------------------------------------------------
		lipids_ff_contacts_u2l_during_tot_nb_by_type[l_index] = np.sum(lipids_ff_contacts_during_nb[l_index], axis = 1)
		lipids_ff_contacts_u2l_during_tot_pc_by_type[l_index] = lipids_ff_contacts_u2l_during_tot_nb_by_type[l_index] / float(lipids_ff_contacts_u2l_during_tot_nb[l_index]) * 100
		lipids_ff_contacts_u2l_outside_tot_nb_by_type[l_index] = np.sum(lipids_ff_contacts_outside_nb[l_index], axis = 1)
		lipids_ff_contacts_u2l_outside_tot_pc_by_type[l_index] = lipids_ff_contacts_u2l_outside_tot_nb_by_type[l_index] / float(lipids_ff_contacts_u2l_outside_tot_nb[l_index]) * 100
		
		#calculate distribution of contacts over sizes/groups for each type
		#------------------------------------------------------------------
		for t in range(0,4):
			lipids_ff_contacts_during_pc[l_index][t,:] = 100 * lipids_ff_contacts_during_nb[l_index][t,:] / float(np.sum(lipids_ff_contacts_during_nb[l_index][t,:]))
			lipids_ff_contacts_outside_pc[l_index][t,:] = 100 * lipids_ff_contacts_outside_nb[l_index][t,:] / float(np.sum(lipids_ff_contacts_outside_nb[l_index][t,:]))
			if args.cluster_groups_file != "no":
				lipids_ff_contacts_during_pc_groups[l_index][t,:] = 100 * lipids_ff_contacts_during_nb_groups[l_index][t,:] / float(np.sum(lipids_ff_contacts_during_nb[l_index][t,:]))
				lipids_ff_contacts_outside_pc_groups[l_index][t,:] = 100 * lipids_ff_contacts_outside_nb_groups[l_index][t,:] / float(np.sum(lipids_ff_contacts_outside_nb[l_index][t,:]))

		#calculate distribution of contacts over z for each type
		#-------------------------------------------------------
		if args.profile: 
			#for each size
			for c_size in protein_sizes_sampled:
				tmp_tot_c_size_during = np.sum(lipids_ff_contacts_during_nb_profile[l_index][:, :, c_size])
				tmp_tot_c_size_outside = np.sum(lipids_ff_contacts_outside_nb_profile[l_index][:, :, c_size])
				#during ff
				if  tmp_tot_c_size_during > 0:
					#calculate distribution
					for t in range(0,4):
						lipids_ff_contacts_during_pc_profile[l_index][t, :, c_size] = 100 * lipids_ff_contacts_during_nb_profile[l_index][t, :, c_size] / float(tmp_tot_c_size_during)
					#store it in a format allowing easy averaging over ff lipids later
					tmp_profile_u2l_during_sizes[l_index] = np.zeros((4,2*bins_nb))
					tmp_profile_u2l_during_sizes[l_index][:,:] = lipids_ff_contacts_during_pc_profile[l_index][:, :, c_size]
				#outside ff
				if  tmp_tot_c_size_outside > 0:
					#calculate distribution
					for t in range(0,4):
						lipids_ff_contacts_outside_pc_profile[l_index][t, :, c_size] = 100 * lipids_ff_contacts_outside_nb_profile[l_index][t, :, c_size] / float(tmp_tot_c_size_outside)
					#store it in a format allowing easy averaging over ff lipids later
					tmp_profile_u2l_outside_sizes[l_index] = np.zeros((4,2*bins_nb))
					tmp_profile_u2l_outside_sizes[l_index][:,:] = lipids_ff_contacts_outside_pc_profile[l_index][:, :, c_size]

			#for each group
			if args.cluster_groups_file != "no":
				for g_index in range(0,group_gmax):
					tmp_tot_g_index_during = np.sum(lipids_ff_contacts_during_nb_groups_profile[l_index][:, :, g_index])
					tmp_tot_g_index_outside = np.sum(lipids_ff_contacts_outside_nb_groups_profile[l_index][:, :, g_index])
					#during ff
					if tmp_tot_g_index_during > 0:
						#calculate distribution
						for t in range(0,4):
							lipids_ff_contacts_during_pc_groups_profile[l_index][t, :, g_index] = 100 * lipids_ff_contacts_during_nb_groups_profile[l_index][t, :, g_index] / float(tmp_tot_g_index_during)
						#store it in a format allowing easy averaging over ff lipids later
						tmp_profile_u2l_during_groups[g_index][l_index] = np.zeros((4,2*bins_nb))
						tmp_profile_u2l_during_groups[g_index][l_index][:,:] = lipids_ff_contacts_during_pc_groups_profile[l_index][:, :, g_index]		
					#outside ff
					if tmp_tot_g_index_outside > 0:
						#calculate distribution
						for t in range(0,4):
							lipids_ff_contacts_outside_pc_groups_profile[l_index][t, :, g_index] = 100 * lipids_ff_contacts_outside_nb_groups_profile[l_index][t, :, g_index] / float(tmp_tot_g_index_outside)
						#store it in a format allowing easy averaging over ff lipids later
						tmp_profile_u2l_outside_groups[g_index][l_index] = np.zeros((4,2*bins_nb))
						tmp_profile_u2l_outside_groups[g_index][l_index][:,:] = lipids_ff_contacts_outside_pc_groups_profile[l_index][:, :, g_index]

	# for each flip-flopping lipids: l2u
	#==============================
	for l_index in lipids_ff_l2u_index:
		#calculate total nb of contacts
		#------------------------------
		lipids_ff_contacts_l2u_during_tot_nb[l_index] = np.sum(lipids_ff_contacts_during_nb[l_index])
		lipids_ff_contacts_l2u_outside_tot_nb[l_index] = np.sum(lipids_ff_contacts_outside_nb[l_index])

		#calculate total nb and distribution of contacts with each cluster size (and group)
		#----------------------------------------------------------------------------------
		lipids_ff_contacts_l2u_during_tot_nb_by_size[l_index] = np.sum(lipids_ff_contacts_during_nb[l_index], axis = 0)
		lipids_ff_contacts_l2u_during_tot_pc_by_size[l_index] = lipids_ff_contacts_l2u_during_tot_nb_by_size[l_index] / float(lipids_ff_contacts_l2u_during_tot_nb[l_index]) * 100
		lipids_ff_contacts_l2u_outside_tot_nb_by_size[l_index] = np.sum(lipids_ff_contacts_outside_nb[l_index], axis = 0)
		lipids_ff_contacts_l2u_outside_tot_pc_by_size[l_index] = lipids_ff_contacts_l2u_outside_tot_nb_by_size[l_index] / float(lipids_ff_contacts_l2u_outside_tot_nb[l_index]) * 100
		if args.cluster_groups_file != "no":
			lipids_ff_contacts_l2u_during_tot_nb_by_size_group[l_index] = np.sum(lipids_ff_contacts_during_nb_groups[l_index], axis = 0)
			lipids_ff_contacts_l2u_during_tot_pc_by_size_group[l_index] = lipids_ff_contacts_l2u_during_tot_nb_by_size_group[l_index] / float(lipids_ff_contacts_l2u_during_tot_nb[l_index]) * 100
			lipids_ff_contacts_l2u_outside_tot_nb_by_size_group[l_index] = np.sum(lipids_ff_contacts_outside_nb_groups[l_index], axis = 0)
			lipids_ff_contacts_l2u_outside_tot_pc_by_size_group[l_index] = lipids_ff_contacts_l2u_outside_tot_nb_by_size_group[l_index] / float(lipids_ff_contacts_l2u_outside_tot_nb[l_index]) * 100

		#calculate total nb and distribution of contacts with each residue type
		#----------------------------------------------------------------------
		lipids_ff_contacts_l2u_during_tot_nb_by_type[l_index] = np.sum(lipids_ff_contacts_during_nb[l_index], axis = 1)
		lipids_ff_contacts_l2u_during_tot_pc_by_type[l_index] = lipids_ff_contacts_l2u_during_tot_nb_by_type[l_index] / float(lipids_ff_contacts_l2u_during_tot_nb[l_index]) * 100
		lipids_ff_contacts_l2u_outside_tot_nb_by_type[l_index] = np.sum(lipids_ff_contacts_outside_nb[l_index], axis = 1)
		lipids_ff_contacts_l2u_outside_tot_pc_by_type[l_index] = lipids_ff_contacts_l2u_outside_tot_nb_by_type[l_index] / float(lipids_ff_contacts_l2u_outside_tot_nb[l_index]) * 100
		
		#calculate distribution of contacts over sizes/groups for each type
		#------------------------------------------------------------------
		for t in range(0,4):
			lipids_ff_contacts_during_pc[l_index][t,:] = 100 * lipids_ff_contacts_during_nb[l_index][t,:] / float(np.sum(lipids_ff_contacts_during_nb[l_index][t,:]))
			lipids_ff_contacts_outside_pc[l_index][t,:] = 100 * lipids_ff_contacts_outside_nb[l_index][t,:] / float(np.sum(lipids_ff_contacts_outside_nb[l_index][t,:]))
			if args.cluster_groups_file != "no":
				lipids_ff_contacts_during_pc_groups[l_index][t,:] = 100 * lipids_ff_contacts_during_nb_groups[l_index][t,:] / float(np.sum(lipids_ff_contacts_during_nb[l_index][t,:]))
				lipids_ff_contacts_outside_pc_groups[l_index][t,:] = 100 * lipids_ff_contacts_outside_nb_groups[l_index][t,:] / float(np.sum(lipids_ff_contacts_outside_nb[l_index][t,:]))

		#calculate distribution of contacts over z for each type
		#-------------------------------------------------------
		if args.profile: 
			#for each size
			for c_size in protein_sizes_sampled:
				tmp_tot_c_size_during = np.sum(lipids_ff_contacts_during_nb_profile[l_index][:, :, c_size])
				tmp_tot_c_size_outside = np.sum(lipids_ff_contacts_outside_nb_profile[l_index][:, :, c_size])
				#during ff
				if  tmp_tot_c_size_during > 0:
					#calculate distribution
					for t in range(0,4):
						lipids_ff_contacts_during_pc_profile[l_index][t, :, c_size] = 100 * lipids_ff_contacts_during_nb_profile[l_index][t, :, c_size] / float(tmp_tot_c_size_during)
					#store it in a format allowing easy averaging over ff lipids later
					tmp_profile_l2u_during_sizes[l_index] = np.zeros((4,2*bins_nb))
					tmp_profile_l2u_during_sizes[l_index][:,:] = lipids_ff_contacts_during_pc_profile[l_index][:, :, c_size]				
				#outside ff
				if  tmp_tot_c_size_outside > 0:
					#calculate distribution
					for t in range(0,4):
						lipids_ff_contacts_outside_pc_profile[l_index][t, :, c_size] = 100 * lipids_ff_contacts_outside_nb_profile[l_index][t, :, c_size] / float(tmp_tot_c_size_outside)
					#store it in a format allowing easy averaging over ff lipids later
					tmp_profile_l2u_outside_sizes[l_index] = np.zeros((4,2*bins_nb))
					tmp_profile_l2u_outside_sizes[l_index][:,:] = lipids_ff_contacts_outside_pc_profile[l_index][:, :, c_size]
			
			#for each group
			if args.cluster_groups_file != "no":
				for g_index in range(0,group_gmax):
					tmp_tot_g_index_during = np.sum(lipids_ff_contacts_during_nb_groups_profile[l_index][:, :, g_index])
					tmp_tot_g_index_outside = np.sum(lipids_ff_contacts_outside_nb_groups_profile[l_index][:, :, g_index])
					#during ff
					if tmp_tot_g_index_during > 0:
						#calculate distribution
						for t in range(0,4):
							lipids_ff_contacts_during_pc_groups_profile[l_index][t, :, g_index] = 100 * lipids_ff_contacts_during_nb_groups_profile[l_index][t, :, g_index] / float(tmp_tot_g_index_during)
						#store it in a format allowing easy averaging over ff lipids later
						tmp_profile_l2u_during_groups[g_index][l_index] = np.zeros((4,2*bins_nb))
						tmp_profile_l2u_during_groups[g_index][l_index][:,:] = lipids_ff_contacts_during_pc_groups_profile[l_index][:, :, g_index]		
					#outside ff
					if tmp_tot_g_index_outside > 0:
						#calculate distribution
						for t in range(0,4):
							lipids_ff_contacts_outside_pc_groups_profile[l_index][t, :, g_index] = 100 * lipids_ff_contacts_outside_nb_groups_profile[l_index][t, :, g_index] / float(tmp_tot_g_index_outside)
						#store it in a format allowing easy averaging over ff lipids later
						tmp_profile_l2u_outside_groups[g_index][l_index] = np.zeros((4,2*bins_nb))
						tmp_profile_l2u_outside_groups[g_index][l_index][:,:] = lipids_ff_contacts_outside_pc_groups_profile[l_index][:, :, g_index]

	# averages
	#=========
	#distribution of contacts over sizes
	global lipids_ff_contacts_u2l_during_by_size_avg, lipids_ff_contacts_u2l_during_by_size_std
	global lipids_ff_contacts_l2u_during_by_size_avg, lipids_ff_contacts_l2u_during_by_size_std
	global lipids_ff_contacts_u2l_outside_by_size_avg, lipids_ff_contacts_u2l_outside_by_size_std
	global lipids_ff_contacts_l2u_outside_by_size_avg, lipids_ff_contacts_l2u_outside_by_size_std
	lipids_ff_contacts_u2l_during_by_size_avg = np.average(lipids_ff_contacts_u2l_during_tot_pc_by_size.values(), axis = 0)
	lipids_ff_contacts_u2l_during_by_size_std = np.std(lipids_ff_contacts_u2l_during_tot_pc_by_size.values(), axis = 0)
	lipids_ff_contacts_u2l_outside_by_size_avg = np.average(lipids_ff_contacts_u2l_outside_tot_pc_by_size.values(), axis = 0)
	lipids_ff_contacts_u2l_outside_by_size_std = np.std(lipids_ff_contacts_u2l_outside_tot_pc_by_size.values(), axis = 0)
	lipids_ff_contacts_l2u_during_by_size_avg = np.average(lipids_ff_contacts_l2u_during_tot_pc_by_size.values(), axis = 0)
	lipids_ff_contacts_l2u_during_by_size_std = np.std(lipids_ff_contacts_l2u_during_tot_pc_by_size.values(), axis = 0)
	lipids_ff_contacts_l2u_outside_by_size_avg = np.average(lipids_ff_contacts_l2u_outside_tot_pc_by_size.values(), axis = 0)
	lipids_ff_contacts_l2u_outside_by_size_std = np.std(lipids_ff_contacts_l2u_outside_tot_pc_by_size.values(), axis = 0)

	#new
	global lipids_ff_contacts_u2l_during_by_size_nb_sum, lipids_ff_contacts_u2l_outside_by_size_nb_sum
	lipids_ff_contacts_u2l_during_by_size_nb_sum = np.sum(lipids_ff_contacts_u2l_during_tot_nb_by_size.values(), axis = 0)
	lipids_ff_contacts_u2l_outside_by_size_nb_sum = np.sum(lipids_ff_contacts_u2l_outside_tot_nb_by_size.values(), axis = 0)
	lipids_ff_contacts_u2l_during_by_size_nb_sum = lipids_ff_contacts_u2l_during_by_size_nb_sum / float(np.sum(lipids_ff_contacts_u2l_during_by_size_nb_sum)) * 100
	lipids_ff_contacts_u2l_outside_by_size_nb_sum = lipids_ff_contacts_u2l_outside_by_size_nb_sum / float(np.sum(lipids_ff_contacts_u2l_outside_by_size_nb_sum)) * 100	
	global lipids_ff_contacts_l2u_during_by_size_nb_sum, lipids_ff_contacts_l2u_outside_by_size_nb_sum
	lipids_ff_contacts_l2u_during_by_size_nb_sum = np.sum(lipids_ff_contacts_l2u_during_tot_nb_by_size.values(), axis = 0)
	lipids_ff_contacts_l2u_outside_by_size_nb_sum = np.sum(lipids_ff_contacts_l2u_outside_tot_nb_by_size.values(), axis = 0)
	lipids_ff_contacts_l2u_during_by_size_nb_sum = lipids_ff_contacts_l2u_during_by_size_nb_sum / float(np.sum(lipids_ff_contacts_l2u_during_by_size_nb_sum)) * 100
	lipids_ff_contacts_l2u_outside_by_size_nb_sum = lipids_ff_contacts_l2u_outside_by_size_nb_sum / float(np.sum(lipids_ff_contacts_l2u_outside_by_size_nb_sum)) * 100

	#distribution of contacts over groups
	if args.cluster_groups_file != "no":
		global lipids_ff_contacts_u2l_during_by_size_group_avg, lipids_ff_contacts_u2l_during_by_size_group_std
		global lipids_ff_contacts_l2u_during_by_size_group_avg, lipids_ff_contacts_l2u_during_by_size_group_std
		global lipids_ff_contacts_u2l_outside_by_size_group_avg, lipids_ff_contacts_u2l_outside_by_size_group_std
		global lipids_ff_contacts_l2u_outside_by_size_group_avg, lipids_ff_contacts_l2u_outside_by_size_group_std
		lipids_ff_contacts_u2l_during_by_size_group_avg = np.average(lipids_ff_contacts_u2l_during_tot_pc_by_size_group.values(), axis = 0)
		lipids_ff_contacts_u2l_during_by_size_group_std = np.std(lipids_ff_contacts_u2l_during_tot_pc_by_size_group.values(), axis = 0)
		lipids_ff_contacts_u2l_outside_by_size_group_avg = np.average(lipids_ff_contacts_u2l_outside_tot_pc_by_size_group.values(), axis = 0)
		lipids_ff_contacts_u2l_outside_by_size_group_std = np.std(lipids_ff_contacts_u2l_outside_tot_pc_by_size_group.values(), axis = 0)
		lipids_ff_contacts_l2u_during_by_size_group_avg = np.average(lipids_ff_contacts_l2u_during_tot_pc_by_size_group.values(), axis = 0)
		lipids_ff_contacts_l2u_during_by_size_group_std = np.std(lipids_ff_contacts_l2u_during_tot_pc_by_size_group.values(), axis = 0)
		lipids_ff_contacts_l2u_outside_by_size_group_avg = np.average(lipids_ff_contacts_l2u_outside_tot_pc_by_size_group.values(), axis = 0)
		lipids_ff_contacts_l2u_outside_by_size_group_std = np.std(lipids_ff_contacts_l2u_outside_tot_pc_by_size_group.values(), axis = 0)

		#new
		global lipids_ff_contacts_u2l_during_by_size_group_nb_sum, lipids_ff_contacts_u2l_outside_by_size_group_nb_sum
		lipids_ff_contacts_u2l_during_by_size_group_nb_sum = np.sum(lipids_ff_contacts_u2l_during_tot_nb_by_size_group.values(), axis = 0)
		lipids_ff_contacts_u2l_outside_by_size_group_nb_sum = np.sum(lipids_ff_contacts_u2l_outside_tot_nb_by_size_group.values(), axis = 0)
		lipids_ff_contacts_u2l_during_by_size_group_nb_sum = lipids_ff_contacts_u2l_during_by_size_group_nb_sum / float(np.sum(lipids_ff_contacts_u2l_during_by_size_group_nb_sum)) * 100
		lipids_ff_contacts_u2l_outside_by_size_group_nb_sum = lipids_ff_contacts_u2l_outside_by_size_group_nb_sum / float(np.sum(lipids_ff_contacts_u2l_outside_by_size_group_nb_sum)) * 100
		global lipids_ff_contacts_l2u_during_by_size_group_nb_sum, lipids_ff_contacts_l2u_outside_by_size_group_nb_sum
		lipids_ff_contacts_l2u_during_by_size_group_nb_sum = np.sum(lipids_ff_contacts_l2u_during_tot_nb_by_size_group.values(), axis = 0)
		lipids_ff_contacts_l2u_outside_by_size_group_nb_sum = np.sum(lipids_ff_contacts_l2u_outside_tot_nb_by_size_group.values(), axis = 0)
		lipids_ff_contacts_l2u_during_by_size_group_nb_sum = lipids_ff_contacts_l2u_during_by_size_group_nb_sum / float(np.sum(lipids_ff_contacts_l2u_during_by_size_group_nb_sum)) * 100
		lipids_ff_contacts_l2u_outside_by_size_group_nb_sum = lipids_ff_contacts_l2u_outside_by_size_group_nb_sum / float(np.sum(lipids_ff_contacts_l2u_outside_by_size_group_nb_sum)) * 100

	#distribution of contacts over residue types
	global lipids_ff_contacts_u2l_during_by_type_avg, lipids_ff_contacts_u2l_during_by_type_std
	global lipids_ff_contacts_l2u_during_by_type_avg, lipids_ff_contacts_l2u_during_by_type_std
	global lipids_ff_contacts_u2l_outside_by_type_avg, lipids_ff_contacts_u2l_outside_by_type_std
	global lipids_ff_contacts_l2u_outside_by_type_avg, lipids_ff_contacts_l2u_outside_by_type_std
	lipids_ff_contacts_u2l_during_by_type_avg = np.average(lipids_ff_contacts_u2l_during_tot_pc_by_type.values(), axis = 0)
	lipids_ff_contacts_u2l_during_by_type_std = np.std(lipids_ff_contacts_u2l_during_tot_pc_by_type.values(), axis = 0)
	lipids_ff_contacts_u2l_outside_by_type_avg = np.average(lipids_ff_contacts_u2l_outside_tot_pc_by_type.values(), axis = 0)
	lipids_ff_contacts_u2l_outside_by_type_std = np.std(lipids_ff_contacts_u2l_outside_tot_pc_by_type.values(), axis = 0)
	lipids_ff_contacts_l2u_during_by_type_avg = np.average(lipids_ff_contacts_l2u_during_tot_pc_by_type.values(), axis = 0)
	lipids_ff_contacts_l2u_during_by_type_std = np.std(lipids_ff_contacts_l2u_during_tot_pc_by_type.values(), axis = 0)
	lipids_ff_contacts_l2u_outside_by_type_avg = np.average(lipids_ff_contacts_l2u_outside_tot_pc_by_type.values(), axis = 0)
	lipids_ff_contacts_l2u_outside_by_type_std = np.std(lipids_ff_contacts_l2u_outside_tot_pc_by_type.values(), axis = 0)
	
	#new
	global lipids_ff_contacts_u2l_during_by_type_nb_sum, lipids_ff_contacts_u2l_outside_by_type_nb_sum
	lipids_ff_contacts_u2l_during_by_type_nb_sum = np.sum(lipids_ff_contacts_u2l_during_tot_nb_by_type.values(), axis = 0)
	lipids_ff_contacts_u2l_outside_by_type_nb_sum = np.sum(lipids_ff_contacts_u2l_outside_tot_nb_by_type.values(), axis = 0)
	lipids_ff_contacts_u2l_during_by_type_nb_sum = lipids_ff_contacts_u2l_during_by_type_nb_sum / float(np.sum(lipids_ff_contacts_u2l_during_by_type_nb_sum)) * 100
	lipids_ff_contacts_u2l_outside_by_type_nb_sum = lipids_ff_contacts_u2l_outside_by_type_nb_sum / float(np.sum(lipids_ff_contacts_u2l_outside_by_type_nb_sum)) * 100
	global lipids_ff_contacts_l2u_during_by_type_nb_sum, lipids_ff_contacts_l2u_outside_by_type_nb_sum
	lipids_ff_contacts_l2u_during_by_type_nb_sum = np.sum(lipids_ff_contacts_l2u_during_tot_nb_by_type.values(), axis = 0)
	lipids_ff_contacts_l2u_outside_by_type_nb_sum = np.sum(lipids_ff_contacts_l2u_outside_tot_nb_by_type.values(), axis = 0)
	lipids_ff_contacts_l2u_during_by_type_nb_sum = lipids_ff_contacts_l2u_during_by_type_nb_sum / float(np.sum(lipids_ff_contacts_l2u_during_by_type_nb_sum)) * 100
	lipids_ff_contacts_l2u_outside_by_type_nb_sum = lipids_ff_contacts_l2u_outside_by_type_nb_sum / float(np.sum(lipids_ff_contacts_l2u_outside_by_type_nb_sum)) * 100

	#distribution of contacts over normal
	if args.profile:
		global lipids_ff_contacts_u2l_during_profile_avg, lipids_ff_contacts_u2l_during_profile_std
		global lipids_ff_contacts_l2u_during_profile_avg, lipids_ff_contacts_l2u_during_profile_std
		global lipids_ff_contacts_u2l_outside_profile_avg, lipids_ff_contacts_u2l_outside_profile_std
		global lipids_ff_contacts_l2u_outside_profile_avg, lipids_ff_contacts_l2u_outside_profile_std
		lipids_ff_contacts_u2l_during_profile_avg = np.average(tmp_profile_u2l_during_sizes.values(), axis = 0)
		lipids_ff_contacts_u2l_during_profile_std = np.std(tmp_profile_u2l_during_sizes.values(), axis = 0)
		lipids_ff_contacts_u2l_outside_profile_avg = np.average(tmp_profile_u2l_outside_sizes.values(), axis = 0)
		lipids_ff_contacts_u2l_outside_profile_std = np.std(tmp_profile_u2l_outside_sizes.values(), axis = 0)
		lipids_ff_contacts_l2u_during_profile_avg = np.average(tmp_profile_l2u_during_sizes.values(), axis = 0)
		lipids_ff_contacts_l2u_during_profile_std = np.std(tmp_profile_l2u_during_sizes.values(), axis = 0)
		lipids_ff_contacts_l2u_outside_profile_avg = np.average(tmp_profile_l2u_outside_sizes.values(), axis = 0)
		lipids_ff_contacts_l2u_outside_profile_std = np.std(tmp_profile_l2u_outside_sizes.values(), axis = 0)
		if args.cluster_groups_file != "no":	
			global lipids_ff_contacts_u2l_during_profile_groups_avg, lipids_ff_contacts_u2l_during_profile_groups_std
			global lipids_ff_contacts_l2u_during_profile_groups_avg, lipids_ff_contacts_l2u_during_profile_groups_std
			global lipids_ff_contacts_u2l_outside_profile_groups_avg, lipids_ff_contacts_u2l_outside_profile_groups_std
			global lipids_ff_contacts_l2u_outside_profile_groups_avg, lipids_ff_contacts_l2u_outside_profile_groups_std
			lipids_ff_contacts_u2l_during_profile_groups_avg = {}
			lipids_ff_contacts_u2l_during_profile_groups_std = {}
			lipids_ff_contacts_l2u_during_profile_groups_avg = {}
			lipids_ff_contacts_l2u_during_profile_groups_std = {}
			lipids_ff_contacts_u2l_outside_profile_groups_avg = {}
			lipids_ff_contacts_u2l_outside_profile_groups_std = {}			
			lipids_ff_contacts_l2u_outside_profile_groups_avg = {}
			lipids_ff_contacts_l2u_outside_profile_groups_std = {}
			for g_index in range(0,group_gmax):
				#u2l
				if len(tmp_profile_u2l_during_groups[g_index].values()) > 0:
					lipids_ff_contacts_u2l_during_profile_groups_avg[g_index] = np.average(tmp_profile_u2l_during_groups[g_index].values(), axis = 0)
					lipids_ff_contacts_u2l_during_profile_groups_std[g_index] = np.std(tmp_profile_u2l_during_groups[g_index].values(), axis = 0)
				else:
					lipids_ff_contacts_u2l_during_profile_groups_avg[g_index] = np.zeros((4,2*bins_nb))
					lipids_ff_contacts_u2l_during_profile_groups_std[g_index] = np.zeros((4,2*bins_nb))
				if len(tmp_profile_u2l_outside_groups[g_index].values()) > 0:
					lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index] = np.average(tmp_profile_u2l_outside_groups[g_index].values(), axis = 0)
					lipids_ff_contacts_u2l_outside_profile_groups_std[g_index] = np.std(tmp_profile_u2l_outside_groups[g_index].values(), axis = 0)
				else:
					lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index] = np.zeros((4,2*bins_nb))
					lipids_ff_contacts_u2l_outside_profile_groups_std[g_index] = np.zeros((4,2*bins_nb))
				#l2u
				if len(tmp_profile_l2u_during_groups[g_index].values()) > 0:
					lipids_ff_contacts_l2u_during_profile_groups_avg[g_index] = np.average(tmp_profile_l2u_during_groups[g_index].values(), axis = 0)
					lipids_ff_contacts_l2u_during_profile_groups_std[g_index] = np.std(tmp_profile_l2u_during_groups[g_index].values(), axis = 0)
				else:
					lipids_ff_contacts_l2u_during_profile_groups_avg[g_index] = np.zeros((4,2*bins_nb))
					lipids_ff_contacts_l2u_during_profile_groups_std[g_index] = np.zeros((4,2*bins_nb))
				if len(tmp_profile_l2u_outside_groups[g_index].values()) > 0:
					lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index] = np.average(tmp_profile_l2u_outside_groups[g_index].values(), axis = 0)
					lipids_ff_contacts_l2u_outside_profile_groups_std[g_index] = np.std(tmp_profile_l2u_outside_groups[g_index].values(), axis = 0)
				else:
					lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index] = np.zeros((4,2*bins_nb))
					lipids_ff_contacts_l2u_outside_profile_groups_std[g_index] = np.zeros((4,2*bins_nb))
			
	return

#=========================================================================================
# outputs
#=========================================================================================

#contacts distribution over residue types
def write_ff_ctcts_by_type():
	
	#averages
	#========
	filename=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_pc.stat'
	output_stat = open(filename, 'w')	
	output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if args.m_algorithm == "density":
		output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
		output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
		output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
	elif args.m_algorithm == "min":
		output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	else:
		output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")

	#caption
	output_stat.write("\n")
	output_stat.write("caption: average distribution of contacts over residue types (%)\n")
	
	#upper to lower
	if np.size(lipids_ff_u2l_index) > 0:
		output_stat.write("\n")
		output_stat.write("upper to lower (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		output_stat.write("		basic	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		output_stat.write("peptide (ref)	"  + str(round(protein_composition[0],1)) + "	" + str(round(protein_composition[1],1)) + "	" + str(round(protein_composition[2],1)) + "	" + str(round(protein_composition[3],1)) + "	" + str(round(np.sum(protein_composition),1)) + "\n")
		output_stat.write("during		" + str(round(lipids_ff_contacts_u2l_during_by_type_avg[0],1)) + "	" + str(round(lipids_ff_contacts_u2l_during_by_type_avg[1],1)) + " 	" + str(round(lipids_ff_contacts_u2l_during_by_type_avg[2],1)) + "	" + str(round(lipids_ff_contacts_u2l_during_by_type_avg[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_u2l_during_by_type_avg),1)) + "\n")
		output_stat.write("before/after	" + str(round(lipids_ff_contacts_u2l_outside_by_type_avg[0],1)) + "	" + str(round(lipids_ff_contacts_u2l_outside_by_type_avg[1],1)) + " 	" + str(round(lipids_ff_contacts_u2l_outside_by_type_avg[2],1)) + "	" + str(round(lipids_ff_contacts_u2l_outside_by_type_avg[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_u2l_outside_by_type_avg),1)) + "\n")
		output_stat.write("\n")

	#lower to upper	
	if np.size(lipids_ff_l2u_index) > 0:
		output_stat.write("\n")
		output_stat.write("lower to upper (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		output_stat.write("		Q+	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		output_stat.write("peptide (ref)	"  + str(round(protein_composition[0],1)) + "	" + str(round(protein_composition[1],1)) + "	" + str(round(protein_composition[2],1)) + "	" + str(round(protein_composition[3],1)) + "	" + str(round(np.sum(protein_composition),1)) + "\n")
		output_stat.write("during		" + str(round(lipids_ff_contacts_l2u_during_by_type_avg[0],1)) + "	" + str(round(lipids_ff_contacts_l2u_during_by_type_avg[1],1)) + " 	" + str(round(lipids_ff_contacts_l2u_during_by_type_avg[2],1)) + "	" + str(round(lipids_ff_contacts_l2u_during_by_type_avg[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_l2u_during_by_type_avg),1)) + "\n")
		output_stat.write("before/after	" + str(round(lipids_ff_contacts_l2u_outside_by_type_avg[0],1)) + "	" + str(round(lipids_ff_contacts_l2u_outside_by_type_avg[1],1)) + " 	" + str(round(lipids_ff_contacts_l2u_outside_by_type_avg[2],1)) + "	" + str(round(lipids_ff_contacts_l2u_outside_by_type_avg[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_l2u_outside_by_type_avg),1)) + "\n")
		output_stat.write("\n")
	output_stat.close()
	
	#details: u2l
	#============
	if np.size(lipids_ff_u2l_index) > 0:
		filename=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_pc_u2l.stat'
		output_stat = open(filename, 'w')	
		output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_stat.write("\n")
	
		#general info
		output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#caption
		output_stat.write("\n")
		output_stat.write("caption: total nb of contacts with each residue types for lipids flip-flopping from the upper to the lower leaflet (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")
	
		#before/after
		output_stat.write("\n")
		output_stat.write("before/after\n")	
		output_stat.write("============\n")
		output_stat.write("		basic	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		for l_index in lipids_ff_u2l_index:
			output_stat.write(str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"  + str(lipids_ff_contacts_u2l_outside_tot_nb_by_type[l_index][0]) + "	" + str(lipids_ff_contacts_u2l_outside_tot_nb_by_type[l_index][1]) + "	" + str(lipids_ff_contacts_u2l_outside_tot_nb_by_type[l_index][2]) + "	" + str(lipids_ff_contacts_u2l_outside_tot_nb_by_type[l_index][3]) + "	" + str(np.sum(lipids_ff_contacts_u2l_outside_tot_nb_by_type[l_index])) + "\n")
		output_stat.write("\n")

		#during
		output_stat.write("\n")
		output_stat.write("during\n")	
		output_stat.write("======\n")
		output_stat.write("		basic	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		for l_index in lipids_ff_u2l_index:
			output_stat.write(str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"  + str(lipids_ff_contacts_u2l_during_tot_nb_by_type[l_index][0]) + "	" + str(lipids_ff_contacts_u2l_during_tot_nb_by_type[l_index][1]) + "	" + str(lipids_ff_contacts_u2l_during_tot_nb_by_type[l_index][2]) + "	" + str(lipids_ff_contacts_u2l_during_tot_nb_by_type[l_index][3]) + "	" + str(np.sum(lipids_ff_contacts_u2l_during_tot_nb_by_type[l_index])) + "\n")
		output_stat.write("\n")
		output_stat.close()
	
	#details: l2u
	#============
	if np.size(lipids_ff_l2u_index) > 0:
		filename=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_pc_l2u.stat'
		output_stat = open(filename, 'w')	
		output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_stat.write("\n")
	
		#general info
		output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#caption
		output_stat.write("\n")
		output_stat.write("caption: total nb of contacts with each residue types for lipids flip-flopping from the lower to the upper leaflet (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")
	
		#before/after
		output_stat.write("\n")
		output_stat.write("before/after\n")	
		output_stat.write("============\n")
		output_stat.write("		basic	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		for l_index in lipids_ff_l2u_index:
			output_stat.write(str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"  + str(lipids_ff_contacts_l2u_outside_tot_nb_by_type[l_index][0]) + "	" + str(lipids_ff_contacts_l2u_outside_tot_nb_by_type[l_index][1]) + "	" + str(lipids_ff_contacts_l2u_outside_tot_nb_by_type[l_index][2]) + "	" + str(lipids_ff_contacts_l2u_outside_tot_nb_by_type[l_index][3]) + "	" + str(np.sum(lipids_ff_contacts_l2u_outside_tot_nb_by_type[l_index])) + "\n")
		output_stat.write("\n")

		#during
		output_stat.write("\n")
		output_stat.write("during\n")	
		output_stat.write("======\n")
		output_stat.write("		basic	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		for l_index in lipids_ff_l2u_index:
			output_stat.write(str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"  + str(lipids_ff_contacts_l2u_during_tot_nb_by_type[l_index][0]) + "	" + str(lipids_ff_contacts_l2u_during_tot_nb_by_type[l_index][1]) + "	" + str(lipids_ff_contacts_l2u_during_tot_nb_by_type[l_index][2]) + "	" + str(lipids_ff_contacts_l2u_during_tot_nb_by_type[l_index][3]) + "	" + str(np.sum(lipids_ff_contacts_l2u_during_tot_nb_by_type[l_index])) + "\n")
		output_stat.write("\n")
		output_stat.close()

	return
def graph_ff_ctcts_by_type():

	#-------------------------------------------------------------------
	#-what: distribution of ff contacts over residue types 
	#-plot: 2 bar charts (ff u2l and l2u) each with 3 bars (peptide, before/after, during) for each residue type
	#-------------------------------------------------------------------
			
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_pc.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_pc.svg'
	
	#create figure
	#-------------
	fig=plt.figure(figsize=(6, 5)) 								#1 column format
	#fig.suptitle("Distribution of peptides orientation")
	xticks_pos=np.arange(1,5)
	xticks_lab=['basic', 'polar', 'hydrophobic', 'backbone']

	#plot data: upper to lower
	#-------------------------
	ax1 = fig.add_subplot(211)
	#peptide
	plt.bar(xticks_pos-0.375, protein_composition, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_u2l_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_u2l_outside_by_type_avg, width=0.25, color='#808080', label="before/after", yerr = lipids_ff_contacts_u2l_outside_by_type_std, ecolor='k')
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_u2l_during_by_type_avg, width=0.25, color='w', label="during", yerr = lipids_ff_contacts_u2l_during_by_type_std, ecolor='k')
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Upper to Lower", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_u2l_index) > 0:
		ax1.legend(prop=fontP)

	#plot data: lower to upper
	#-------------------------
	ax2 = fig.add_subplot(212)
	#peptide
	plt.bar(xticks_pos-0.375, protein_composition, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_l2u_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_l2u_outside_by_type_avg, width=0.25, color='#808080', label="before/after", yerr = lipids_ff_contacts_l2u_outside_by_type_std, ecolor='k')
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_l2u_during_by_type_avg, width=0.25, color='w', label="during", yerr = lipids_ff_contacts_l2u_during_by_type_std, ecolor='k')
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Lower to Upper", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_l2u_index) > 0:
		ax2.legend(prop=fontP)

	#legend, labels, etc
	#-------------------
	#ax1.set_ylabel('nb of peptides', fontsize="small")
	#ax2.set_ylabel('% of peptides', fontsize="small")
	ax1.tick_params(axis='both', direction='out')
	ax2.tick_params(axis='both', direction='out')
	ax1.spines["right"].set_visible(False)								# remove unneeded axes
	ax1.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["top"].set_visible(False)
	ax1.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax1.get_yaxis().tick_left()
	ax2.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax2.get_yaxis().tick_left()
	ax1.set_xlim(0.5, 4.5)
	ax2.set_xlim(0.5, 4.5)
	ax1.set_ylim(ymin=0)
	ax2.set_ylim(ymin=0)
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="x-small")
	#plt.setp(ax1.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")	
	#plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")

	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax1.yaxis.get_major_formatter().set_powerlimits((-3, 3))

	#save figure
	#-----------
	plt.subplots_adjust(left=0.1, right=0.96, hspace=0.3, top=0.92, bottom=0.07)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def write_ff_ctcts_by_type_sum():
	
	#averages
	#========
	filename=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_sum.stat'
	output_stat = open(filename, 'w')	
	output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if args.m_algorithm == "density":
		output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
		output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
		output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
	elif args.m_algorithm == "min":
		output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	else:
		output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")

	#caption
	output_stat.write("\n")
	output_stat.write("caption: average distribution of contacts over residue types (%)\n")
	
	#upper to lower
	if np.size(lipids_ff_u2l_index) > 0:
		output_stat.write("\n")
		output_stat.write("upper to lower (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		output_stat.write("		basic	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		output_stat.write("peptide (ref)	"  + str(round(protein_composition[0],1)) + "	" + str(round(protein_composition[1],1)) + "	" + str(round(protein_composition[2],1)) + "	" + str(round(protein_composition[3],1)) + "	" + str(round(np.sum(protein_composition),1)) + "\n")
		output_stat.write("during		" + str(round(lipids_ff_contacts_u2l_during_by_type_nb_sum[0],1)) + "	" + str(round(lipids_ff_contacts_u2l_during_by_type_nb_sum[1],1)) + " 	" + str(round(lipids_ff_contacts_u2l_during_by_type_nb_sum[2],1)) + "	" + str(round(lipids_ff_contacts_u2l_during_by_type_nb_sum[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_u2l_during_by_type_nb_sum),1)) + "\n")
		output_stat.write("before/after	" + str(round(lipids_ff_contacts_u2l_outside_by_type_nb_sum[0],1)) + "	" + str(round(lipids_ff_contacts_u2l_outside_by_type_nb_sum[1],1)) + " 	" + str(round(lipids_ff_contacts_u2l_outside_by_type_nb_sum[2],1)) + "	" + str(round(lipids_ff_contacts_u2l_outside_by_type_nb_sum[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_u2l_outside_by_type_nb_sum),1)) + "\n")
		output_stat.write("\n")

	#lower to upper	
	if np.size(lipids_ff_l2u_index) > 0:
		output_stat.write("\n")
		output_stat.write("lower to upper (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		output_stat.write("		Q+	polar	hphobic	BB	total\n")
		output_stat.write("-----------------------------------------------------\n")
		output_stat.write("peptide (ref)	"  + str(round(protein_composition[0],1)) + "	" + str(round(protein_composition[1],1)) + "	" + str(round(protein_composition[2],1)) + "	" + str(round(protein_composition[3],1)) + "	" + str(round(np.sum(protein_composition),1)) + "\n")
		output_stat.write("during		" + str(round(lipids_ff_contacts_l2u_during_by_type_nb_sum[0],1)) + "	" + str(round(lipids_ff_contacts_l2u_during_by_type_nb_sum[1],1)) + " 	" + str(round(lipids_ff_contacts_l2u_during_by_type_nb_sum[2],1)) + "	" + str(round(lipids_ff_contacts_l2u_during_by_type_nb_sum[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_l2u_during_by_type_nb_sum),1)) + "\n")
		output_stat.write("before/after	" + str(round(lipids_ff_contacts_l2u_outside_by_type_nb_sum[0],1)) + "	" + str(round(lipids_ff_contacts_l2u_outside_by_type_nb_sum[1],1)) + " 	" + str(round(lipids_ff_contacts_l2u_outside_by_type_nb_sum[2],1)) + "	" + str(round(lipids_ff_contacts_l2u_outside_by_type_nb_sum[3],1)) + "	" + str(round(np.sum(lipids_ff_contacts_l2u_outside_by_type_nb_sum),1)) + "\n")
		output_stat.write("\n")
	output_stat.close()
	
	return
def graph_ff_ctcts_by_type_sum():

	#-------------------------------------------------------------------
	#-what: distribution of ff contacts over residue types 
	#-plot: 2 bar charts (ff u2l and l2u) each with 3 bars (peptide, before/after, during) for each residue type
	#-------------------------------------------------------------------
			
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_sum.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/by_type/ff_ctcts_by_type_sum.svg'
	
	#create figure
	#-------------
	fig=plt.figure(figsize=(6, 5)) 								#1 column format
	#fig.suptitle("Distribution of peptides orientation")
	xticks_pos=np.arange(1,5)
	xticks_lab=['basic', 'polar', 'hydrophobic', 'backbone']

	#plot data: upper to lower
	#-------------------------
	ax1 = fig.add_subplot(211)
	#peptide
	plt.bar(xticks_pos-0.375, protein_composition, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_u2l_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_u2l_outside_by_type_nb_sum, width=0.25, color='#808080', label="before/after")
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_u2l_during_by_type_nb_sum, width=0.25, color='w', label="during")
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Upper to Lower", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_u2l_index) > 0:
		ax1.legend(prop=fontP)

	#plot data: lower to upper
	#-------------------------
	ax2 = fig.add_subplot(212)
	#peptide
	plt.bar(xticks_pos-0.375, protein_composition, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_l2u_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_l2u_outside_by_type_nb_sum, width=0.25, color='#808080', label="before/after")
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_l2u_during_by_type_nb_sum, width=0.25, color='w', label="during")
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Lower to Upper", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_l2u_index) > 0:
		ax2.legend(prop=fontP)

	#legend, labels, etc
	#-------------------
	#ax1.set_ylabel('nb of peptides', fontsize="small")
	#ax2.set_ylabel('% of peptides', fontsize="small")
	ax1.tick_params(axis='both', direction='out')
	ax2.tick_params(axis='both', direction='out')
	ax1.spines["right"].set_visible(False)								# remove unneeded axes
	ax1.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["top"].set_visible(False)
	ax1.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax1.get_yaxis().tick_left()
	ax2.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax2.get_yaxis().tick_left()
	ax1.set_xlim(0.5, 4.5)
	ax2.set_xlim(0.5, 4.5)
	ax1.set_ylim(ymin=0)
	ax2.set_ylim(ymin=0)
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="x-small")
	#plt.setp(ax1.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")	
	#plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")

	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax1.yaxis.get_major_formatter().set_powerlimits((-3, 3))

	#save figure
	#-----------
	plt.subplots_adjust(left=0.1, right=0.96, hspace=0.3, top=0.92, bottom=0.07)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return

#contacts distribution over cluster sizes
def write_ff_ctcts_by_size():
	
	#averages
	#========
	filename=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_pc.stat'
	output_stat = open(filename, 'w')	
	output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if args.m_algorithm == "density":
		output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
		output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
		output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
	elif args.m_algorithm == "min":
		output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	else:
		output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")

	#caption
	output_stat.write("\n")
	output_stat.write("caption: average distribution of contacts over cluster sizes (%)\n")
	
	#upper to lower
	if np.size(lipids_ff_u2l_index) > 0:
		output_stat.write("\n")
		output_stat.write("upper to lower (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_sizes[c],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_u2l_during_by_size_avg[c],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_u2l_outside_by_size_avg[c],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_sizes),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_u2l_during_by_size_avg),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_u2l_outside_by_size_avg),1)) + ")\n")
		output_stat.write("\n")

	#lower to upper	
	if np.size(lipids_ff_l2u_index) > 0:
		output_stat.write("\n")
		output_stat.write("lower to upper (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_sizes[c],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_l2u_during_by_size_avg[c],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_l2u_outside_by_size_avg[c],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_sizes),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_l2u_during_by_size_avg),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_l2u_outside_by_size_avg),1)) + ")\n")
		output_stat.write("\n")
	output_stat.close()

	#details: u2l
	#============
	if np.size(lipids_ff_u2l_index) > 0:
		filename=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_pc_u2l.stat'
		output_stat = open(filename, 'w')	
		output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_stat.write("\n")
	
		#general info
		output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#caption
		output_stat.write("\n")
		output_stat.write("caption: total nb of contacts with each cluster size for lipids flip-flopping from the upper to the lower leaflet (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")
	
		#before/after
		output_stat.write("\n")
		output_stat.write("before/after\n")	
		output_stat.write("============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_u2l_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for c in range(0, protein_max_size_sampled):
				tmp_lip += str(lipids_ff_contacts_u2l_outside_tot_nb_by_size[l_index][c]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_u2l_outside_tot_nb_by_size[l_index])) + ")\n")
		output_stat.write("\n")

		#during
		output_stat.write("\n")
		output_stat.write("during\n")	
		output_stat.write("======\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_u2l_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for c in range(0, protein_max_size_sampled):
				tmp_lip += str(lipids_ff_contacts_u2l_during_tot_nb_by_size[l_index][c]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_u2l_during_tot_nb_by_size[l_index])) + ")\n")
		output_stat.write("\n")
		output_stat.close()
	
	#details: l2u
	#============
	if np.size(lipids_ff_l2u_index) > 0:
		filename=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_pc_l2u.stat'
		output_stat = open(filename, 'w')	
		output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_stat.write("\n")
	
		#general info
		output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#caption
		output_stat.write("\n")
		output_stat.write("caption: total nb of contacts with each cluster size for lipids flip-flopping from the lower to the upper leaflet (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")
	
		#before/after
		output_stat.write("\n")
		output_stat.write("before/after\n")	
		output_stat.write("============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_l2u_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for c in range(0, protein_max_size_sampled):
				tmp_lip += str(lipids_ff_contacts_l2u_outside_tot_nb_by_size[l_index][c]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_l2u_outside_tot_nb_by_size[l_index])) + ")\n")
		output_stat.write("\n")

		#during
		output_stat.write("\n")
		output_stat.write("during\n")	
		output_stat.write("======\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_l2u_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for c in range(0, protein_max_size_sampled):
				tmp_lip += str(lipids_ff_contacts_l2u_during_tot_nb_by_size[l_index][c]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_l2u_during_tot_nb_by_size[l_index])) + ")\n")
		output_stat.write("\n")
		output_stat.close()

	return
def graph_ff_ctcts_by_size():

	#-------------------------------------------------------------------
	#-what: distribution of ff contacts over cluster size
	#-plot: 2 bar charts (ff u2l and l2u) each with 3 bars (peptide, before/after, during) for each cluster size
	#-------------------------------------------------------------------
			
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_pc.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_pc.svg'
		
	#create figure
	#-------------
	fig=plt.figure(figsize=(6, 5)) 								#1 column format
	#fig.suptitle("Distribution of peptides orientation")
	xticks_pos=np.arange(1,protein_max_size_sampled + 1)
	xticks_lab=[str(c_size) for c_size in range(1, protein_max_size_sampled+1)]

	#plot data: upper to lower
	#-------------------------
	ax1 = fig.add_subplot(211)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_sizes, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_u2l_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_u2l_outside_by_size_avg[:protein_max_size_sampled], width=0.25, color='#808080', label="before/after", yerr = lipids_ff_contacts_u2l_outside_by_size_std[:protein_max_size_sampled], ecolor='k')
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_u2l_during_by_size_avg[:protein_max_size_sampled], width=0.25, color='w', label="during", yerr = lipids_ff_contacts_u2l_during_by_size_std[:protein_max_size_sampled], ecolor='k')
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Upper to Lower", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_u2l_index) > 0:
		ax1.legend(prop=fontP)

	#plot data: lower to upper
	#-------------------------
	ax2 = fig.add_subplot(212)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_sizes, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_l2u_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_l2u_outside_by_size_avg[:protein_max_size_sampled], width=0.25, color='#808080', label="before/after", yerr = lipids_ff_contacts_l2u_outside_by_size_std[:protein_max_size_sampled], ecolor='k')
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_l2u_during_by_size_avg[:protein_max_size_sampled], width=0.25, color='w', label="during", yerr = lipids_ff_contacts_l2u_during_by_size_std[:protein_max_size_sampled], ecolor='k')
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Lower to Upper", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_l2u_index) > 0:
		ax2.legend(prop=fontP)

	#legend, labels, etc
	#-------------------
	#ax1.set_ylabel('nb of peptides', fontsize="small")
	#ax2.set_ylabel('% of peptides', fontsize="small")
	ax1.tick_params(axis='both', direction='out')
	ax2.tick_params(axis='both', direction='out')
	ax1.spines["right"].set_visible(False)								# remove unneeded axes
	ax1.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["top"].set_visible(False)
	ax1.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax1.get_yaxis().tick_left()
	ax2.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax2.get_yaxis().tick_left()
	ax1.set_xlim(0.5, protein_max_size_sampled + 0.5)
	ax2.set_xlim(0.5, protein_max_size_sampled + 0.5)
	ax1.set_ylim(ymin=0)
	ax2.set_ylim(ymin=0)
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="x-small")
	#plt.setp(ax1.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")	
	#plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")

	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax1.yaxis.get_major_formatter().set_powerlimits((-3, 3))

	#save figure
	#-----------
	plt.subplots_adjust(left=0.1, right=0.96, hspace=0.3, top=0.92, bottom=0.07)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def write_ff_ctcts_by_size_sum():
	
	#averages
	#========
	filename=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_sum.stat'
	output_stat = open(filename, 'w')	
	output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if args.m_algorithm == "density":
		output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
		output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
		output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
	elif args.m_algorithm == "min":
		output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	else:
		output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")

	#caption
	output_stat.write("\n")
	output_stat.write("caption: average distribution of contacts over cluster sizes (%)\n")
	
	#upper to lower
	if np.size(lipids_ff_u2l_index) > 0:
		output_stat.write("\n")
		output_stat.write("upper to lower (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_sizes[c],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_u2l_during_by_size_nb_sum[c],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_u2l_outside_by_size_nb_sum[c],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_sizes),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_u2l_during_by_size_nb_sum),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_u2l_outside_by_size_nb_sum),1)) + ")\n")
		output_stat.write("\n")

	#lower to upper	
	if np.size(lipids_ff_l2u_index) > 0:
		output_stat.write("\n")
		output_stat.write("lower to upper (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for c in range(0, protein_max_size_sampled):
			tmp_title1 += str(c+1) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_sizes[c],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_l2u_during_by_size_nb_sum[c],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_l2u_outside_by_size_nb_sum[c],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_sizes),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_l2u_during_by_size_nb_sum),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_l2u_outside_by_size_nb_sum),1)) + ")\n")
		output_stat.write("\n")
	output_stat.close()

	return
def graph_ff_ctcts_by_size_sum():

	#-------------------------------------------------------------------
	#-what: distribution of ff contacts over cluster size
	#-plot: 2 bar charts (ff u2l and l2u) each with 3 bars (peptide, before/after, during) for each cluster size
	#-------------------------------------------------------------------
			
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_sum.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/by_size/ff_ctcts_by_size_sum.svg'
		
	#create figure
	#-------------
	fig=plt.figure(figsize=(6, 5)) 								#1 column format
	#fig.suptitle("Distribution of peptides orientation")
	xticks_pos=np.arange(1,protein_max_size_sampled + 1)
	xticks_lab=[str(c_size) for c_size in range(1, protein_max_size_sampled+1)]

	#plot data: upper to lower
	#-------------------------
	ax1 = fig.add_subplot(211)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_sizes, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_u2l_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_u2l_outside_by_size_nb_sum[:protein_max_size_sampled], width=0.25, color='#808080', label="before/after")
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_u2l_during_by_size_nb_sum[:protein_max_size_sampled], width=0.25, color='w', label="during")
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Upper to Lower", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_u2l_index) > 0:
		ax1.legend(prop=fontP)

	#plot data: lower to upper
	#-------------------------
	ax2 = fig.add_subplot(212)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_sizes, width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_l2u_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_l2u_outside_by_size_nb_sum[:protein_max_size_sampled], width=0.25, color='#808080', label="before/after")
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_l2u_during_by_size_nb_sum[:protein_max_size_sampled], width=0.25, color='w', label="during")
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Lower to Upper", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_l2u_index) > 0:
		ax2.legend(prop=fontP)

	#legend, labels, etc
	#-------------------
	#ax1.set_ylabel('nb of peptides', fontsize="small")
	#ax2.set_ylabel('% of peptides', fontsize="small")
	ax1.tick_params(axis='both', direction='out')
	ax2.tick_params(axis='both', direction='out')
	ax1.spines["right"].set_visible(False)								# remove unneeded axes
	ax1.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["top"].set_visible(False)
	ax1.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax1.get_yaxis().tick_left()
	ax2.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax2.get_yaxis().tick_left()
	ax1.set_xlim(0.5, protein_max_size_sampled + 0.5)
	ax2.set_xlim(0.5, protein_max_size_sampled + 0.5)
	ax1.set_ylim(ymin=0)
	ax2.set_ylim(ymin=0)
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="x-small")
	#plt.setp(ax1.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")	
	#plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")

	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax1.yaxis.get_major_formatter().set_powerlimits((-3, 3))

	#save figure
	#-----------
	plt.subplots_adjust(left=0.1, right=0.96, hspace=0.3, top=0.92, bottom=0.07)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return

#contacts distribution over cluster groups
def write_ff_ctcts_by_group():
	
	#averages
	#========
	filename=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_pc.stat'
	output_stat = open(filename, 'w')	
	output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if args.m_algorithm == "density":
		output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
		output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
		output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
	elif args.m_algorithm == "min":
		output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	else:
		output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")

	#caption
	output_stat.write("\n")
	output_stat.write("caption: average distribution of contacts over cluster size groups (%)\n")
	
	#upper to lower
	if np.size(lipids_ff_u2l_index) > 0:
		output_stat.write("\n")
		output_stat.write("upper to lower (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_groups[g],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_u2l_during_by_size_group_avg[g],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_u2l_outside_by_size_group_avg[g],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_groups),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_u2l_during_by_size_group_avg),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_u2l_outside_by_size_group_avg),1)) + ")\n")
		output_stat.write("\n")

	#lower to upper	
	if np.size(lipids_ff_l2u_index) > 0:
		output_stat.write("lower to upper (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		output_stat.write("\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_groups[g],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_l2u_during_by_size_group_avg[g],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_l2u_outside_by_size_group_avg[g],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_groups),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_l2u_during_by_size_group_avg),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_l2u_outside_by_size_group_avg),1)) + ")\n")
		output_stat.write("\n")
	output_stat.close()

	#details: u2l
	#============
	if np.size(lipids_ff_u2l_index) > 0:
		filename=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_pc_u2l.stat'
		output_stat = open(filename, 'w')	
		output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_stat.write("\n")
	
		#general info
		output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#caption
		output_stat.write("\n")
		output_stat.write("caption: total nb of contacts with each cluster size group for lipids flip-flopping from the upper to the lower leaflet (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")
	
		#before/after
		output_stat.write("\n")
		output_stat.write("before/after\n")	
		output_stat.write("============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_u2l_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for g in range(0, group_gmax):
				tmp_lip += str(lipids_ff_contacts_u2l_outside_tot_nb_by_size_group[l_index][g]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_u2l_outside_tot_nb_by_size_group[l_index])) + ")\n")
		output_stat.write("\n")

		#during
		output_stat.write("\n")
		output_stat.write("during\n")	
		output_stat.write("======\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_u2l_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for g in range(0, group_gmax):
				tmp_lip += str(lipids_ff_contacts_u2l_during_tot_nb_by_size_group[l_index][g]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_u2l_during_tot_nb_by_size_group[l_index])) + ")\n")
		output_stat.write("\n")
		output_stat.close()
	
	#details: l2u
	#============
	if np.size(lipids_ff_l2u_index) > 0:
		filename=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_pc_l2u.stat'
		output_stat = open(filename, 'w')	
		output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_stat.write("\n")
	
		#general info
		output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#caption
		output_stat.write("\n")
		output_stat.write("caption: total nb of contacts with each cluster size group for lipids flip-flopping from the lower to the upper leaflet (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")
	
		#before/after
		output_stat.write("\n")
		output_stat.write("before/after\n")	
		output_stat.write("============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_l2u_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for g in range(0, group_gmax):
				tmp_lip += str(lipids_ff_contacts_l2u_outside_tot_nb_by_size_group[l_index][g]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_l2u_outside_tot_nb_by_size_group[l_index])) + ")\n")
		output_stat.write("\n")

		#during
		output_stat.write("\n")
		output_stat.write("during\n")	
		output_stat.write("======\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		for l_index in lipids_ff_l2u_index:
			tmp_lip = str(lipids_ff_info[l_index][0]) + " " + str(str(lipids_ff_info[l_index][1])) + "	"	
			for g in range(0, group_gmax):
				tmp_lip += str(lipids_ff_contacts_l2u_during_tot_nb_by_size_group[l_index][g]) + "	"
			output_stat.write(tmp_lip + "(" + str(np.sum(lipids_ff_contacts_l2u_during_tot_nb_by_size_group[l_index])) + ")\n")
		output_stat.write("\n")
		output_stat.close()

	return
def graph_ff_ctcts_by_group():

	#-------------------------------------------------------------------
	#-what: distribution of ff contacts over cluster size group
	#-plot: 2 bar charts (ff u2l and l2u) each with 3 bars (peptide, before/after, during) for each cluster group
	#-------------------------------------------------------------------
			
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_pc.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_pc.svg'
		
	#create figure
	#-------------
	fig=plt.figure(figsize=(6, 5)) 								#1 column format
	#fig.suptitle("Distribution of peptides orientation")
	xticks_pos=np.arange(1,group_gmax+1)
	xticks_lab=[str(groups_labels[g]) for g in range(0, group_gmax)]

	#plot data: upper to lower
	#-------------------------
	ax1 = fig.add_subplot(211)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_groups[:group_gmax], width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_u2l_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_u2l_during_by_size_group_avg[:group_gmax], width=0.25, color='#808080', label="before/after", yerr = lipids_ff_contacts_u2l_during_by_size_group_std[:group_gmax], ecolor='k')
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_u2l_outside_by_size_group_avg[:group_gmax], width=0.25, color='w', label="during", yerr = lipids_ff_contacts_u2l_outside_by_size_group_std[:group_gmax], ecolor='k')
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Upper to Lower", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_u2l_index) > 0:
		ax1.legend(prop=fontP)

	#plot data: lower to upper
	#-------------------------
	ax2 = fig.add_subplot(212)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_groups[:group_gmax], width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_l2u_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_l2u_during_by_size_group_avg[:group_gmax], width=0.25, color='#808080', label="before/after", yerr = lipids_ff_contacts_l2u_during_by_size_group_std[:group_gmax], ecolor='k')
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_l2u_outside_by_size_group_avg[:group_gmax], width=0.25, color='w', label="during", yerr = lipids_ff_contacts_l2u_outside_by_size_group_std[:group_gmax], ecolor='k')
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Lower to Upper", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_l2u_index) > 0:
		ax2.legend(prop=fontP)

	#legend, labels, etc
	#-------------------
	#ax1.set_ylabel('nb of peptides', fontsize="small")
	#ax2.set_ylabel('% of peptides', fontsize="small")
	ax1.tick_params(axis='both', direction='out')
	ax2.tick_params(axis='both', direction='out')
	ax1.spines["right"].set_visible(False)								# remove unneeded axes
	ax1.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["top"].set_visible(False)
	ax1.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax1.get_yaxis().tick_left()
	ax2.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax2.get_yaxis().tick_left()
	ax1.set_xlim(0.5, group_gmax + 0.5)
	ax2.set_xlim(0.5, group_gmax + 0.5)
	ax1.set_ylim(ymin=0)
	ax2.set_ylim(ymin=0)
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="x-small")
	#plt.setp(ax1.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")	
	#plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")

	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax1.yaxis.get_major_formatter().set_powerlimits((-3, 3))

	#save figure
	#-----------
	plt.subplots_adjust(left=0.1, right=0.96, hspace=0.3, top=0.92, bottom=0.07)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def write_ff_ctcts_by_group_sum():
	
	#averages
	#========
	filename=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_sum.stat'
	output_stat = open(filename, 'w')	
	output_stat.write("[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("-nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if args.m_algorithm == "density":
		output_stat.write("-method cluster: density based algorithm using distances between proteins COGs\n")	
		output_stat.write(" -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
		output_stat.write(" -> nb neighbours: " + str(args.dbscan_nb) + "\n")
	elif args.m_algorithm == "min":
		output_stat.write("-method cluster: connectivity algorithm using minimum distance between proteins\n")
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	else:
		output_stat.write("-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
		output_stat.write(" -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
	output_stat.write("-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")

	#caption
	output_stat.write("\n")
	output_stat.write("caption: average distribution of contacts over cluster size groups (%)\n")
	
	#upper to lower
	if np.size(lipids_ff_u2l_index) > 0:
		output_stat.write("\n")
		output_stat.write("upper to lower (" + str(np.size(lipids_ff_u2l_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_groups[g],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_u2l_during_by_size_group_nb_sum[g],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_u2l_outside_by_size_group_nb_sum[g],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_groups),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_u2l_during_by_size_group_nb_sum),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_u2l_outside_by_size_group_nb_sum),1)) + ")\n")
		output_stat.write("\n")

	#lower to upper	
	if np.size(lipids_ff_l2u_index) > 0:
		output_stat.write("lower to upper (" + str(np.size(lipids_ff_l2u_index)) + " lipids)\n")	
		output_stat.write("==============\n")
		output_stat.write("\n")
		tmp_title1 = "		"
		tmp_title2 = "----------------"
		tmp_pep = "peptide (ref)	"
		tmp_dur = "during		"
		tmp_out = "before/after	"
		for g in range(0, group_gmax):
			tmp_title1 += str(groups_labels[g]) + "	"
			tmp_title2 += "--------"
			tmp_pep += str(round(protein_TM_distribution_groups[g],1)) + "	"
			tmp_dur += str(round(lipids_ff_contacts_l2u_during_by_size_group_nb_sum[g],1)) + "	"
			tmp_out += str(round(lipids_ff_contacts_l2u_outside_by_size_group_nb_sum[g],1)) + "	"
		output_stat.write(tmp_title1 + "\n")
		output_stat.write(tmp_title2 + "\n")
		output_stat.write(tmp_pep + "(" + str(round(np.sum(protein_TM_distribution_groups),1)) + ")\n")
		output_stat.write(tmp_dur + "(" + str(round(np.sum(lipids_ff_contacts_l2u_during_by_size_group_nb_sum),1)) + ")\n")
		output_stat.write(tmp_out + "(" + str(round(np.sum(lipids_ff_contacts_l2u_outside_by_size_group_nb_sum),1)) + ")\n")
		output_stat.write("\n")
	output_stat.close()

	return
def graph_ff_ctcts_by_group_sum():

	#-------------------------------------------------------------------
	#-what: distribution of ff contacts over cluster size group
	#-plot: 2 bar charts (ff u2l and l2u) each with 3 bars (peptide, before/after, during) for each cluster group
	#-------------------------------------------------------------------
			
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_sum.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/by_group/ff_ctcts_by_group_sum.svg'
		
	#create figure
	#-------------
	fig=plt.figure(figsize=(6, 5)) 								#1 column format
	#fig.suptitle("Distribution of peptides orientation")
	xticks_pos=np.arange(1,group_gmax+1)
	xticks_lab=[str(groups_labels[g]) for g in range(0, group_gmax)]

	#plot data: upper to lower
	#-------------------------
	ax1 = fig.add_subplot(211)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_groups[:group_gmax], width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_u2l_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_u2l_outside_by_size_group_nb_sum[:group_gmax], width=0.25, color='#808080', label="before/after")
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_u2l_during_by_size_group_nb_sum[:group_gmax], width=0.25, color='w', label="during")
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Upper to Lower", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_u2l_index) > 0:
		ax1.legend(prop=fontP)

	#plot data: lower to upper
	#-------------------------
	ax2 = fig.add_subplot(212)
	#peptide
	plt.bar(xticks_pos-0.375, protein_TM_distribution_groups[:group_gmax], width=0.25, color='#C0C0C0', label="peptide")	
	if np.size(lipids_ff_l2u_index) > 0:
		#outside
		plt.bar(xticks_pos-0.125, lipids_ff_contacts_l2u_outside_by_size_group_nb_sum[:group_gmax], width=0.25, color='#808080', label="before/after")
		#during
		plt.bar(xticks_pos+0.125, lipids_ff_contacts_l2u_during_by_size_group_nb_sum[:group_gmax], width=0.25, color='w', label="during")
	plt.xticks(xticks_pos, xticks_lab, size='x-small')
	plt.yticks(range(int(min(plt.yticks()[0])), int(math.ceil(max(plt.yticks()[0])))+1))
	plt.title("Lower to Upper", size='small')
	fontP.set_size("small")
	if np.size(lipids_ff_l2u_index) > 0:
		ax2.legend(prop=fontP)

	#legend, labels, etc
	#-------------------
	#ax1.set_ylabel('nb of peptides', fontsize="small")
	#ax2.set_ylabel('% of peptides', fontsize="small")
	ax1.tick_params(axis='both', direction='out')
	ax2.tick_params(axis='both', direction='out')
	ax1.spines["right"].set_visible(False)								# remove unneeded axes
	ax1.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["top"].set_visible(False)
	ax1.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax1.get_yaxis().tick_left()
	ax2.get_xaxis().tick_bottom()  										# remove unneeded ticks 
	ax2.get_yaxis().tick_left()
	ax1.set_xlim(0.5, group_gmax + 0.5)
	ax2.set_xlim(0.5, group_gmax + 0.5)
	ax1.set_ylim(ymin=0)
	ax2.set_ylim(ymin=0)
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="x-small")
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="x-small")
	#plt.setp(ax1.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")	
	#plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60, fontsize="x-small")

	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
	ax1.yaxis.get_major_formatter().set_powerlimits((-3, 3))

	#save figure
	#-----------
	plt.subplots_adjust(left=0.1, right=0.96, hspace=0.3, top=0.92, bottom=0.07)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return

#contacts distribution along local normal: all sizes
def write_ff_ctcts_profile_during_all():
	
	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		#create file
		filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/xvg/ff_ctcts_profile_u2l_during_all.xvg'
		output_xvg = open(filename, 'w')	

		#general info as comment
		output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
		output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#xvg metadata
		output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with all TM cluster sizes\"\n")
		output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
		output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
		output_xvg.write("@ s0 legend \"basic\"\n")
		output_xvg.write("@ s1 legend \"polar\"\n")
		output_xvg.write("@ s2 legend \"hydrophobic\"\n")
		output_xvg.write("@ s3 legend \"bb_only\"\n")
		output_xvg.write("@ s4 legend \"total (avg)\"\n")
		output_xvg.write("@ s5 legend \"basic (std)\"\n")
		output_xvg.write("@ s6 legend \"polar (std)\"\n")
		output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
		output_xvg.write("@ s8 legend \"bb_only (std)\"\n")

		#data
		for n in range(0,2*bins_nb):
			results = str(bins_labels[n]) 
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_during_profile_avg[t,n])
			results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_u2l_during_profile_avg[:,n]))
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_during_profile_std[t,n])
			output_xvg.write(results + "\n")	
		output_xvg.close()

	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		#create file
		filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/xvg/ff_ctcts_profile_l2u_during_all.xvg'
		output_xvg = open(filename, 'w')	

		#general info as comment
		output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
		output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#xvg metadata
		output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with all TM cluster sizes\"\n")
		output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
		output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
		output_xvg.write("@ s0 legend \"basic\"\n")
		output_xvg.write("@ s1 legend \"polar\"\n")
		output_xvg.write("@ s2 legend \"hydrophobic\"\n")
		output_xvg.write("@ s3 legend \"bb_only\"\n")
		output_xvg.write("@ s4 legend \"total (avg)\"\n")
		output_xvg.write("@ s5 legend \"basic (std)\"\n")
		output_xvg.write("@ s6 legend \"polar (std)\"\n")
		output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
		output_xvg.write("@ s8 legend \"bb_only (std)\"\n")

		#data
		for n in range(0,2*bins_nb):
			results = str(bins_labels[n]) 
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_during_profile_avg[t,n])
			results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_l2u_during_profile_avg[:,n]))
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_during_profile_std[t,n])
			output_xvg.write(results + "\n")	
		output_xvg.close()

	return
def write_ff_ctcts_profile_outside_all():
	
	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		#create file
		filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/xvg/ff_ctcts_profile_u2l_outside_all.xvg'
		output_xvg = open(filename, 'w')	
	
		#general info as comment
		output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
		output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#xvg metadata
		output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with all TM cluster sizes\"\n")
		output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
		output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
		output_xvg.write("@ s0 legend \"basic\"\n")
		output_xvg.write("@ s1 legend \"polar\"\n")
		output_xvg.write("@ s2 legend \"hydrophobic\"\n")
		output_xvg.write("@ s3 legend \"bb_only\"\n")
		output_xvg.write("@ s4 legend \"total (avg)\"\n")
		output_xvg.write("@ s5 legend \"basic (std)\"\n")
		output_xvg.write("@ s6 legend \"polar (std)\"\n")
		output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
		output_xvg.write("@ s8 legend \"bb_only (std)\"\n")

		#data
		for n in range(0,2*bins_nb):
			results = str(bins_labels[n]) 
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_outside_profile_avg[t,n])
			results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_u2l_outside_profile_avg[:,n]))
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_outside_profile_std[t,n])
			output_xvg.write(results + "\n")	
		output_xvg.close()

	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		#create file
		filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/xvg/ff_ctcts_profile_l2u_outside_all.xvg'
		output_xvg = open(filename, 'w')	
	
		#general info as comment
		output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
		output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		if args.m_algorithm == "density":
			output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
			output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
			output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
		elif args.m_algorithm == "min":
			output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		else:
			output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
			output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
		output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
	
		#xvg metadata
		output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with all TM cluster sizes\"\n")
		output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
		output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
		output_xvg.write("@ s0 legend \"basic\"\n")
		output_xvg.write("@ s1 legend \"polar\"\n")
		output_xvg.write("@ s2 legend \"hydrophobic\"\n")
		output_xvg.write("@ s3 legend \"bb_only\"\n")
		output_xvg.write("@ s4 legend \"total (avg)\"\n")
		output_xvg.write("@ s5 legend \"basic (std)\"\n")
		output_xvg.write("@ s6 legend \"polar (std)\"\n")
		output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
		output_xvg.write("@ s8 legend \"bb_only (std)\"\n")

		#data
		for n in range(0,2*bins_nb):
			results = str(bins_labels[n]) 
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_outside_profile_avg[t,n])
			results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_l2u_outside_profile_avg[:,n]))
			for t in range(0,4):
				results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_outside_profile_std[t,n])
			output_xvg.write(results + "\n")	
		output_xvg.close()

	return
def graph_ff_ctcts_profile_during_all():

	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		#create filenames
		filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/ff_ctcts_profile_u2l_during_all.svg'
		filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/png/ff_ctcts_profile_u2l_during_all.png'

		#create figure
		fig = plt.figure(figsize=(8, 6.2))
		fig.suptitle("average distribution of contacts along the local normal to the bilayer")
		
		#plot data: % of contacts
		ax = fig.add_subplot(111)
		plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_avg[0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_avg[0,:] - lipids_ff_contacts_u2l_during_profile_std[0,:], lipids_ff_contacts_u2l_during_profile_avg[0,:] + lipids_ff_contacts_u2l_during_profile_std[0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_avg[1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_avg[1,:] - lipids_ff_contacts_u2l_during_profile_std[1,:], lipids_ff_contacts_u2l_during_profile_avg[1,:] + lipids_ff_contacts_u2l_during_profile_std[1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_avg[2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_avg[2,:] - lipids_ff_contacts_u2l_during_profile_std[2,:], lipids_ff_contacts_u2l_during_profile_avg[2,:] + lipids_ff_contacts_u2l_during_profile_std[2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_avg[3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_avg[3,:] - lipids_ff_contacts_u2l_during_profile_std[3,:], lipids_ff_contacts_u2l_during_profile_avg[3,:] + lipids_ff_contacts_u2l_during_profile_std[3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)

		#plot data: limits
		tmp_max = np.max(np.sum(lipids_ff_contacts_u2l_during_profile_avg, axis = 0)) + 0.5
		plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
		ax.set_ylim(0, tmp_max)

		#formatting
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		fontP.set_size("small")
		ax.legend(prop=fontP)
		plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
		plt.ylabel('contacts distribution (%)')
		
		#save figure
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		#create filenames
		filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/ff_ctcts_profile_l2u_during_all.svg'
		filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/png/ff_ctcts_profile_l2u_during_all.png'

		#create figure
		fig = plt.figure(figsize=(8, 6.2))
		fig.suptitle("average distribution of contacts along the local normal to the bilayer")
		
		#plot data: % of contacts
		ax = fig.add_subplot(111)
		plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_avg[0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_avg[0,:] - lipids_ff_contacts_l2u_during_profile_std[0,:], lipids_ff_contacts_l2u_during_profile_avg[0,:] + lipids_ff_contacts_l2u_during_profile_std[0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_avg[1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_avg[1,:] - lipids_ff_contacts_l2u_during_profile_std[1,:], lipids_ff_contacts_l2u_during_profile_avg[1,:] + lipids_ff_contacts_l2u_during_profile_std[1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_avg[2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_avg[2,:] - lipids_ff_contacts_l2u_during_profile_std[2,:], lipids_ff_contacts_l2u_during_profile_avg[2,:] + lipids_ff_contacts_l2u_during_profile_std[2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_avg[3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_avg[3,:] - lipids_ff_contacts_l2u_during_profile_std[3,:], lipids_ff_contacts_l2u_during_profile_avg[3,:] + lipids_ff_contacts_l2u_during_profile_std[3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)

		#plot data: limits
		tmp_max = np.max(np.sum(lipids_ff_contacts_l2u_during_profile_avg, axis = 0)) + 0.5
		plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
		ax.set_ylim(0, tmp_max)

		#formatting
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		fontP.set_size("small")
		ax.legend(prop=fontP)
		plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
		plt.ylabel('contacts distribution (%)')
		
		#save figure
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	return
def graph_ff_ctcts_profile_outside_all():

	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		#create filenames
		filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/ff_ctcts_profile_u2l_outside_all.svg'
		filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/png/ff_ctcts_profile_u2l_outside_all.png'

		#create figure
		fig = plt.figure(figsize=(8, 6.2))
		fig.suptitle("average distribution of contacts along the local normal to the bilayer")
		
		#plot data: % of contacts
		ax = fig.add_subplot(111)
		plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[0,:] - lipids_ff_contacts_u2l_outside_profile_std[0,:], lipids_ff_contacts_u2l_outside_profile_avg[0,:] + lipids_ff_contacts_u2l_outside_profile_std[0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[1,:] - lipids_ff_contacts_u2l_outside_profile_std[1,:], lipids_ff_contacts_u2l_outside_profile_avg[1,:] + lipids_ff_contacts_u2l_outside_profile_std[1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[2,:] - lipids_ff_contacts_u2l_outside_profile_std[2,:], lipids_ff_contacts_u2l_outside_profile_avg[2,:] + lipids_ff_contacts_u2l_outside_profile_std[2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
		plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_avg[3,:] - lipids_ff_contacts_u2l_outside_profile_std[3,:], lipids_ff_contacts_u2l_outside_profile_avg[3,:] + lipids_ff_contacts_u2l_outside_profile_std[3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)

		#plot data: limits
		tmp_max = np.max(np.sum(lipids_ff_contacts_u2l_outside_profile_avg, axis = 0)) + 0.5
		plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
		ax.set_ylim(0, tmp_max)

		#formatting
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		fontP.set_size("small")
		ax.legend(prop=fontP)
		plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
		plt.ylabel('contacts distribution (%)')
		
		#save figure
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		#create filenames
		filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/ff_ctcts_profile_l2u_outside_all.svg'
		filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/png/ff_ctcts_profile_l2u_outside_all.png'

		#create figure
		fig = plt.figure(figsize=(8, 6.2))
		fig.suptitle("average distribution of contacts along the local normal to the bilayer")
		
		#plot data: % of contacts
		ax = fig.add_subplot(111)
		plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[0,:] - lipids_ff_contacts_l2u_outside_profile_std[0,:], lipids_ff_contacts_l2u_outside_profile_avg[0,:] + lipids_ff_contacts_l2u_outside_profile_std[0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[1,:] - lipids_ff_contacts_l2u_outside_profile_std[1,:], lipids_ff_contacts_l2u_outside_profile_avg[1,:] + lipids_ff_contacts_l2u_outside_profile_std[1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[2,:] - lipids_ff_contacts_l2u_outside_profile_std[2,:], lipids_ff_contacts_l2u_outside_profile_avg[2,:] + lipids_ff_contacts_l2u_outside_profile_std[2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
		plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
		plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_avg[3,:] - lipids_ff_contacts_l2u_outside_profile_std[3,:], lipids_ff_contacts_l2u_outside_profile_avg[3,:] + lipids_ff_contacts_l2u_outside_profile_std[3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)

		#plot data: limits
		tmp_max = np.max(np.sum(lipids_ff_contacts_l2u_outside_profile_avg, axis = 0)) + 0.5
		plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
		plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
		ax.set_ylim(0, tmp_max)

		#formatting
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		fontP.set_size("small")
		ax.legend(prop=fontP)
		plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
		plt.ylabel('contacts distribution (%)')
		
		#save figure
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	return

#contacts distribution along local normal: by groups
def write_ff_ctcts_profile_during_groups():
	
	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		for g_index in range(0,group_gmax):

			#create file
			filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/xvg/ff_ctcts_profile_u2l_during_' + str(groups_labels[g_index]) + '.xvg'
			output_xvg = open(filename, 'w')	
			output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		
			#general info as comment
			output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
			output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
			if args.m_algorithm == "density":
				output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
				output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
				output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
			elif args.m_algorithm == "min":
				output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			else:
				output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
		
			#xvg metadata
			output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with clusters of size " + str(groups_labels[g_index]) + "\"\n")
			output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
			output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
			output_xvg.write("@ autoscale ONREAD xaxes\n")
			output_xvg.write("@ TYPE XY\n")
			output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
			output_xvg.write("@ legend on\n")
			output_xvg.write("@ legend box on\n")
			output_xvg.write("@ legend loctype view\n")
			output_xvg.write("@ legend 0.98, 0.8\n")
			output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
			output_xvg.write("@ s0 legend \"basic\"\n")
			output_xvg.write("@ s1 legend \"polar\"\n")
			output_xvg.write("@ s2 legend \"hydrophobic\"\n")
			output_xvg.write("@ s3 legend \"bb_only\"\n")
			output_xvg.write("@ s4 legend \"total (avg)\"\n")
			output_xvg.write("@ s5 legend \"basic (std)\"\n")
			output_xvg.write("@ s6 legend \"polar (std)\"\n")
			output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
			output_xvg.write("@ s8 legend \"bb_only (std)\"\n")
	
			#data
			for n in range(0,2*bins_nb):
				results = str(bins_labels[n]) 
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][t,n])
				results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][:,n]))
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_during_profile_groups_std[g_index][t,n])
				output_xvg.write(results + "\n")	
			output_xvg.close()

	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		for g_index in range(0,group_gmax):

			#create file
			filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/xvg/ff_ctcts_profile_l2u_during_' + str(groups_labels[g_index]) + '.xvg'
			output_xvg = open(filename, 'w')	
			output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		
			#general info as comment
			output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
			output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
			if args.m_algorithm == "density":
				output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
				output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
				output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
			elif args.m_algorithm == "min":
				output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			else:
				output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
		
			#xvg metadata
			output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with clusters of size " + str(groups_labels[g_index]) + "\"\n")
			output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
			output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
			output_xvg.write("@ autoscale ONREAD xaxes\n")
			output_xvg.write("@ TYPE XY\n")
			output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
			output_xvg.write("@ legend on\n")
			output_xvg.write("@ legend box on\n")
			output_xvg.write("@ legend loctype view\n")
			output_xvg.write("@ legend 0.98, 0.8\n")
			output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
			output_xvg.write("@ s0 legend \"basic\"\n")
			output_xvg.write("@ s1 legend \"polar\"\n")
			output_xvg.write("@ s2 legend \"hydrophobic\"\n")
			output_xvg.write("@ s3 legend \"bb_only\"\n")
			output_xvg.write("@ s4 legend \"total (avg)\"\n")
			output_xvg.write("@ s5 legend \"basic (std)\"\n")
			output_xvg.write("@ s6 legend \"polar (std)\"\n")
			output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
			output_xvg.write("@ s8 legend \"bb_only (std)\"\n")
	
			#data
			for n in range(0,2*bins_nb):
				results = str(bins_labels[n]) 
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][t,n])
				results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][:,n]))
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_during_profile_groups_std[g_index][t,n])
				output_xvg.write(results + "\n")	
			output_xvg.close()

	return
def write_ff_ctcts_profile_outside_groups():
	
	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		for g_index in range(0,group_gmax):

			#create file
			filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/xvg/ff_ctcts_profile_u2l_outside_' + str(groups_labels[g_index]) + '.xvg'
			output_xvg = open(filename, 'w')	
			output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		
			#general info as comment
			output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
			output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
			if args.m_algorithm == "density":
				output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
				output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
				output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
			elif args.m_algorithm == "min":
				output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			else:
				output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
		
			#xvg metadata
			output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with clusters of size " + str(groups_labels[g_index]) + "\"\n")
			output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
			output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
			output_xvg.write("@ autoscale ONREAD xaxes\n")
			output_xvg.write("@ TYPE XY\n")
			output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
			output_xvg.write("@ legend on\n")
			output_xvg.write("@ legend box on\n")
			output_xvg.write("@ legend loctype view\n")
			output_xvg.write("@ legend 0.98, 0.8\n")
			output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
			output_xvg.write("@ s0 legend \"basic\"\n")
			output_xvg.write("@ s1 legend \"polar\"\n")
			output_xvg.write("@ s2 legend \"hydrophobic\"\n")
			output_xvg.write("@ s3 legend \"bb_only\"\n")
			output_xvg.write("@ s4 legend \"total (avg)\"\n")
			output_xvg.write("@ s5 legend \"basic (std)\"\n")
			output_xvg.write("@ s6 legend \"polar (std)\"\n")
			output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
			output_xvg.write("@ s8 legend \"bb_only (std)\"\n")
	
			#data
			for n in range(0,2*bins_nb):
				results = str(bins_labels[n]) 
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][t,n])
				results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][:,n]))
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][t,n])
				output_xvg.write(results + "\n")	
			output_xvg.close()

	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		for g_index in range(0,group_gmax):

			#create file
			filename = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/xvg/ff_ctcts_profile_l2u_outside_' + str(groups_labels[g_index]) + '.xvg'
			output_xvg = open(filename, 'w')	
			output_xvg.write("#[flipflopping lipids contact statistics - written by ff_contacts v" + str(version_nb) +"]\n")
		
			#general info as comment
			output_xvg.write("#-nb of proteins: " + str(proteins_nb) + "\n")
			output_xvg.write("#-nb frames read: " + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
			if args.m_algorithm == "density":
				output_xvg.write("#-method cluster: density based algorithm using distances between proteins COGs\n")	
				output_xvg.write("# -> radius search: " + str(args.dbscan_dist) + " Angstrom\n")
				output_xvg.write("# -> nb neighbours: " + str(args.dbscan_nb) + "\n")
			elif args.m_algorithm == "min":
				output_xvg.write("#-method cluster: connectivity algorithm using minimum distance between proteins\n")
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			else:
				output_xvg.write("#-method cluster: connectivity algorithm using distance between the center of geometry of proteins\n")	
				output_xvg.write("# -> connect cutoff: " + str(args.nx_cutoff) + " Angstrom\n")
			output_xvg.write("#-cutoff distance for protein-lipid contact: " + str(args.cutoff_pl) + " Angstrom\n")
		
			#xvg metadata
			output_xvg.write("@ title \"Average distribution along the local normal to the bilayer of contacts with clusters of size " + str(groups_labels[g_index]) + "\"\n")
			output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
			output_xvg.write("@ yaxis label \"contacts distribution (% of total contacts)\"\n")
			output_xvg.write("@ autoscale ONREAD xaxes\n")
			output_xvg.write("@ TYPE XY\n")
			output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
			output_xvg.write("@ legend on\n")
			output_xvg.write("@ legend box on\n")
			output_xvg.write("@ legend loctype view\n")
			output_xvg.write("@ legend 0.98, 0.8\n")
			output_xvg.write("@ legend length " + str(2*len(["basic","polar","hydrophobic","bb_only"])+1) + "\n")
			output_xvg.write("@ s0 legend \"basic\"\n")
			output_xvg.write("@ s1 legend \"polar\"\n")
			output_xvg.write("@ s2 legend \"hydrophobic\"\n")
			output_xvg.write("@ s3 legend \"bb_only\"\n")
			output_xvg.write("@ s4 legend \"total (avg)\"\n")
			output_xvg.write("@ s5 legend \"basic (std)\"\n")
			output_xvg.write("@ s6 legend \"polar (std)\"\n")
			output_xvg.write("@ s7 legend \"hydrophobic (std)\"\n")
			output_xvg.write("@ s8 legend \"bb_only (std)\"\n")
	
			#data
			for n in range(0,2*bins_nb):
				results = str(bins_labels[n]) 
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][t,n])
				results += "	" + "{:.6e}".format(np.sum(lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][:,n]))
				for t in range(0,4):
					results += "	" + "{:.6e}".format(lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][t,n])
				output_xvg.write(results + "\n")	
			output_xvg.close()

	return
def graph_ff_ctcts_profile_during_groups():

	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		for g_index in range(0,group_gmax):
			#create filenames
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/ff_ctcts_profile_u2l_during_' + str(groups_labels[g_index]) + '.svg'
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/png/ff_ctcts_profile_u2l_during_' + str(groups_labels[g_index]) + '.svg'
	
			#create figure
			fig = plt.figure(figsize=(8, 6.2))
			fig.suptitle("average distribution of contacts along the local normal to the bilayer")
						
			#plot data: % of contacts
			ax = fig.add_subplot(111)
			plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][0,:] - lipids_ff_contacts_u2l_during_profile_groups_std[g_index][0,:], lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][0,:] + lipids_ff_contacts_u2l_during_profile_groups_std[g_index][0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][1,:] - lipids_ff_contacts_u2l_during_profile_groups_std[g_index][1,:], lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][1,:] + lipids_ff_contacts_u2l_during_profile_groups_std[g_index][1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][2,:] - lipids_ff_contacts_u2l_during_profile_groups_std[g_index][2,:], lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][2,:] + lipids_ff_contacts_u2l_during_profile_groups_std[g_index][2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][3,:] - lipids_ff_contacts_u2l_during_profile_groups_std[g_index][3,:], lipids_ff_contacts_u2l_during_profile_groups_avg[g_index][3,:] + lipids_ff_contacts_u2l_during_profile_groups_std[g_index][3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)
	
			#plot data: limits
			tmp_max = np.max(np.sum(lipids_ff_contacts_u2l_during_profile_groups_avg[g_index], axis = 0)) + 0.5
			plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
			ax.set_ylim(0, tmp_max)
	
			#formatting
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('left')
			fontP.set_size("small")
			ax.legend(prop=fontP)
			plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
			plt.ylabel('contacts distribution (%)')
			
			#save figure
			fig.savefig(filename_png)
			fig.savefig(filename_svg)
			plt.close()
		
	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		for g_index in range(0,group_gmax):
			#create filenames
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/ff_ctcts_profile_l2u_during_' + str(groups_labels[g_index]) + '.svg'
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/png/ff_ctcts_profile_l2u_during_' + str(groups_labels[g_index]) + '.svg'
	
			#create figure
			fig = plt.figure(figsize=(8, 6.2))
			fig.suptitle("average distribution of contacts along the local normal to the bilayer")
			
			#plot data: % of contacts
			ax = fig.add_subplot(111)
			plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][0,:] - lipids_ff_contacts_l2u_during_profile_groups_std[g_index][0,:], lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][0,:] + lipids_ff_contacts_l2u_during_profile_groups_std[g_index][0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][1,:] - lipids_ff_contacts_l2u_during_profile_groups_std[g_index][1,:], lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][1,:] + lipids_ff_contacts_l2u_during_profile_groups_std[g_index][1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][2,:] - lipids_ff_contacts_l2u_during_profile_groups_std[g_index][2,:], lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][2,:] + lipids_ff_contacts_l2u_during_profile_groups_std[g_index][2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][3,:] - lipids_ff_contacts_l2u_during_profile_groups_std[g_index][3,:], lipids_ff_contacts_l2u_during_profile_groups_avg[g_index][3,:] + lipids_ff_contacts_l2u_during_profile_groups_std[g_index][3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)
	
			#plot data: limits
			tmp_max = np.max(np.sum(lipids_ff_contacts_l2u_during_profile_groups_avg[g_index], axis = 0)) + 0.5
			plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
			ax.set_ylim(0, tmp_max)
	
			#formatting
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('left')
			fontP.set_size("small")
			ax.legend(prop=fontP)
			plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
			plt.ylabel('contacts distribution (%)')
			
			#save figure
			fig.savefig(filename_png)
			fig.savefig(filename_svg)
			plt.close()	

	return
def graph_ff_ctcts_profile_outside_groups():

	#upper to lower
	#==============
	if np.size(lipids_ff_u2l_index) > 0:
		for g_index in range(0,group_gmax):
			#create filenames
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/ff_ctcts_profile_u2l_outside_' + str(groups_labels[g_index]) + '.svg'
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/png/ff_ctcts_profile_u2l_outside_' + str(groups_labels[g_index]) + '.svg'
	
			#create figure
			fig = plt.figure(figsize=(8, 6.2))
			fig.suptitle("average distribution of contacts along the local normal to the bilayer")
			
			#plot data: % of contacts
			ax = fig.add_subplot(111)
			plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][0,:] - lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][0,:], lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][0,:] + lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][1,:] - lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][1,:], lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][1,:] + lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][2,:] - lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][2,:], lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][2,:] + lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
			plt.fill_between(bins_range, lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][3,:] - lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][3,:], lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index][3,:] + lipids_ff_contacts_u2l_outside_profile_groups_std[g_index][3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)
	
			#plot data: limits
			tmp_max = np.max(np.sum(lipids_ff_contacts_u2l_outside_profile_groups_avg[g_index], axis = 0)) + 0.5
			plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
			ax.set_ylim(0, tmp_max)
	
			#formatting
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('left')
			fontP.set_size("small")
			ax.legend(prop=fontP)
			plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
			plt.ylabel('contacts distribution (%)')
			
			#save figure
			fig.savefig(filename_png)
			fig.savefig(filename_svg)
			plt.close()
		
	#lower to upper
	#==============
	if np.size(lipids_ff_l2u_index) > 0:
		for g_index in range(0,group_gmax):
			#create filenames
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/ff_ctcts_profile_l2u_outside_' + str(groups_labels[g_index]) + '.svg'
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/TM_profile/groups/png/ff_ctcts_profile_l2u_outside_' + str(groups_labels[g_index]) + '.svg'
	
			#create figure
			fig = plt.figure(figsize=(8, 6.2))
			fig.suptitle("average distribution of contacts along the local normal to the bilayer")
			
			#plot data: % of contacts
			ax = fig.add_subplot(111)
			plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][0,:], color = residues_types_colours["basic"], linewidth = 2.0, label = 'basic')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][0,:] - lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][0,:], lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][0,:] + lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][0,:], color = residues_types_colours["basic"], edgecolor = residues_types_colours["basic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][1,:], color = residues_types_colours["polar"], linewidth = 2.0, label = 'polar')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][1,:] - lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][1,:], lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][1,:] + lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][1,:], color = residues_types_colours["polar"], edgecolor = residues_types_colours["polar"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][2,:], color = residues_types_colours["hydrophobic"], linewidth = 2.0, label = 'hydrophobic')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][2,:] - lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][2,:], lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][2,:] + lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][2,:], color = residues_types_colours["hydrophobic"], edgecolor = residues_types_colours["hydrophobic"], linewidth = 0, alpha = 0.2)
			plt.plot(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][3,:], color = residues_types_colours["bb_only"], linewidth = 2.0, label = 'bb_only')
			plt.fill_between(bins_range, lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][3,:] - lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][3,:], lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index][3,:] + lipids_ff_contacts_l2u_outside_profile_groups_std[g_index][3,:], color = residues_types_colours["bb_only"], edgecolor = residues_types_colours["bb_only"], linewidth = 0, alpha = 0.2)
	
			#plot data: limits
			tmp_max = np.max(np.sum(lipids_ff_contacts_l2u_outside_profile_groups_avg[g_index], axis = 0)) + 0.5
			plt.vlines(z_lower_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(z_upper_avg, 0, tmp_max, linestyles = 'dashed')
			plt.vlines(0, 0, tmp_max, linestyles = 'dashdot')
			ax.set_ylim(0, tmp_max)
	
			#formatting
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('left')
			fontP.set_size("small")
			ax.legend(prop=fontP)
			plt.xlabel('z coordinate centered on middle of bilayer [$\AA$]')
			plt.ylabel('contacts distribution (%)')
			
			#save figure
			fig.savefig(filename_png)
			fig.savefig(filename_svg)
			plt.close()	

	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
# process inputs
#=========================================================================================
#data loading
set_lipids_beads()
load_MDA_universe()
if args.selection_file_ff != "no":
	identify_ff()
identify_proteins()
if args.cluster_groups_file != "no":
	initialise_groups()
identify_leaflets()

#create data structures
print "\nInitialising data structures..."
data_struct_time()
data_ff_contacts()

#=========================================================================================
# process frames
#=========================================================================================
print "\nCalculating sizes sampled by flip-flopping lipids..."

for f_index in range(0,nb_frames_to_process):
	#frame properties
	ts = U.trajectory[frames_to_process[f_index]]
	f_time = ts.time/float(1000)
	f_nb = ts.frame
	frames_nb[f_index] = f_nb
	frames_time[f_index] = f_time
	box_dim = U.trajectory.ts.dimensions
	
	#process ff lipids
	identify_ff_contacts(box_dim, f_time, f_index)
	
print ''

#=========================================================================================
# process data
#=========================================================================================
print "\nCalculating statistics..."
calc_stats_ctcts()
		
#=========================================================================================
# produce outputs
#=========================================================================================
print "\nWriting outputs..."
#distribution over types
write_ff_ctcts_by_type()
graph_ff_ctcts_by_type()
write_ff_ctcts_by_type_sum()
graph_ff_ctcts_by_type_sum()
#distribution over sizes
write_ff_ctcts_by_size()
graph_ff_ctcts_by_size()
write_ff_ctcts_by_size_sum()
graph_ff_ctcts_by_size_sum()

if args.cluster_groups_file != "no":
	write_ff_ctcts_by_group()
	graph_ff_ctcts_by_group()
	write_ff_ctcts_by_group_sum()
	graph_ff_ctcts_by_group_sum()
if args.profile: 
	write_ff_ctcts_profile_during_all()
	write_ff_ctcts_profile_outside_all()
	graph_ff_ctcts_profile_during_all()
	graph_ff_ctcts_profile_outside_all()
	if args.cluster_groups_file != "no":
		write_ff_ctcts_profile_during_groups()
		write_ff_ctcts_profile_outside_groups()
		graph_ff_ctcts_profile_during_groups()
		graph_ff_ctcts_profile_outside_groups()


#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
