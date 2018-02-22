#!/sansom/s137/bioc1535/anaconda2/bin/python
import sys,os
import numpy 
import MDAnalysis
import MDAnalysis.analysis.leaflet
import matplotlib.pyplot as plt
import numpy.ma
import matplotlib.cm

## This script needs a cenetered fitted trajectory 
## use trjconv -fit transxy+rotxy
##

####################################################
##### Variables you can change #####################
####################################################

#### frames from the trajexctory to analyse
start = 1
finish = 5000

#### species to use for bilayer thickness calculations - for atomistic change to P8
n = "No longer needed"

thickstring = "name PO4 or name PO1 or name PO2 or name PO3"

#### these are the species (atoms or CG particles) for which we find densities

### Uncomment the line below for an atomistic system  #####
#Species = ["P8", "C13",  "CA2",  "C50" , "C43", "C3", "N4", "O4", "O1", "C3", "C39" , "O16", "O35" , "C44" , "C20", "C26",  "C42" , "C41" , "C45" , "C46" , "C47", "C48" ,"C49", "C50", "C4*" , "C3*" ]

### Uncomment the line below for a coarse grained system   #####
Species = ["PO4","NC3","GLH","NH3","GL1","GL2", "C1A", "PO1", "PO2","NCO","AM1","AM2","PO3","ROH" ]

#### Flags for ploting residues on 2D plots
plot_basic = True
plot_acidic = True
plot_side = False
plot_CA = True

### Colorscheme see matplotlib.cm ?? for option or create your own thick is for thickness CM is fro densities
#CM = matplotlib.cm.gray
#CM = matplotlib.cm.summer
### define your own colormap
cdict = {'red': (#(0.0, 0.0, 0.5),
                 (0.0, 1.0, 1.0),
                 (0.25, 1.0, 1.0),
                 (0.35, 0.0, 0.0),
                 (0.6, 0.0, 0.0),
                 (1.0, 1.0, 0.0)),
         'green': (#(0.0, 0.0, 0.0),
                 (0.0, 0.75, 0.75),
                 (0.25, 1.0, 1.0),
                 (0.35, 1.0, 1.0),
                 (0.6, 0.5, 0.5),
                 (1.0, 1.0, 0.0)),
         'blue': (#(0.0, 0.0, 0.5),
                 (0.0, 0.79, 0.79),
                 (0.25, 1.0, 1.0),
                 (0.35, 1.0, 1.0),
                 (0.6, 0.0, 0.0),
                 (1.0, 0.0, 0.0))}
CM = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

####################################################
##### Variables you shouldn't need to change #######
####################################################

##### only find leaflets for the first 'leafletts' frames
leafletts = 10

##### parameters for discretising the data
xrange = [-50,50] ### Angstroms around protein
yrange = [-50,50] ### Angstroms around protein 
nbins=(100,100)   #### no. of bins in which to discretise data
center=[0,0,0]

### selection strings to select different types of residues
sel_basic = "(resname ARG or resname LYS or resname HIS)"   ## basic residues
sel_acidic = "(resname ASP or resname GLU)"                 ## acidic residues

####################################################
##### Load the trajectory into the universe ########
####################################################

#### load the trajectory
M = MDAnalysis.Universe(sys.argv[1],sys.argv[2])
print "loaded trajectory\n"
print M.trajectory.dt, "ps between frames" 

##### Define seperate atomgropups for protein, caplhas and side chains
Protein = M.selectAtoms("protein")
CA = M.selectAtoms("name B* or name 0* or name 5B* or name 4B*")
Side = M.selectAtoms("protein and not (name B* or name 0* or name 5B* or name 4B*)")

####################################################
##### Initialise the system and directories ########
####################################################

### create directory for files/images
dir = 'analysis/'
if not os.path.exists(dir):
            os.makedirs(dir)

#### create the dictionaries to store the data
species={} ##dictionary containing species atomgroups both leaflet
species0={} ##dictionary containing species atomgroups outer leaflet
species1={} ##dictionary containing species atomgroups inner leaflet
speciesPop0={} ##dictionary containing species populations outer leaflet
speciesPop1={} ##dictionary containing species populations inner leaflet
speciesZ0={} ##dictionary containing species average height of protein in proximity of lipids
speciesZ1={} ##dictionary containing species  average height of protein in proximity of lipids
for s in Species:
    species[s] = M.selectAtoms("name " + s) ## initializ(s)e the atom group dict arrays
    speciesPop0[s] = numpy.zeros(nbins) ## initializ(s)e the dict arrays
    speciesPop1[s] = numpy.zeros(nbins) ## initializ(s)e the dict arrays
    speciesZ0[s] = () ### list of minimum corrdinates of each species
    speciesZ1[s] = ()

Species.append("head")
speciesPop0["head"] = numpy.zeros(nbins) ## initializ(s)e the dict arrays
speciesPop1["head"] = numpy.zeros(nbins) ## initializ(s)e the dict arrays
speciesZ0["head"] = () ### list of minimum corrdinates of each species
speciesZ1["head"] = ()



##### initialize the arrays for the data and histogram
####leaflet 0
thick0 = numpy.zeros(nbins)   #### bilayer thickness based on minimum PO4 distance outer
thick1 = numpy.zeros(nbins)   #### bilayer thickness based on minimum PO4 distance  inner
z0 = numpy.zeros(nbins)   #### bilayer zcoordinate based on PO4 distance outer
z1 = numpy.zeros(nbins)   #### bilayer zcoordinate based on PO4 distance inner



####################################################
#### functions                              ########
####################################################


def get_dist(atomgroup0, atomgroup1, cutoff = 5):
    """distance array between atomgroups"""
    dist = MDAnalysis.analysis.distances.distance_array(atomgroup0.coordinates(), atomgroup1.coordinates())
    min0 = numpy.sort(dist,1)[:,:cutoff]
    min1 = numpy.transpose(numpy.sort(dist,0)[:cutoff])
    mean0 = numpy.mean(min0,1)
    mean1 = numpy.mean(min1,1)
    dist0 = numpy.hstack((atomgroup0.coordinates() , mean0.reshape(-1,1) , min0[:,0].reshape(-1,1)))
    dist1 = numpy.hstack((atomgroup1.coordinates() , mean1.reshape(-1,1) , min1[:,0].reshape(-1,1)))
    return dist0, dist1


def get_leaflets():
    """populates species dictionaries with atomgroups separated into upper and lower leaflets based on the z coordinates """
    #### Run the leaflet finder using species[n] (PO4 for coarse grained sims).
    L = MDAnalysis.analysis.leaflet.LeafletFinder(M, thickstring)
    ###### Determine which leaflet is upper from the C.O.M 	
    if (L.groups(0).centerOfmass()[2] > L.groups(1).centerOfmass()[2]):
        species0["head"] = L.groups(0) #### upper   
        species1["head"] = L.groups(1) #### lower
    else:
        species0["head"] = L.groups(1) #### upper
        species1["head"] = L.groups(0) #### lower
    #### and create the appropriate atom groups
    for s in Species:
        if s != "head":
            species0[s] = species0["head"].residues.selectAtoms("name " + s)
            species1[s] = species1["head"].residues.selectAtoms("name " + s)

def update_2Dhistogram(coordinates, histogram, center, weights=None):
    """ Updates histograms of species populations """
    temp = coordinates - center
    #print center
    ### special case for z-coordinates....
    if weights == "Z":
        weights = temp[:,2]
        #print weights
    #### bin the results
    histogram = histogram + numpy.histogram2d(temp[:,0], temp[:,1] , bins = nbins , range=[xrange, yrange] , weights = weights)[0]
    return histogram


def update_zrange(atomgroup, zrange, center):
    """Appends the maximum and minimum z coordinates of each species to the appropriate dictionary """
    try:
        closeAtoms = (Protein + atomgroup).selectAtoms("protein and around 7 name " + s)
        zrange = numpy.append(zrange , (numpy.min(closeAtoms.coordinates()[:,2] - center[2]) , numpy.max(closeAtoms.coordinates()[:,2] -center[2])))
    except:
        pass
    return zrange

def normalise_histogram(histogram, frequency):
    """Normalises histograms based on species populations/frequencies"""
    zeros = frequency < 1
    histogram = histogram / (frequency + zeros)
    return histogram

def write_z_pdb(pdbfile, leaflet, populations, data, data1):
    """Writes out pdb of z coordinates data"""
    ### dummy histogram
    dummy = numpy.histogram2d([0,0], [0,0] , bins = nbins, range=[xrange, yrange], weights = [0,0])
    #### Correctly formatted 1st part of each line in pdb
    head = "ATOM     01  SS  SUR     1    "
    fname = dir + pdbfile
    with file(fname, 'w') as outfile:
        outfile.write("TITLE  bilayer z-coordinates as b-factors\n" )
        outfile.write("COMMENT lower leaflet has  negative z-coordinate as b-factors so coloring is consistent\n" )
    for i in numpy.arange(len(dummy[1])):
        for j in numpy.arange(len(dummy[2])):
            if populations[i-1][j-1] > 3:
                oneline =  head + "%8.3f" *3  %(dummy[1][i] , dummy[2][j] , data[i-1][j-1]) + " %5.2f" *2 %(1.00, data1[i-1][j-1])
                with file(fname, 'a') as outfile:
                    outfile.write(oneline +"\n")

####################################################
##### Loop through the trajectory and       ########
##### collect data                          ########
####################################################
for ts in M.trajectory[start:finish]:
    print ts
    #### for the first "leafletts" frames return the atomgroup selections for each leflet for all species
    if ((ts.frame - start) < leafletts):
        get_leaflets()

    ##### Calculation of Blayer thickness for each lipid
    ##### and density of NH3, PO4 and GLH     
    dist0, dist1 = get_dist(species0["head"],species1["head"])
    ##### center data around protein
    center = Protein.centroid()
    #center[0] = M.dimensions[0]/2 
    #center[1] = M.dimensions[1]/2
    #center[2] = M.dimensions[2]/2
    #print center
    #### loop through all the species in the dictionary
    ### update population histograms
    for s in Species:
        if species0[s].numberOfAtoms() > 0:
            speciesPop0[s] = update_2Dhistogram(species0[s].coordinates(), speciesPop0[s], center)
            speciesZ0[s] = update_zrange(species0[s], speciesZ0[s], center)
        if species1[s].numberOfAtoms() > 0:       
            speciesPop1[s] = update_2Dhistogram(species1[s].coordinates(), speciesPop1[s], center)
            speciesZ1[s] = update_zrange(species1[s], speciesZ1[s], center)
    ### update thickness histograms
    thick0 = update_2Dhistogram(dist0[:,0:3], thick0, center, weights = dist0[:,4])
    thick1 = update_2Dhistogram(dist1[:,0:3], thick1, center, weights = dist1[:,4])
    ### update zcoordinate histograms
    z0 = update_2Dhistogram(species0["head"].coordinates(), z0, center, weights = "Z")
    z1 = update_2Dhistogram(species1["head"].coordinates(), z1, center, weights = "Z")

####################################################
##### Normalize the histograms and center   ########
##### the z coordinates about the bilayer   ########
##### center                                ########
####################################################
thick0 = normalise_histogram(thick0, speciesPop0["head"])
z0 = normalise_histogram(z0, speciesPop0["head"])
thick1 = normalise_histogram(thick1, speciesPop1["head"])
z1 = normalise_histogram(z1, speciesPop1["head"])

### center z coordinates histograms so the center of the bilayer is at z=0  
meanz =  0.5 * (z0.mean() + z1.mean())
z1 = z1 - meanz
z0 = z0 - meanz

####################################################
##### Write out the pdb zcoordinate data    ########
##### and protein                           ########
####################################################
write_z_pdb("z0xyz.pdb", "upper", speciesPop0["head"], z0, z0 )
write_z_pdb("z1xyz.pdb", "lower", speciesPop1["head"], z1,-z1 ) #### NEEEDS changing....
##### Use the protein from the final frame and translate it relative to bilayer in z0,z1
Protein.set_positions(Protein.coordinates() - center - [0,0,meanz])
##### Write out as pdb using MDA writer
Protein.write( dir + "protein.pdb")
##### Translate it to the middle of the box
Protein.set_positions(Protein.coordinates() + [0,0,meanz])


####################################################
##### Write out the 2D histograms -         ########
##### populations, thickness, zcoordinates  ########
####################################################

#### density scaling factor Angstroms^2 per bin
#### each point represents and instantaneous density of 1 species per bin area
scale = float((xrange[1]-xrange[0])*(yrange[1]-yrange[0]))/(nbins[0]*nbins[1])  ### scale for 1A^2 area
scale = scale * (finish - start)

for data, name in [thick0, "thick0"], [thick1,"thick1"], [z0,"z0"], [z1,"z1"]:
    numpy.savetxt(dir + name + ".txt", data)

for s in Species:
    if species0[s].numberOfAtoms() > 0:
        speciesPop0[s] = speciesPop0[s]/scale
        numpy.savetxt( dir + s + "upper.txt" , speciesPop0[s])
    if species1[s].numberOfAtoms() > 0:
        speciesPop1[s] = speciesPop1[s]/scale
        numpy.savetxt( dir + s + "lower.txt" , speciesPop1[s])

def nearestProt(zmin,zmax):
    """Returns atomgroups nearest the protein"""
    zmin = str(zmin)
    #print zmin , "zmin"
    zmax = str(zmax)
    #print zmax, "zmax"
    try:
        sideatoms = Side.selectAtoms("prop  z <= " + zmax + " and prop z >= " + zmin  )
    except:
        sideatoms = []
    try: 
        basic = sideatoms.selectAtoms(sel_basic)
    except:
        basic = []
    try:
        acidic = sideatoms.selectAtoms(sel_acidic)
    except:
        acidic =[]
    #print basic, acidic, sideatoms
    #### return centered coordinates
    if sideatoms:
        sideatoms = sideatoms.coordinates() 
        if basic:
            basic = basic.coordinates() 
        if acidic:
            acidic = acidic.coordinates()
    return basic, acidic, sideatoms

def get_average(data):
    """returns values for color map based on the average values in data"""
    masked_values = (data <= 0)
    temp = numpy.ma.array(data, mask = masked_values)
    mean = numpy.ma.mean(temp)
    median = numpy.ma.median(temp)
    print median
    return mean, median

def plot_prot(CA, sideatoms, acidic,basic, \
              plot_basic=plot_basic, plot_acidic=plot_acidic, \
              plot_side = plot_side, plot_CA=plot_CA):
    """plot protein coordinates and basic/acid/side particles"""
    ### set fontsize
    plt.rcParams.update({'font.size': 22})
    #### add colorbar to imshow plot
    cb = plt.colorbar()
    ### add labels
    plt.xlabel(r'$\AA$')
    plt.ylabel(r'$\AA$')
    if plot_CA:
        plt.plot(CA.coordinates()[:,0], -CA.coordinates()[:,1], "k-", alpha = 0.5 )
    if (len(sideatoms) and plot_side):
        plt.plot(sideatoms[:,0], -sideatoms[:,1], "ko", alpha=0.2)
    if (len(basic) and plot_basic):
            plt.plot(basic[:,0], -basic[:,1], "ro", alpha = 0.5 )
    if (len(acidic) and plot_acidic):
            plt.plot(acidic[:,0], -acidic[:,1], "bo" , alpha =0.5 )
    return cb 

def plot2Ddens(s, data0, data1, zrange0, zrange1, title ="density", CM=CM):
    """plots the 2D map of properties (densities,thicknesses etc)  around the protein"""
    plt.figure(figsize=(8,16))
    CM.set_under(color = 'k' , alpha = 0.4)
    for i, leaflet, data, zrange in [0,"upper", data0, zrange0], [1,"lower", data1, zrange1]:
        plt.subplot(2,1,i+1)
        if len(zrange) > 0:
            zmin,zmax = numpy.mean(zrange.reshape(-1,2),axis =0)
        else:
            zmin,zmax = 0,0
        basic, acidic, sideatoms = nearestProt(zmin ,zmax)
        mean, median = get_average(data.T)
        #print "mean", mean, "\nmedian" , median
        if title =="distortion":
            data = data - median.T
        plt.imshow(data.T, extent = extent , cmap = CM )
        cb = plot_prot(CA, sideatoms, acidic, basic)
        plt.title(s + " " + title + ": " + leaflet + " leaflet")
        if title == "density":
            plt.clim(vmin=median/3, vmax=median*4)
            cb.set_label(r'species density [$\AA^{-2}$]')
        elif title =="distortion":
            CM.set_over(color = 'w' , alpha = 0.0)
            plt.clim(vmin=-15, vmax=10)
            cb.set_label(r'bilayer distortion [$\AA$]')
    fname =  dir + s + title + '.png' 
    print fname
    plt.savefig(fname)
    plt.close()


def plot2Ddistortion(s, data0, data1, CM=CM):
    plt.figure(figsize=(8,8))
    plt.title(s + " distortion")
    data = (data0 + data1)/2
    mean, median = get_average(data.T)
    data = data - median.T
    plt.imshow(data.T, extent = extent , cmap = CM )
    basic, acidic, sideatoms = nearestProt(-100 ,100)
    cb = plot_prot(CA, sideatoms, acidic, basic)
    plt.clim(vmin=-15, vmax=10)
    cb.set_label(r'bilayer distortion [$\AA$]')
    CM.set_over(color = 'w' , alpha = 0.0)
    CM.set_under(color = 'k' , alpha = 0.4)
    ## save it!
    fname =  dir + s + 'overalldistortion.png' 
    print fname
    plt.savefig(fname)
    plt.close()

extent = [xrange[0], xrange[1] , yrange[0], yrange[1]]
for s in Species:
    plot2Ddens(s, speciesPop0[s], speciesPop1[s], speciesZ0[s], speciesZ1[s], title ="density" )
plot2Ddens("head", thick0, thick1, speciesZ0["head"], speciesZ1["head"], title ="distortion" )
plot2Ddistortion("head", thick0, thick1)

####################################################
##### Map the "n" density onto the bilayer  ########
##### and write out to pdb                  ########
####################################################
for s in Species:
	write_z_pdb(s + "0dens.pdb", "upper", speciesPop0["head"]*scale, z0, speciesPop0[s]*1000 )
for s in Species:
	write_z_pdb(s + "1dens.pdb", "lower", speciesPop1["head"]*scale, z1, speciesPop1[s]*1000 ) 

