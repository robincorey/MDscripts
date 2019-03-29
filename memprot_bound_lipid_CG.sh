#!/usr/bin/perl 

# Location for Gromacs
my $gromacs = "/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin";
my $gmx = "gmx_sse";
my $gmxd = "gmx_sse";

####Location of the shared files
my $shared = "/sansom/s46/stansfeld/MPIMD/";

my $force = "yes" ; # yes or no
my $omp = "no" ; # yes if outer membrane protein or no if alpha-helical

my $netupper = 1.0; # also try 1.0 1.1 1.2?
my $netlower = 0.5;

my $try = 50;

our %lipidhead = ( 0 => 'NC3', 1 => 'NH3', 2 => 'GLH', 3 => 'B1', 4 => 'CNO', 5 => 'PO1', 6 => 'ROH', 7 => 'GL0', 8 => 'B6', 9 => 'S1', 10 => 'B0', 11 => 'NCO', 12 => 'DOH', 13 => 'PO3', 14 => 'PO0', 15 => 'BAS', 16 => 'GM1'); #The idenfity of the atom to use to count the number of lipid particles
our %lipidkey = ( 'DSPC' => 0, 'dspc' => 0, 'DPPC' => 0, 'dppc' => 0, 'DOPC' => 0, 'dopc' => 0,'DHPC' => 0, 'dhpc' => 0, 'DLPC' => 0, 'dlpc' => 0, 'popc' => 0, POPC => 0, 'pope' => 1, 'POPE' => 1, 'pvpe' => 1, 'PVPE' => 1,  'dppe' => 1, 'DPPE' => 1,'dspe' => 1, 'DSPE' => 1, 'dlpe' => 1, 'DLPE' => 1,'POPG' => 2, 'popg' => 2, 'PVPG' => 2, 'pvpg' => 2, 'DSPG' => 2, 'dspg' => 2, 'DPPG' => 2, 'dppg' => 2, 'DLPG' => 2, 'dlpg' => 2, 'DDM' => 3, 'ddm' => 3, 'BOG' => 3, 'bog' => 3, 'POPS' => 4, 'pops' => 4, 'PIP2' => 5, 'pip2' => 5, 'CHOL' => 6, 'chol' => 6, 'CARD' => 7, 'card' => 7, 'DGD' => 8, 'dgd' => 8, 'sqd' => 9, 'SQD' => 9, 'LMG' => 10, 'lmg' => 10, 'POPS' => 11, 'pops' => 11, 'MAG' => 12, 'mag' => 12, 'LPA' => 13, 'lpa' => 13, 'upp' => 14, 'UPP' => 14,'cysd' => 15, 'CYSD' => 15,'cyst' => 15, 'CYST' => 15,'pam' => 15, 'PAM' => 15, 'lipa' => 16, 'LIPA' => 16, 'ramp' => 16, 'RAMP' => 16,'remp' => 16, 'REMP' => 16);

unless ( $ARGV[2] )
{
	die "Usage: <protein-ready.pdb> <number of lipid types to be added> <name of lipid 1> <% of lipid 1> <name of lipid 2> <% of lipids 2>\ne.g. setup.pl 3JYC 1 dppc\ne.g. setup.pl 3JYC 3 card 10 popg 20 pope 70\nN.B. lipid names MUST be be the same as in the directory storing the shared files...\n";
}

my $pdbid = $ARGV[0];
my $lipidtypes = $ARGV[1];

print $pdbid,"\n";

system("rm Protein*.itp");
mkdir "temp-pdb";

if ($lipidtypes == 1){
$ARGV[3] = 100;
}

my @lipidnames;
my @lipidnumbers;

for (my $i = 0; $i < $lipidtypes; $i++)
{
	$lipidnames[$i] = $ARGV[2+(2*$i)];
	$lipidnumbers[$i] = $ARGV[3+(2*$i)];
	$percentage+=$lipidnumbers[$i];
}

unless ($percentage == 100){
	die "\nPercentage of lipids does not equal 100\n"; 
}

for (my $i = 0; $i < $lipidtypes; $i++)
{
	if (($lipidnames[$i] eq CARD)||($lipidnames[$i] eq card)){
	$lipidnumbers[$i] = $lipidnumbers[$i]/2;
	}
}

########  Bondini End ##############

#system("python2.7 $shared/martinize.py -f temp-pdb/chains.pdb -merge all -dssp /sbcb/packages/opt/Linux_x86_64/dssp/bin/dsspcmbi -ff martini22 -v -x protein-cg.pdb -o protein-cg.top -elastic -ef 1000 -el $netlower -eu $netupper -ea 0 -ep 0");

#Copy the required files
my $copyfile = $shared."/mdp_files/*em.mdp";
my $command = "cp ".$copyfile." .";
system $command;

$copyfile = $shared."/wat.pdb";
$command = "cp ".$copyfile." ./wat.pdb";
system $command;

$copyfile = $shared."/martini_v2.2.itp";
$command = "cp ".$copyfile." ./martini_v2.2.itp";
system $command;

$copyfile = $shared."/mdp_files/cgmd-gpu-5.0.mdp";
$command = "cp ".$copyfile." ./md.mdp";
system $command;

#Make the output directory if it does not exist
my $outpath = "MD/";
if (! -d $outpath )
{
        my $command = "mkdir ".$outpath;
        system $command;
}


my $seed = 9999;

open (INB, $pdbid) || die "Not possible to open protein-cg.pdb\n";
while (<INB>){
	chomp;
	($spa, $spb, $spc, $spd, $spe, $spf, $spg, $sph, $spi,$spj) = split(/\S+/, $_);
	($woa, $wob, $woc, $wod, $woe, $wof, $wog, $woh, $woi, $woj) = split;
	if ($woa eq CRYST1){
	$xbox = $wob/10;
	$ybox = $woc/10;
	$zbox = $wod/10;
	}
}

if ($xbox < 12){
	$xbox = 12;
}
if ($ybox < 12){
	$ybox = 12;
}
if ($zbox < 12){
	$zbox = 12;
}

system("mv Protein*.itp protein-cg.itp");
system("mv $pdbid prot.pdb");

$command = "$gromacs/$gmx editconf -f prot.pdb -o cg-temp.pdb -box ".$xbox." ".$ybox." 50 -c";
system $command;
$command = "$gromacs/$gmx grompp -f cgem.mdp -maxwarn 5 -c cg-temp.pdb -o temp_em -p protein.top";
system $command;

$command = "$gromacs/$gmxd mdrun -deffnm temp_em -c temp_em.pdb";
system $command;
$command = "$gromacs/$gmx editconf -f temp_em.pdb -o lipids.pdb -box ".$xbox." ".$ybox." 8 -c";
system $command;

if(($lipidnames[0] eq mag)||($lipidnames[0] eq lpa)||($lipidnames[0] eq upp)){
	$lipidconc = (($xbox*$ybox)*3.5)/100;
}
elsif (($lipidnames[0] ne mag)&&($lipidnames[0] ne lpa)&&($lipidnames[0] ne upp)){
	$lipidconc = (($xbox*$ybox)*2.5)/100;
}

for (my $i = 0; $i < $lipidtypes; $i++ )
{
my $lipidnum = sprintf("%.0f", (($lipidnumbers[$i])*($lipidconc)));
my $command = "$gromacs/$gmx insert-molecules -f lipids.pdb -rot xyz -nmol ".$lipidnum." -try $try -ci ".$shared."/lipids/".$lipidnames[$i].".pdb -o lipids.pdb -radius 0.25";
system $command;
}
#check that nothing has gone wrong with the preparation process - the first residue must be the same in temp.pdb and the original pdb file.

my $originalfile = "protein-cg.pdb";
my $firstresidue;
open (ORIG, $originalfile ) || die "Cannot open $originalfile\n";
while ( my $line = <ORIG> )
{
	next unless ($line =~ m/^ATOM/ );
	$firstresidue = substr($line, 17, 3 );
	last;
}
#print "$firstresidue\n"; 

my $newfile = "lipids.pdb";
my $newfirstresidue;
open (NEW, $newfile ) || die "Cannot open $newfile\n";
while ( my $line = <NEW> )
{
        next unless ($line =~ m/^ATOM/ );
        $newfirstresidue = substr($line, 17, 3 );
        last;
}

unless ( $newfirstresidue eq $firstresidue )
{
	print "Something has gone wrong with lipid addition\n";
	exit;
}

$command = "$gromacs/$gmx editconf -f lipids.pdb -o prot+lipid.pdb -box ".$xbox." ".$ybox." ".$zbox." -c";
system $command;

$command = "$gromacs/$gmx solvate -cp prot+lipid.pdb -o prot+lipid+wat.pdb -cs ".$shared."/wat.pdb -radius 0.24";
system $command;

#### Create the topology file
#### Read the prot+lipid+wat file to work out how many of everything we have
my $file = "prot+lipid+wat.pdb";
my $water = 0;
my @lipids;
for ( my $i = 0; $i < $lipidtypes; $i++ )
{
	$lipids[$i] = 0;
}

open (FILE, $file) || die "Cannot open $file\n";
while (my $line = <FILE> )
{
	next unless ($line =~ m/^ATOM/ );
	my $residue = substr($line, 13, 3);

	for ( my $i = 0; $i < $lipidtypes; $i++ )
	{
		my $key = $lipidkey{$lipidnames[$i]};
		my $match = $lipidhead{$key};
		if ( $residue =~m/$match/ )
		{
			$lipids[$i]++;	
		}
	}
	if ( $residue =~ m/W/ )
	{
		$water++;
	}
}


my $temp = "temp.top";
open (OUT, ">$temp" ) || die "Cannot open $temp for writing\n";
printf (OUT "#include \"martini_v2.2.itp\"\n#include \"protein-cg.itp\"\n\n[ system ]\nSELF-ASSEMBLY\n\n[ molecules ]\nProtein 1\n");
for (my $i = 0; $i < $lipidtypes; $i++ )
{
        my $name = uc ($lipidnames[$i]);
        printf (OUT "%s %d\n", $name, $lipids[$i] );
}
printf (OUT "W %d\n", $water );

$command = "$gromacs/$gmx grompp -maxwarn 5  -f em.mdp -c prot+lipid+wat.pdb -p temp.top -o em.tpr";
system $command;

system("$gromacs/$gmx make_ndx -f prot+lipid+wat.pdb -o SOL_ION.ndx << EOD
del 0-30
aW
q
EOD");

system("$gromacs/$gmx genion -s em.tpr -o prot+lipid+wat+ion.pdb -neutral -conc 0.15 -n SOL_ION.ndx");

my $na = 0;
my $sol = 0;
my $cl = 0;

open (INP, "prot+lipid+wat+ion.pdb") || die "Something went wrong with the ion addition - cannot open prot+lipid+wat+ion.pdb\n";
$atom = 0;
while (<INP>){
	chomp;
	($woa, $wob, $woc, $wod, $woe, $wof, $wog) = split;
	if (($woc eq NA)){
	$na++;
	}
	if (($woc eq CL)){
	$cl++;
	}
	if (($woc eq W)){
	$w++;
	}
}

close (INP);
my $topol = "topol.top";
open (OUT, ">$topol" ) || die "Cannot open $temp for writing\n";
printf (OUT "#include \"martini_v2.2.itp\"\n#define RUBBER_BANDS\n\n#include \"protein-cg.itp\"\n\n[ system ]\nSELF-ASSEMBLY\n\n[ molecules ]\nProtein 1\n");
for (my $i = 0; $i < $lipidtypes; $i++ )
{
	my $name = uc ($lipidnames[$i]);
	printf (OUT "%s %d\n", $name, $lipids[$i] );
}

print OUT "\nW $w";
print OUT "\nNA+ $na";
print OUT "\nCL- $cl";
close(OUT);

#Run the energy minimization
$command = "$gromacs/$gmx grompp -maxwarn 5  -f em.mdp -c prot+lipid+wat+ion.pdb -p topol.top -o em.tpr";
system $command;

$command = "$gromacs/$gmxd mdrun -v -nice 0 -s em.tpr -deffnm out-em -c CG-system.pdb";
system $command;

system("$gromacs/$gmx make_ndx -f CG-system.pdb -o sys.ndx << EOD
del 0
del 1-40
0|rPOP*|rPIP*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDSP*|rTOC*|rCAR*|rDLP*|rSQD*|rDGD*|rLMG*|rrMAG*|rLPA*|rUPP*|rCYST*|rCYSD*|rPAM*|rLIP*|rREM*|rRAM*
1&!0
!1
del 1
name 1 Lipid
name 2 SOL_ION
q
EOD");

$command = "$gromacs/$gmx grompp -maxwarn 5  -f md.mdp -c CG-system.pdb -p topol.top -n sys.ndx -o MD/md.tpr";
system $command;

open (SUB, ">MD/submit.sh" ) || die "Cannot open submit.sh for writing\n";

print SUB "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -N md\n#\$ -r n -j y -cwd  -q all.q\n#\$ -pe openmpi_singlenode 4\n#\$ -l mem_free=1G\nsource /sbcb/packages/modules/init\nmodule add gromacs openmpi\nmpirun gmx mdrun_mpi -v -stepout 1000 -deffnm md\n";

close(SUB);

system ("rm temp* *.list *\#* *.txt ss.xpm scount.xvg posre* mdout.mdp conf.gro");
`rm -rf  temp* *#* structure.txt *.list *.txt out-em.* charge.* *.log prot+* em.tpr protein.top`;
`rm -rf MD/*#* tf*.dat sequ.txts split_results *.dat tm.pdb list_protein_names Yoursequenc princ*.pdb topol_Protein* temp-pdb`;
open (MD, "MD/md.tpr" ) || die "Something went wrong\n";
close (MD);

print "\n Ready for CG Self-Assembly\n\n";
