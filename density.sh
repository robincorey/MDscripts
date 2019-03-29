#!/bin/bash

#-------------------------------------------------------------
# Notes
#-------------------------------------------------------------

#how to use: ./thickness -h

#TO EDIT: location of python script
path_to_py_file_NIS=/sansom/s98/orie2528/phd/SBCB/scripts/density/ 	#user NIS path to python file
path_to_py_file_local=/Users/Hal/phd/SBCB/scripts/density/			#user customised path (e.g. to use on your laptop)

#-------------------------------------------------------------
# retrieve user defined parameters
#-------------------------------------------------------------

#set default values
xtcfilename="no"
tprfilename="no"
time_start=0
time_end=0
name_specie="none"
name_particle="none"
nb_bins=40
plot_protein="no"
name_water="W"
color_map=1
options=""									#to append to filename to create output folder name
path_to_py_file=$path_to_py_file_NIS		#assumes NIS run by default

#retrieve menu options
while getopts ":a:b:c:e:g:p:r:s:x:w:hl" opt; do
	case $opt in
    a)
	name_particle=$OPTARG
	  ;;
    b)
	time_start=$OPTARG
	  ;;
    c)
	color_map=$OPTARG
	options=${options}_c$color_map
	  ;;
    e)
	time_end=$OPTARG
      ;;
    g)
	nb_bins=$OPTARG
      ;;
    l)
	path_to_py_file=$path_to_py_file_local
	  ;;
    r)
	name_specie=$OPTARG
	  ;;
    p)
	plot_protein="yes"
      ;;
    s)
	tprfilename=$OPTARG
      ;;
    x)
	xtcfilename=$OPTARG
      ;;
    w)
	name_water=$OPTARG
      ;;
    h)
	echo ""	
	echo "*******************"
	echo "v1.0.1"
	echo "Author: Jean Helie"
	echo "*******************"
	echo ""
	echo "This script plots a specie (either a residue or a particle) density in each leaflet as the"
	echo "presence probability per A^2/ns."
	echo ""	
	echo "[Example]"
	echo "PPCS density, using frames every 1 ns between 50 and 500 ns:"
	echo "> density -s tprfile -x xtcfile -b 50 -e 500 -r PPCS"
	echo ""
	echo "TO DO: choose leaflet"
	echo ""
	echo "Option    Default  	Description                    "
	echo "-----------------------------------------------------"
	echo "-x		: trajectory file [.xtc]"
	echo "-s		: tpr file [.tpr]"
	echo "-b		: starting time (ns) for analysis"
	echo "-e		: ending time (ns) for analysis"	
	echo "-r		: name of residue for which to calculate density"	
	echo "-a		: name of atom/particle for which to calculate density (ignored if -r specified)"	
	echo "-g [40]		: nb of bins (per dimension) in the grid"	
	echo "-l [false]	: local run (i.e. not on NIS folder)"
	echo "-p [no]		: plot protein [TO DO]"
	echo "-w [W]		: name of water residue"
	echo "-c [1]		: color map, gnuplot colours by default"
	echo "    2 		: cyan"
	echo "    3 		: olive"
	echo "    4 		: red"
	echo "    5 		: green"
	echo "    6 		: purple"
	echo "    7 		: blue"
	echo "    8 		: orange"
	echo "    9 		: black & white"
	echo ""	
	exit
	  ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
	esac
done

#-------------------------------------------------------------
# process parameters
#-------------------------------------------------------------

echo "Processing input parameters..."

#sanity check
if [ ! -f "$xtcfilename" ]
then
	echo "Error: xtc file not found."
	exit
elif [ ! -f "$tprfilename" ]
then
	echo "Error: tpr file not found."
	exit
elif [ "$time_start" -eq "0" -o "$time_end" -eq "0" ]
then
	echo "Error: both -b and -e must be specified"
	exit
elif [ "$name_specie" = "none" -a "$name_particle" = "none" ]
then
	echo "Error: either -r or -a must be specified."
	exit
fi

#set ouput name
if [ "$name_specie" = "none" ]
then
	output_name=${xtcfilename%.*}_density_${time_start}-${time_end}ns_${name_particle}$options
else
	output_name=${xtcfilename%.*}_density_${time_start}-${time_end}ns_${name_specie}$options
fi

#extract first frame (NB: gromacs -dump option seems utterly useless if the first frame isn't at t=0 - which seems really fucking retarded?).
filename=${xtcfilename%.*}_tmp.gro
echo 0 | gmx_avx trjconv -f $xtcfilename -s $tprfilename -dump $[$time_start*1000] -o first_$filename &>/dev/null

#case particle: exit if selected particle is not present
if [ "$name_specie" = "none" -a $(grep $name_particle first_$filename | wc -l) -eq "0" ]
then
	echo "Error: no such particle found, check option -a."
	rm first_$filename
	exit
fi

#case residue: exit if selected residue not present
if [ "$name_specie" != "none" -a $(grep $name_specie first_$filename | wc -l) -eq "0" ]
then
	echo "Error: no such residue found, check option -r."
	rm first_$filename
	exit
fi 

#otherwise create an index file with the system minus water molecules	
ndxfilename=${filename%.*}_no_$name_water.ndx
echo q | gmx_avx make_ndx -f first_$filename -o tmp_$ndxfilename
tmp_ndx=$(grep '\[' tmp_$ndxfilename | wc -l)
echo -e "! r $name_water \\n del 0-$[$tmp_ndx-1] \\n q \\n" | gmx_avx make_ndx -f first_$filename -n tmp_$ndxfilename -o $ndxfilename
rm tmp_$ndxfilename

#extract gro file of the system minus waters
gmx_avx trjconv -f first_$filename -s $tprfilename -n $ndxfilename -o first_noW_$filename &>/dev/null
rm first_$filename

#extract xtc of the system minus waters from -b to -e
gmx_avx trjconv -f $xtcfilename -s $tprfilename -n $ndxfilename -b $[$time_start*1000] -e $[$time_end*1000] -o noW_$xtcfilename &>/dev/null
xtcfilename=noW_$xtcfilename
	
#-------------------------------------------------------------
# create log file with user's input
#-------------------------------------------------------------

#echo "[command]" >> ${output_name}.log
#echo "python ${path_to_py_file}density.py $xtcfilename $tprfilename $time_start $time_end $nb_snapshots $name_specie $leaflet_out $output_name $ndxfilename $name_particle $plot_protein" >> ${output_name}.log
#echo "" >> ${output_name}.log
#echo "[Options]" >> ${output_name}.log
#echo "xtc 				: $xtcfilename" >> ${output_name}.log
#echo "start 			: ${time_start} ns" >> ${output_name}.log
#echo "end 				: ${time_end} ns" >> ${output_name}.log
#echo "number of bins	: $nb_bins" >> ${output_name}.log
#echo "specie			: $name_specie" >> ${output_name}.log

#-------------------------------------------------------------
# call script
#-------------------------------------------------------------

python ${path_to_py_file}density.py $xtcfilename $tprfilename $time_start $time_end $name_specie $output_name $ndxfilename $name_particle $plot_protein $nb_bins first_noW_$filename $color_map

return_status=$?

#----------------------------------------------------------------
# after successful run
#----------------------------------------------------------------
if [ "$return_status" = "0" ]
then
	echo ""
	echo "Finished successfully! Check output in ./$output_name/" 
else
	echo ""
	echo "Error, early termination! (must implement clean-up when failed...)"
	#rm files...
fi
