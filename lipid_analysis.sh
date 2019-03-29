#!/bin/sh

#
#Author: Heidi Koldsoe
#Date: Nov2013
#Description: 


#defaults:
lipid="all"
sel1offset=0
forcefield="martini22"
proteinsel="name@BB"
cutoff=6
#retrieve menu options
while getopts ":g:x:n:t:o:f:a:l:r:h" opt; do
        case $opt in
        g)
                grofile=$OPTARG
                ;;
        x)
                xtcfile=$OPTARG
                ;;
        n)
                numframes=$OPTARG
                ;;
        t)
                simtime=$OPTARG
                ;;
        o)
                sel1offset=$OPTARG
                ;;
        f)      
                forcefield=$OPTARG
                ;;

	a)
		proteinsel=$OPTARG
                ;;
	l)
                lipid=$OPTARG
                ;;
	r)
		cutoff=$OPTARG
		;;


    h)
        echo ""
        echo "v1.0.0"
        echo " Number of contacs over time protein+lipid"
        echo "Author: Heidi Koldsø"
        echo "                          "
        echo "                          "
        echo ""
        echo "---------------------------------------------------------------------------------------------"
        echo "This script calculate the number of contacts with cutoff between protein and lipid over time    "
        echo "---------------------------------------------------------------------------------------------"
        echo ""
        echo "[Example]"
        echo " USAGE: sh run_prot-lip_nbcontacts.sh -g <grofile>"
        echo "  -x <xtcfile> -n <number of frames > -t <simulation time in ns > -o <residue number offset> [-a <protein seletion separated by @]"
        echo "  [-o <residue number offsert>] [-f <force field >]  [-r <contact cutoff>] [-l <lipid types>]"

        echo ""
        echo ""
        echo ""
        echo " EXAMPLE:     sh run_prot-lip_fractional_interaction.sh -g system.pdb -x xtcfile -n 15000 -t 3000 "
	echo " -o 612 [-a name@BB] [-l all] [-f martini22] [-r 6]"
	echo "			"
        echo ""
        echo ""
        echo ""
        echo "Option    Default          Description                    "
        echo "-----------------------------------------------------"
        echo "-g                : grofile "
        echo "-x                : xtcfile"
        echo "-n                : number of frames "
        echo "-t                : simtime in ns"
        echo "-f    (opt)       : forcefield (default: martini22)"
        echo "-o    (opt)       : resoffset (default 1)"
        echo "-a    (opt)       : protein selection (default BB - Change for other FF (e.g. AT -a protein)"
        echo "-l    (opt)       : lipid (default all)"
	echo "-r    (opt)       : cutoff distance for a contact (default 6)"
      
echo ""
        exit
          ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
        esac
done

ratio=`echo ${numframes}/${simtime}|bc -l`

echo "ratio=$ratio"

#ratio=$(( ${numframes}.00 / ${simtime}.00 ))

#################################
#
#       Make submit.out file with input
#
#############################

if [ -f "protein_lipid_nbcontacts.out" ]; then
        mv protein_lipid_nbcontacts.out protein_lipid_nbcontacts.out.old
fi

echo "input was: sh run_prot-lip_nbcontacts.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}" > protein_lipid_nbcontacts.out



#########################################################
#                                                       #
#   Create the new directories				#
#                                                       #
#########################################################
                
                dir=$PWD
		if [ ! -d "$dir/PROTEIN_LIPID_TIME_cutoff${cutoff}" ]; then
         	mkdir $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}
		fi
		if [ ! -d "$dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/PLOTS" ]; then
		mkdir $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/PLOTS
		fi
		if [ ! -d "$dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA" ]; then
		mkdir $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA
		fi
		cd $dir
		

#########################################################
#                                                       #
#       Get machine info                                #
#                                                       #
#########################################################
machine=`uname -s`
architecture=`uname -m`

 if ([ "$machine" = "Linux" ] && [ "$architecture" = "x86_64" ])
                        then
                        vmddir=/sbcb/packages/opt/Linux_x86_64/vmd/1.9/bin/vmd
			grodir=/sbcb/packages/opt/Linux_x86_64/gromacs/4.5.4/bin
                        echo "vmddir is: $vmddir"
 fi

 if ([ "$machine" = "Linux" ] && [ "$architecture" = "i686" ])
                        then
                        vmddir=/sbcb/packages/opt/Linux_x86/vmd/1.9/bin/vmd
			grodir=/sbcb/packages/opt/Linux_x86/gromacs/4.5.4/bin
                        echo "vmddir is: $vmddir"
 fi

 if ([ "$machine" = "Darwin" ] && [ "$architecture" = "x86_64" ])
                        then
                        vmddir=/sbcb/packages/opt/Darwin_x86_64/vmd/1.9/vmd/vmd_MACOSXX86
			grodir=/sbcb/packages/opt/Darwin_x86_64/gromacs/4.5.4/bin
                        echo "vmddir is: $vmddir"
 fi

 if ([ "$machine" = "Darwin" ] && [ "$architecture" = "i386" ])
                        then
                        vmddir=/sbcb/packages/opt/Darwin_x86/vmd/1.9/vmd/vmd_MACOSXX86
                        grodir=/sbcb/packages/opt/Darwin_x86/gromacs/4.5.4/bin
			echo "vmddir is: $vmddir"
 fi



                #########################################################
                #                                                       #
                #       Run the tcl script making all frames 		#
		#	to be used for movie thorugh vmd        	#
                #                                                       #
                #########################################################

	           if [ "$forcefield" = "atomistic" -o "$forcefield" = "AA" -o "$forcefield" = "AT" -o "$forcefield" = "Atomistic" -o "$forcefield" = "ATOMISTIC" ]
                       then
                        forcefield="AT"
                fi
                if  [ "$forcefield" = "AT" ]
                        then
                        $vmddir -dispdev text -e /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/prot_lip_frames_AT.tcl -args $proteinsel $lipid $grofile $xtcfile $forcefield $sel1offset $cutoff
                fi

                if  [ "$forcefield" != "AT" ]
                        then

                        $vmddir -dispdev text -e /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/prot_lip_frames.tcl -args $proteinsel $lipid $grofile $xtcfile $forcefield $sel1offset $cutoff
                fi


				

#########################################################
#                                                       #
#  Do gnuplot figures of all data files                 #
#                                                       #
#                                                       #
#                                                       #
#########################################################

                #########################################################
                #                                                       #
                #       Test if gnuplot and fonts are installed         #
                #                                                       #
                #########################################################
                if which gnuplot >/dev/null; then
                echo "gnuplot exists"
                else
                echo "gnuplot does not exist"
                echo "please install (e.g. 'sudo apt-get install gnuplot')"
                echo ""
                exit 0
                fi


                if [ ! -f /usr/share/fonts/truetype/msttcorefonts/times.ttf ]
                then
                    echo "Times fonts not installed"
                    echo " Please install (e.g use 'sudo apt-get ttf-mscorefonts-installer')"
                    echo " Set the GDFONTPATH:"
                echo ".BASHRC export GDFONTPATH=/usr/share/fonts/truetype/msttcorefonts"
                echo ".TCSH setenv GDFONTPATH /usr/share/ fonts/truetype/msttcorefonts/"
                fi


                #########################################################
                #                                                       #
                #       Run the gnuplot script with the                 #
                #       number of frames and cutoff as input            #
                #                                                       #
                #########################################################

       lipidlist=`cat $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/liplistout_frames.dat`		


               for lipidtype in $lipidlist
                       do
PROTlines=`awk 'BEGIN {IFS=","} END {print NF}' ${dir}/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/${lipidtype}_num.dat`
                         PROTNUM=$(( $PROTlines + 1 ))

		#### make pdb-files with beta factor as average number:
			                      resnum=1
                        cd $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA
                                        while [ "$resnum" -lt "$PROTNUM" ]; do
        #                                       echo "$resnum = resnum"
                                                avgnum=`cat ${lipidtype}_num.dat | awk -v R=$resnum '{sum+=$R} END {print sum/NR}'`
                                                avgnum_hg=`cat ${lipidtype}_num_hg.dat | awk -v R=$resnum '{sum+=$R} END {print sum/NR}'`
                                                avgrelnum=`cat ${lipidtype}_relnum.dat | awk -v R=$resnum '{sum+=$R} END {print sum/NR}'`
                                                avgrelnum_hg=`cat ${lipidtype}_relnum_hg.dat | awk -v R=$resnum '{sum+=$R} END {print sum/NR}'`
                                                        if [ "$resnum" -eq "1" ]; then
                                                                echo "$avgnum" > ${lipidtype}_mapping_num_beta.dat
                                                                echo "$avgnum_hg" > ${lipidtype}_mapping_num_hg_beta.dat
                                                                echo "$avgrelnum" > ${lipidtype}_mapping_relnum_beta.dat
                                                                echo "$avgrelnum_hg" > ${lipidtype}_mapping_relnum_hg_beta.dat

                                                        fi
                                                        if [ "$resnum" -gt "1" ]; then
                                                                echo "$avgnum" >> ${lipidtype}_mapping_num_beta.dat
                                                                echo "$avgnum_hg" >> ${lipidtype}_mapping_num_hg_beta.dat
                                                                echo "$avgrelnum" >> ${lipidtype}_mapping_relnum_beta.dat
                                                                echo "$avgrelnum_hg" >> ${lipidtype}_mapping_relnum_hg_beta.dat
                                                        fi
                                                resnum=$(($resnum + 1))
        #                                               echo "$resnum = resnum"

                                        done
                                if [ ! -d "$dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/PDBs" ]; then
                                        mkdir $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/PDBs
                                fi
                                cd  $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/PDBs

                                        $vmddir -dispdev text -e /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/Bfactor.tcl -args ../../$grofile $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/${lipidtype}_mapping_num_beta.dat $forcefield Prot_${lipidtype}_mapping_num_beta.pdb $proteinsel
                                        $vmddir -dispdev text -e /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/Bfactor.tcl -args ../../$grofile $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/${lipidtype}_mapping_num_hg_beta.dat $forcefield Prot_${lipidtype}_mapping_num_hg_beta.pdb $proteinsel
                                        $vmddir -dispdev text -e /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/Bfactor.tcl -args ../../$grofile $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/${lipidtype}_mapping_relnum_beta.dat $forcefield Prot_${lipidtype}_mapping_relnum_beta.pdb $proteinsel
                                        $vmddir -dispdev text -e /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/Bfactor.tcl -args ../../$grofile $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/${lipidtype}_mapping_relnum_hg_beta.dat $forcefield Prot_${lipidtype}_mapping_relnum_hg_beta.pdb $proteinsel


		cd  $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA
		#### normalise the average number:
			 if [ -f "${lipidtype}_mapping_relnum_hg_beta.dat" ]
                                then
                                        
                                        summer=`awk '{sum=sum+$1} END {print sum}' ${lipidtype}_mapping_relnum_hg_beta.dat`
                                        echo "summer is $summer"
                                        awk -v summer=$summer '{print $1/summer}' ${lipidtype}_mapping_relnum_hg_beta.dat > ${lipidtype}_mapping_relnum_hg_beta_norm.dat
				fi
			if [ -f "${lipidtype}_mapping_num_hg_beta.dat" ]
                                then

                                        summer=`awk '{sum=sum+$1} END {print sum}' ${lipidtype}_mapping_num_hg_beta.dat`
                                        echo "summer is $summer"
                                        awk -v summer=$summer '{print $1/summer}' ${lipidtype}_mapping_num_hg_beta.dat > ${lipidtype}_mapping_num_hg_beta_norm.dat
                                fi

				if [ -f "${lipidtype}_mapping_relnum_beta.dat" ]
                                then

                                        summer=`awk '{sum=sum+$1} END {print sum}' ${lipidtype}_mapping_relnum_beta.dat`
                                        echo "summer is $summer"
                                        awk -v summer=$summer '{print $1/summer}' ${lipidtype}_mapping_relnum_beta.dat > ${lipidtype}_mapping_relnum_beta_norm.dat
                                fi

				if [ -f "${lipidtype}_mapping_num_beta.dat" ]
                                then

                                        summer=`awk '{sum=sum+$1} END {print sum}' ${lipidtype}_mapping_num_beta.dat`
                                        echo "summer is $summer"
                                        awk -v summer=$summer '{print $1/summer}' ${lipidtype}_mapping_num_beta.dat > ${lipidtype}_mapping_num_beta_norm.dat
                                fi


		  cd ${dir}/PROTEIN_LIPID_TIME_cutoff${cutoff}/PLOTS/
                         cp ${dir}/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/reslistfile_frames.dat ${dir}/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/tmplist_frames.dat
                       #  PROTlines=`awk 'BEGIN {IFS=","} END {print NF}' ${dir}/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/${lipidtype}_num.dat`
                        # PROTNUM=$(( $PROTlines + 1 ))
                        sed -e "s|*|\"|g;s| \"|, \"|g" $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/tmplist_frames.dat > $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/reslabellist_frames.dat
                        PROTLIST=`cat $dir/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/reslabellist_frames.dat`

                        NUMFRAMES=`awk '{print $1}' ${dir}/PROTEIN_LIPID_TIME_cutoff${cutoff}/DATA/${lipidtype}_num.dat | tail -1`
                        cp /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/nb_contact_time_simple.gnu tmp.gnu

                        fulldir=$dir/PROTEIN_LIPID_TIME_cutoff${cutoff}

                         sed -e "s|PROTNUM|$PROTNUM|g;s|SEL1LIST|$PROTLIST|g;s|LIPPI|$lipidtype|g;s|DIR|$fulldir|g;s|starttime|0|g;s|endtime|$simtime|g;s|ratio|$ratio|g" tmp.gnu > nb_contact_time_${lipidtype}.gnu
                                rm tmp.gnu

                         cp /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/nb_lip_time.gnu tmp.gnu
                         sed -e "s|LIPPI|$lipidtype|g;s|DIR|$fulldir|g;s|starttime|0|g;s|endtime|$simtime|g;s|ratio|$ratio|g;s|CUTOFF|$cutoff|g" tmp.gnu > nb_${lipidtype}_time.gnu

				cp /sansom/s91/bioc1127/gp130/SCRIPTS/JM_SCRIPTS/nb_contact_time_histogram.gnu tmp.gnu
				sed -e "s|SEL1LIST|$PROTLIST|g;s|PROTNUM|$PROTNUM|g;s|LIPPI|$lipidtype|g;s|DIR|$fulldir|g;s|starttime|0|g;s|endtime|$simtime|g;s|ratio|$ratio|g;s|CUTOFF|$cutoff|g" tmp.gnu > nb_${lipidtype}_histogram.gnu
                                gnuplot nb_contact_time_${lipidtype}.gnu
                                gnuplot nb_${lipidtype}_time.gnu
				gnuplot nb_${lipidtype}_histogram.gnu
				


	done
				


 echo "		The scripts is done!"
echo " "
echo "		The files can be found in PROTEIN_LIPID_TIME_cutoff${cutoff}/PLOTS"
		    echo ""
		
