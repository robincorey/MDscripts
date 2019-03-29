#!/bin/sh

#########################################
#
# Author:       Heidi Koldsoe
# Date:         16.07.12
# Description:  Exchange PC lipids with PE and/or PS lipids
#
#
#########################################
#REMARK 070113: van der waals pf water during solvation has been changes from -vdwd 0.25 to -vdwd 0.21
#REMARK 170413: DPPG addded
### RTEMARk: GM3 is now negatively charged



#########################################################
#                                                       #
# Check for proper number of command line args.         #
#                                                       #
#########################################################
echo "$# arguments"
if [ "$#" -lt 10 ] 
        then
        echo ""
        echo ""
        echo ""
        echo "			*****************  Exchange lipid Script   ************************"
        echo "Author: Heidi Koldsø"
	echo "Date: 27 July 2012"
        echo "                          TOO FEW ARGUMENTS"
        echo "                          "
        echo ""
        echo " USAGE:       sh run_exchange_lipids.sh -np <number of proteins> -top <protein.itp> -f <gro-file> -ut <upper-topology> -uc <upper-composition> -lt <lower-topology> -lc <lower composition> "
        echo ""
	echo ""
	echo ""
        echo "v1.0.0"
        echo ""
        echo "Exchange lipid Script"
        echo "Author: Heidi Koldsø"
	echo ""
	echo ""
	echo " DESCRIPTION: This script will transform a POPC (with cholesterol) to a membrane with the specified composition of the upper and lower leaflet "
	echo ""
        echo ""
	echo " USAGE:       sh run_exchange_lipids.sh -np <number of proteins> -top <protein.itp> -f <gro-file> -ut <upper-topology> -uc <upper-composition> -lt <lower-topology> -lc <lower composition> "

        echo ""
	echo "EXAMPLE with protein: sh run_exchange_lipids.sh -np 1 -top protwin-cg.itp -f system.gro -ut POPC:POPE:DPPC -uc 80:10:10 -lt POPC:POPS:POPG:PPCS -lc 50:10:20:20"
	echo ""
 	echo ""
        echo "EXAMPLE without protein: sh run_exchange_lipids.sh -np 0 -top none -f system.gro -ut POPC:POPE:DPPC -uc 80:10:10 -lt POPC:POPS:POPG:PPCS -lc 50:10:20:20"
        echo ""

        echo "Option    	          Description                    "
        echo "-----------------------------------------------------"
        echo "-np               : Number of proteins"
        echo "-top              : Protein itp file [protein-cg.itp]. If no protein type [none]"
        echo "-f                : System gro file "
        echo "-ut               : Topology of upper leaflet  [e.g. POPC:PPCS:DPPC]"       
        echo "-uc               : Composition of upper leaflet  [e.g. 50:20:30]"      
	echo "-lt               : Topology of lower leaflet  [e.g. POPC:POPG:POPS]"       
        echo "-lc               : Composition of lower leaflet  [e.g. 70:10:20]"      
	echo "-lip (optional)	: Lipid input (Can be POPC (default), POPS, POPE, POPG, PVPG, PVPE  "
	echo ""
	echo ""
        echo ""
	echo " LIPIDS currently available: DPPC, DHPC, DLPC, POPC, DUPC, PPCS, DPPE, DHPE, DLPE, DLPG, DPPG, POPE, POPG, POPS, DOPC, DSPC, DAPC, DSPE, DOPE, DOPG, DOPS, GM3, GM3P (uncharged GM3) PIP2, PIP3, CHOL, POPA, PVPG, PVPE, DAG, PVPA, CERA"
	echo ""
        echo " "
  exit 0
fi
#########################################################
#                                                       #
#       Set the selections based on arguments           #
#                                                       #
#########################################################

protnum=$2
protein_itp=$4
grofile=$6
#liplist="$7 $8 $9 ${10} ${11} ${12} ${13} ${14}"


	#################################################################
	#								#
	#	1) Test if both lower and upper leaflet are specified	#
	#	2) Set the lip composition and percentage  		#
	#	3) Substitute : with space				#
	#								#
	#################################################################
	### the ut flag
	if [ "$7" = "-ut" ]
		then
		up_liplist=$8
		 up_liplist=`echo $up_liplist | sed "s|:| |g"`
	fi
	if [ "$9" = "-ut" ]
                then
                up_liplist=${10}
                 up_liplist=`echo $up_liplist | sed "s|:| |g"`
        fi
	if [ "${11}" = "-ut" ]
                then
                up_liplist=${12}
                 up_liplist=`echo $up_liplist | sed "s|:| |g"`
        fi
	if [ "${13}" = "-ut" ]
                then
                up_liplist=${14}
                 up_liplist=`echo $up_liplist | sed "s|:| |g"`
        fi
	if [ "${15}" = "-ut" ]
                then
                up_liplist=${16}
                 up_liplist=`echo $up_liplist | sed "s|:| |g"`
        fi
	##### the uc flag
	if [ "$7" = "-uc" ]
                then
                up_liplist_percentage=$8
                 up_liplist_percentage=`echo $up_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "$9" = "-uc" ]
                then
                up_liplist_percentage=${10}
                 up_liplist_percentage=`echo $up_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "${11}" = "-uc" ]
                then
                up_liplist_percentage=${12}
                 up_liplist_percentage=`echo $up_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "${13}" = "-uc" ]
                then
                up_liplist_percentage=${14}
                 up_liplist_percentage=`echo $up_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "${15}" = "-uc" ]
                then
                up_liplist_percentage=${16}
                 up_liplist_percentage=`echo $up_liplist_percentage | sed "s|:| |g"`
        fi
	### the lt flag
        if [ "$7" = "-lt" ]
                then
                down_liplist=$8
                 down_liplist=`echo $down_liplist | sed "s|:| |g"`
        fi
        if [ "$9" = "-lt" ]
                then
                down_liplist=${10}
                 down_liplist=`echo $down_liplist | sed "s|:| |g"`
        fi
        if [ "${11}" = "-lt" ]
                then
                down_liplist=${12}
                 down_liplist=`echo $down_liplist | sed "s|:| |g"`
        fi
        if [ "${13}" = "-lt" ]
                then
                down_liplist=${14}
                 down_liplist=`echo $down_liplist | sed "s|:| |g"`
        fi
        if [ "${15}" = "-lt" ]
                then
                down_liplist=${16}
                 down_liplist=`echo $down_liplist | sed "s|:| |g"`
        fi
        ##### the uc flag
        if [ "$7" = "-lc" ]
                then
                down_liplist_percentage=$8
                 down_liplist_percentage=`echo $down_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "$9" = "-lc" ]
                then
                down_liplist_percentage=${10}
                 down_liplist_percentage=`echo $down_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "${11}" = "-lc" ]
                then
                down_liplist_percentage=${12}
                 down_liplist_percentage=`echo $down_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "${13}" = "-lc" ]
                then
                down_liplist_percentage=${14}
                 down_liplist_percentage=`echo $down_liplist_percentage | sed "s|:| |g"`
        fi
        if [ "${15}" = "-lc" ]
                then
                down_liplist_percentage=${16}
                 down_liplist_percentage=`echo $down_liplist_percentage | sed "s|:| |g"`
        fi

	###the lip flag
        if [ "$7" = "-lip" ]
                then
                inlip=$8
        fi
        if [ "$9" = "-lip" ]
                then
                inlip=${10}
        fi
        if [ "${11}" = "-lip" ]
                then
                inlip=${12}
        fi
        if [ "${13}" = "-lip" ]
                then
                inlip=${14}
        fi
        if [ "${15}" = "-lip" ]
                then
                inlip=${16}
        fi

#################################
#
#       Make submit.out file with input
#
#############################

if [ -f "exchange_lipids.out" ]; then
        mv exchange_lipids.out exchange_lipids.out.old
fi

echo "input was: sh run_exchange_lipids.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17}" > exchange_lipids.out




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
                        echo "vmddir is: $vmddir"
 fi

 if ([ "$machine" = "Linux" ] && [ "$architecture" = "i686" ])
                        then
                        vmddir=/sbcb/packages/opt/Linux_x86/vmd/1.9/bin/vmd
                        echo "vmddir is: $vmddir"
 fi

 if ([ "$machine" = "Darwin" ] && [ "$architecture" = "x86_64" ])
                        then
                        vmddir=/sbcb/packages/opt/Darwin_x86/vmd/1.9/vmd/vmd_MACOSXX86
                        echo "vmddir is: $vmddir"
 fi

 if ([ "$machine" = "Darwin" ] && [ "$architecture" = "i386" ])
                        then
                        vmddir=/sbcb/packages/opt/Darwin_x86/vmd/1.9/vmd/vmd_MACOSXX86
                        echo "vmddir is: $vmddir"
 fi



		#########################################################
                #                                                       #
                #       Run the tcl script thorugh vmd                  #
                #                                                       #
                #########################################################
 		echo "vmdinput is: $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16}"               
                $vmddir -dispdev text -e /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/exchange_lipids.tcl -args $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16}


		#########################################################
                #                                                       #
                #       Cat the pdb-files from vmd together dependent 	#
		# 	and get number of residues of each type		#
                #                                                       #
                #########################################################

		#################################################
		# Get the complete list of lipid types		#
		################################################
		
		updownlist="$up_liplist $down_liplist"
	
		#################################################
                # Get the complete list of lipid types          #
                ################################################

                tail -1 $grofile > boxdim.dat
		
		#########################################################
		#							#
                # 	If protein is present put the protein in the 	#
		#	start of the pdb-file and put the number of 	#
		#	proteins in the top file. Else use nonprotein 	#
		#	file						#
		#							#
                #########################################################
	#	if [ -f PROTEIN.pdb ]; then
		if [ "$protnum" -gt "0" ]; then 
		cp PROTEIN.pdb exchange_lipids.pdb
		sed -e "s|PROTEIN-ITP|$protein_itp|g;s|PROTNUM|$protnum|g" /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/topol_protein.top > exchange_lipids.top
		fi
		if [ "$protnum" -eq "0" ]; then
#		if [ ! -f PROTEIN.pdb ]; then
		echo "CRYST1  225.221  225.221  110.042  90.00  90.00  90.00 P 1           1" > exchange_lipids.pdb
		cp /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/topol_noprotein.top exchange_lipids.top
		fi
		
                cp -r /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/mdp-files/*.itp .

		#########################################################
                #                                                       #
                #      	For all lipids get number of lipid and print to	# 
		#	topology file and cat the pdb-files together    #
                #       and put the number of proteins in the top file  #
                #                                                       #
                #########################################################
		i=1
#		mv exchange_lipdids.pdb tmp.pdb
		for lip in $updownlist; do
			#########################################################
	                #                                                       #
        	        #  Cat together pdb-files from lipids which have been 	#
			#  aligned instead of exchanged within the tcl script 	# 
        	        #                                                       #
	                #########################################################
		if [ -f ${lip}_0001.pdb ]; then
                     if [ "$lip"  = "DOPC" -o "$lip"  = "DSPC" -o "$lip"  = "DAPC" -o "$lip"  = "DSPE" -o "$lip"  = "DOPE" -o "$lip"  = "DOPG" -o "$lip"  = "DOPS" -o "$lip"  = "GM1" -o "$lip"  = "GM3" -o "$lip"  = "GM3P"  -o "$lip"  = "GM3C" -o "$lip"  = "PIP2" -o "$lip"  = "PIP3" -o "$lip"  = "CHOL" ]; then
			echo "trying to cat"
			cat ${lip}_*.pdb | grep ATOM > ${lip}.pdb
			cp ${lip}.pdb TEST.pdb
			rm ${lip}_*.pdb
		     fi	
		   fi
			#########################################################
                        #                                                       #
                        #  Get number of lipids and add to topology file   	#
			#  1) Phospholipids					#
			#  2) Chol, PIP2, PIP3 and GM3				#
                        #	                                                #
                        #########################################################
		  	if [ -f ${lip}.pdb ] && [ "$lip" != "CHOL" -a "$lip" != "PIP2" -a "$lip" != "PIP3" -a "$lip" != "GM3" -a "$lip" != "GM3P" -a "$lip" != "GM3C" -a "$lip" != "CERA" -a "$lip" != "GM1" ]; then
				echo "putting lipid=$lip"
				lipnum=`grep "PO4 ${lip}" ${lip}.pdb | wc -l`
				echo "$lip $lipnum" >> exchange_lipids.top			
 		  	fi
			#################################################################
                        #                                                               #
                        #       Test if cholesterol is present?                         # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
		  	if [ "$lip"  = "CHOL" ] && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "ROH ${lip}" ${lip}.pdb | wc -l`
                                echo "CHOL $lipnum" >> exchange_lipids.top
		  	fi

			#################################################################
			#                                                               #
                        #       Test if ceramide is present?                         # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
                        if [ "$lip"  = "CERA" ] && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "GL1 ${lip}" ${lip}.pdb | wc -l`
                                echo "CERA $lipnum" >> exchange_lipids.top
                        fi

			#################################################################
                        #                                                               #
                        #       Test if gm3 is present? 	                        # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
			if [ "$lip"  = "GM3" ]  && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "B5  ${lip}" ${lip}.pdb | wc -l`
                                echo "GM3 $lipnum" >> exchange_lipids.top
				mv martini_v2.0_lipids.itp martini_v2.0_lipids_wogm3.itp
				cat martini_v2.0_lipids_wogm3.itp gm3.itp > martini_v2.0_lipids.itp
                  	fi

			#################################################################
                        #                                                               #
                        #       Test if gm3p (uncharged gm3) is present?                                 # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
                        if [ "$lip"  = "GM3P" ]  && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "B5  ${lip}" ${lip}.pdb | wc -l`
                                echo "GM3P $lipnum" >> exchange_lipids.top
				 mv martini_v2.0_lipids.itp martini_v2.0_lipids_wogm3.itp
                                cat martini_v2.0_lipids_wogm3.itp gm3.itp > martini_v2.0_lipids.itp

                        fi

			#################################################################
                        #                                                               #
                        #       Test if gm3c (same as gm3) is present?                                 # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
                        if [ "$lip"  = "GM3C" ]  && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "B5  ${lip}" ${lip}.pdb | wc -l`
                                echo "GM3C $lipnum" >> exchange_lipids.top
				mv martini_v2.0_lipids.itp martini_v2.0_lipids_wogm3.itp
                                cat martini_v2.0_lipids_wogm3.itp gm3.itp > martini_v2.0_lipids.itp

                        fi

			 #################################################################
                        #                                                               #
                        #       Test if gm1 is present?                                 # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
                        if [ "$lip"  = "GM1" ]  && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "AM1 ${lip}" ${lip}.pdb | wc -l`
                                echo "GM1 $lipnum" >> exchange_lipids.top
                                mv martini_v2.0_lipids.itp martini_v2.0_lipids_wogm1.itp
                                cat martini_v2.0_lipids_wogm1.itp gm1.itp > martini_v2.0_lipids.itp
                        fi


			#################################################################
                        #                                                               #
                        #       Test if pip2 is present?                                 # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
			if [ "$lip"  = "PIP2" ]  && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "PO3 PIP" ${lip}.pdb | wc -l`
                                echo "PIP2 $lipnum" >> exchange_lipids.top
                        fi	
			#################################################################
                        #                                                               #
                        #       Test if pip3 is present?                                 # 
                        #       -if yes: count number of molecules and cat together     #
                        #                                                               #
                        #################################################################
			if [ "$lip"  = "PIP3" ]  && [ -f ${lip}.pdb ]; then
                                lipnum=`grep "PO3 PI3" ${lip}.pdb | wc -l`
                                echo "PI3 $lipnum" >> exchange_lipids.top
                        fi
		
		  	### move exchange_lipid file to a tmp file
                       mv exchange_lipids.pdb tmp.pdb

			#########################################################
	                #                                                       #
        	        #       Write out pdb files in correct order for some   #
	                #       of the lipids (POPS, POPG, POPE are okay)       #
	                #                                                       #
	                #########################################################
				#########################################
                                #                                       #
                                #       Test for Ceramide     #
                                #       and write put in correct order  #
                                #########################################
				if [ "$lip"  = "CERA" ]; then
					mv ${lip}.pdb ${lip}_tmp.pdb
                                        end=`cat ${lip}_tmp.pdb | wc -l`
                                        end=$(( $end - 9 ))
                                         k=1
                                              echo "`cat ${lip}_tmp.pdb | head -$(($k )) | tail -$(($k))`" > ${lip}.pdb
                                        while [ $k -lt "$end" ]
                                                do
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 2 )) | tail -2`" >> ${lip}.pdb
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 10 )) | tail -4`" >> ${lip}.pdb
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 6 )) | tail -4`" >> ${lip}.pdb
                                                k=$(($k + 10))

                                        done
					rm ${lip}_tmp.pdb
				fi


				#########################################
				#					#
				#	Test for 8 bead long lipids	#
				#		(DHPE+DHPC)		#
                                #       and write put in correct order  #
				#########################################
		if [ -f ${lip}.pdb ]; then
					if [ "$lip"  = "DHPC" -o "$lip"  = "DHPE" ]; then
        	        	                mv ${lip}.pdb ${lip}_tmp.pdb
                	        	        end=`cat ${lip}_tmp.pdb | wc -l`
                        	        	end=$(( $end - 7 ))
                               			 k=1
	         	                	      echo "`cat ${lip}_tmp.pdb | head -$(($k )) | tail -$(($k))`" > ${lip}.pdb
	        	                        while [ $k -lt "$end" ]
        	                                	do
							echo "`cat ${lip}_tmp.pdb | head -$(($k + 4 )) | tail -4`" >> ${lip}.pdb
							echo "`cat ${lip}_tmp.pdb | head -$(($k + 8 )) | tail -2`" >> ${lip}.pdb
							echo "`cat ${lip}_tmp.pdb | head -$(($k + 6 )) | tail -2`" >> ${lip}.pdb 
	                                                k=$(($k + 8))
							
	                        	        done
						 rm ${lip}_tmp.pdb

					fi
	
				#########################################
        	                #                                       #
                                #       Test for 10 bead long lipids    #
                                #             (DLPE+DLPC+DLPG)          #
                                #       and write put in correct order  #
                                #########################################
				if [ "$lip"  = "DLPC" -o "$lip"  = "DLPE" -o "$lip"  = "DLPG" ]; then
                                        mv ${lip}.pdb ${lip}_tmp.pdb
                                        end=`cat ${lip}_tmp.pdb | wc -l`
                                        end=$(( $end - 9 ))
                                         k=1
                                              echo "`cat ${lip}_tmp.pdb | head -$(($k )) | tail -$(($k))`" > ${lip}.pdb
                                        while [ $k -lt "$end" ]
                                                do
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 4 )) | tail -4`" >> ${lip}.pdb
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 10 )) | tail -3`" >> ${lip}.pdb
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 7 )) | tail -3`" >> ${lip}.pdb
                                                k=$(($k + 10))

                                        done
                                        rm ${lip}_tmp.pdb

                                fi

				#########################################
                                #                                       #
                                #       Test for 12 bead long lipids    #
                                #       (DUPC+PPCS+PPCS+DPPC+DPPE)      #
				# 	and write put in correct order	#
                                #########################################
				if [ "$lip"  = "DUPC" -o "$lip"  = "PPCS" -o "$lip"  = "DPPC" -o "$lip" = "DPPE" -o "$lip" = "DPPG" ]; then
                                        mv ${lip}.pdb ${lip}_tmp.pdb
                                        end=`cat ${lip}_tmp.pdb | wc -l`
                                        end=$(( $end - 11 ))
                                         k=1
                                              echo "`cat ${lip}_tmp.pdb | head -$(($k )) | tail -$(($k))`" > ${lip}.pdb
                                        while [ $k -lt "$end" ]
                                                do
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 4 )) | tail -4`" >> ${lip}.pdb
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 12 )) | tail -4`" >> ${lip}.pdb
                                                echo "`cat ${lip}_tmp.pdb | head -$(($k + 8 )) | tail -4`" >> ${lip}.pdb
                                                k=$(($k + 12))
                                        done
				 rm ${lip}_tmp.pdb
	                 fi
		fi

			#### cat all the pdb-files together
		#	if [ -f ${lip}.pdb ]; then
			cat tmp.pdb ${lip}.pdb > exchange_lipids.pdb

		#	 rm ${lip}.pdb
			##### rm the lipid pdb files to avoid some being writtes twice if they are present both in the upper and lower leaflet
			if [ -f ${lip}.pdb ]
			then
				echo "lip=$lip and i=$i"
				rm ${lip}.pdb
				i=$(($i +1))

			fi	
		done	

		gmx_sse editconf -f exchange_lipids.pdb -c -o exchange_lipids.gro
		line=`cat exchange_lipids.gro | wc -l`
		line=$(($line - 1))
		head -$line exchange_lipids.gro > exchange_lipids_short.gro
		cat exchange_lipids_short.gro boxdim.dat > prot+lipid.gro
		cp /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/vdwradii.dat .
		gmx_sse solvate -cp prot+lipid.gro -cs /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/wat.pdb -radius 0.22 -o prot+lipid+wat.pdb

		rm vdwradii.dat
		#################################################################
                #                                                               #
                #       Count number of water and cat together     		#
                #                                                               #
                #################################################################
		watnum=`grep "W     W"  prot+lipid+wat.pdb | wc -l `
		echo "W $watnum" >> exchange_lipids.top 
			
		cp exchange_lipids.top prot+wat.top
		#########################################################
                #                                                       #
                #      Make a grofile from the pdbfile                  #
                #                                                       #
                #########################################################
	
		gmx_sse editconf -f prot+lipid+wat.pdb -c -o prot+lipid+wat.gro 
		
		#########################################################
                #                                                       #
                #       Copy the em.mdp to this directory and 		#
		#	generate a tpr file               		#
                #                                                       #
                #########################################################

		cp /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/em.mdp .
#		cp -r /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/mdp-files/*.itp .
		gmx_sse grompp -f em.mdp -maxwarn 5 -c prot+lipid+wat.gro -p exchange_lipids.top -o em.tpr
		
		#########################################################
                #                                                       #
                #       Generate random number for ioniziation          #
                #                                                       #
                #########################################################

		RANDOM=$(( 1000+(`od -An -N2 -i /dev/random` )%(9999-100+1)))

		#########################################################
                #                                                       #
                #       Ionize the system to a 0.15 M NaCl		# 
		#	concentration 					#
                #                                                       #
                #########################################################
		### get water group number depending on if protein is present or not
		if [ "$protnum" > 0 ]; then
			watergrp=$((12 + $i))
		 fi
		if [ "$protnum" = 0 ]; then
			watergrp=$((1 + $i))
		fi
		gmx_sse  genion -s em.tpr -seed $RANDOM -o tmp1.gro -neutral -conc 0.15 << EOF
$watergrp
EOF
	
		#########################################################
                #                                                       #
                #       Change the ion residue names to ION             #
                #                                                       #
                #########################################################

		sed -e "s|NA      NA|ION    NA+|g;s|CL      CL|ION    CL-|g" tmp1.gro > exchange_lipids_sys.gro
		
		#########################################################
                #                                                       #
                #       Get the number of NA+ and CL- and the new 	#
		#  	number of waters and exchange the values in 	#
		#	the topology file          			#
                #                                                       #
                #########################################################

		newwat=`grep "W " exchange_lipids_sys.gro | wc -l `
                nanum=`grep "NA+" exchange_lipids_sys.gro | wc -l `
                clnum=`grep "CL-" exchange_lipids_sys.gro | wc -l`
		echo "NA+ $nanum" >> exchange_lipids.top		
		echo "CL- $clnum" >> exchange_lipids.top

		cp  exchange_lipids.top tmp.top
		sed -e "s|W $watnum|W $newwat|g" tmp.top > system.top

		#########################################################
                #                                                       #
                #       Dependent on lipids present, make ndx file with #
		#	the lipids combined and ions and water combined	#
                #                                                       #
                #########################################################
		 if [ "$protnum" > 0 ]; then
			mergegrp=1

		fi
		if [ "$protnum" = 0 ]; then

			mergegrp=rPOPC
                fi

	###########
	#
	#	test if lip=PIP3 because then renaming is necessary:
	#	
	###############################
	if [ "$lip"  = "PIP3" ]; then
		lip=PI3
	fi	
	
	gmx_sse  make_ndx -f exchange_lipids_sys.gro -o sys.ndx << EOF
del 3-30
r${lip}|rPOP*|rPI*|rDHP*|rDPP*|rDMP*|rDOP*|rBOG*|rCHO*|rDDM*|rDSP*|rTOC*|rCAR*|rDLP*|rSQD*|rDGD*|rLMG*|rGM*|rPPCS|rDUPC|rPVP*|rDAP*|rDAG|rCERA|rPCS
rW|rION
name 3 LIPID
name 4 SOL_ION
q
EOF

		#########################################################
                #                                                       #
                #       Make new tpr file and run energy minimzation    #
                #                                                       #
                #########################################################
		gmx_sse  grompp -f em.mdp -maxwarn 5 -c exchange_lipids_sys.gro -n sys.ndx -p system.top -o em.tpr


		gmx_sse  mdrun -v -deffnm exchange_lipid_system -s em.tpr
		info=`grep "inf " exchange_lipid_system.log`
	##### Test if minimization was okay?
		if [ "$info" != "" ]
			then 
			echo ""
			echo "Sorry!!"
			echo ""
	                echo "The minimization was UNSUCCESSFUL!!!"
	                echo ""
	                echo "	Please try again "
	                echo " "
		fi


	if [ "$info" = "" ]
		then
		#cp /sansom/n06/stansfeld/CGDB/scripts/EXCHANGE_LIPIDS/eq.mdp .
		#		grompp -f eq.mdp -po eqout.mdp -c exchange_lipid_system.gro -n sys.ndx -p system.top -o eq.tpr -maxwarn 5
	gmx_sse editconf -f exchange_lipid_system.gro -o exchange_lipid_system.pdb	
#### Clean up
rm  W.pdb tmp.top tmp.pdb tmp1.gro PROTEIN.pdb mdout.mdp genion.log exchange_lipids.gro exchange_lipids.pdb exchange_lipids.top exchange_lipids_sys.gro
		echo ""
		echo "         The scripts is done!"
                echo ""
                echo "The final gro-file is exchange_lipid_system.gro "
		echo "The final pdb-file is exchange_lipid_system.pdb "
		echo "The final top-file is system.top "
                echo " "
                echo "   The final file contains the following composition:"
                echo "		$7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}"
                echo " "

                    echo "Enjoy!"

	fi
