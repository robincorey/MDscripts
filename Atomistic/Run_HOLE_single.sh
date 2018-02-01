#!/bin/bash
# do HOLE

NAME=${1//.pdb}

echo "
! Comment
coord ${NAME}.pdb                  ! Co-ordinates in pdb format
radius ${MM}/simple.rad   ! Use simple AMBER vdw radii
! now optional cards
sphpdb hole_out_${NAME}.sph     ! pdb format output of hole sphere centre info
" > ${NAME}.inp

hole < ${NAME}.inp > ${NAME}.out

if [[ ! -f hole_out_${NAME}.sph ]]; then
	echo hole failed
	exit 0
fi

# make dots
sph_process -dotden 15 -color hole_out_${NAME}.sph dotsurface_${NAME}.qpt
echo -e D '\n' dotsurface_${NAME}.qpt '\n' dotsurface_${NAME}2.qpt | qpt_conv 

# make solid
sph_process -sos -dotden 15 -color hole_out_${NAME}.sph solid_${NAME}.sos
sos_triangle -s < solid_${NAME}.sos > solid_${NAME}.vmd_plot

# Get data
grep sampled ${NAME}.out | awk '{print $1","$2}' > ${NAME}.data

# Plot
$MM/PLOT_hole.py ${NAME}.data

# Make image
echo "# load data
mol load pdb ${NAME}.pdb
#mol delrep 0 top

# color accordingly
mol selection "resid 1 to 756"
mol addrep 0
mol modstyle 1 0 Trace
mol modcolor 1 0 ColorID 15

mol selection "resid 757 to 1174"
mol addrep 0
mol modstyle 2 0 Trace
mol modcolor 2 0 ColorID 9

mol selection "resid 1175 to 1230"
mol addrep 0
mol modstyle 3 0 Trace
mol modcolor 3 0 ColorID 32

mol selection "resid 1231 to 1304"
mol addrep 0
mol modstyle 4 0 Trace
mol modcolor 4 0 ColorID 0

# Hole stuff
source solid_${NAME}.vmd_plot

mol delrep 0 top

#orient
rotate x by 90
rotate z by 90
rotate y by -10

#make image
render TachyonInternal ${NAME}_view.tga

quit
" > make_hole.scr

vmd -e make_hole.scr -dispdev text

convert ${NAME}_view.tga ${NAME}_view.png 
