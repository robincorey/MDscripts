#!/bin/bash
# this is the magic combo for editing CG traj for analysable data

TPR=unbiased
XTC=192_small

## firstly, to prevent protein from jumping
echo System | trjconv -s $TPR -f ${XTC}.xtc -pbc nojump -o ${XTC}.nojump.xtc 

## secondly, fit to SecY w/regards rot and trans in x and y direction
echo -e a_1843-2717 '\n' System | trjconv -s $TPR -f ${XTC}.nojump.xtc -fit rotxy+transxy -o ${XTC}.rot.xtc -n chains

## then, use mol whole to keep square shape
echo System | trjconv -s $TPR -f ${XTC}.rot.xtc -pbc mol -o ${XTC}.pbc.xtc 

