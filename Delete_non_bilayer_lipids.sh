#!/bin/bash

GRO=$1
TPR=$2

gro=${GRO%.*}

# edit this for lipids and axis
echo 'delete = (resname POPE or resname POPC or resname CHOL) and z < 4 ;
name "*" and not delete' > delete_lipids.dat

gmx_sse select -f $GRO -s $TPR -sf delete_lipids.dat -selrpos whole_res_com -on delete_these
gmx_sse editconf -f $GRO -o ${gro}_trimmed.gro -n delete_these.ndx
$CG/build_functions.sh $GRO
