#!/bin/bash

DATA=$1
NAME=${DATA::-10}

rm -f time_points.xvg

for i in {47..447..10}
do
	dist=`printf '%.3f\n' "$(echo "scale=3;${i}/100" | bc)"` 
	time=`grep "[[:space:]]$dist" $DATA | sort -k2 | head -n 1 | awk '{print $1}'`
	time2=`echo "$time - 5" | bc`
	actual_dist=`grep "[[:space:]]$dist" $DATA | sort -k2 | head -n 1 | awk '{print $2}'`
	echo $actual_dist $time >> time_points.xvg
	echo -e Protein '\n' System | gmx_avx trjconv -f $NAME.xtc -s $NAME.tpr -dump $time -b $time2 -o $dist.gro -pbc mol -center
	gmx_avx editconf -f $dist.gro -o $dist.pdb
	gmx grompp -f ../cg-pull-dir.mdp -c 3.47.gro -p ../V*top -o del -n memsol_pull.ndx -maxwarn 2
done
