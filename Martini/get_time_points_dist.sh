#!/bin/bash

DATA=$1
NAME=${DATA::-10}

rm -f time_points.xvg

init=`printf "%0.2f\n" $(grep -v \# pull_dir_pullx.xvg | grep -v \@ | awk '{print $2}' | head -n 1)`
end=`printf "%0.2f\n" $(grep -v \# pull_dir_pullx.xvg | grep -v \@ | awk '{print $2}' | tail -n 1)`

#init="-0.56" # 0.2 less than end
#end="-0.86" # make sure this is 0.1 less than init in pull

for i in $(seq $init -0.1 $end)
do
	time=`grep "[[:space:]]$i" $DATA | sort -k2 | head -n 1 | awk '{print $1}'`
	time2=`echo "$time - 5" | bc`
	dist=`grep "[[:space:]]$i" $DATA | sort -k2 | head -n 1 | awk '{print $2}'`
	echo $dist $time >> time_points.xvg
	value=`echo "scale = 0; $dist * 100" | bc | cut -f1 -d'.'`
	echo -e Protein '\n' System | gmx_avx trjconv -f $NAME.xtc -s $NAME.tpr -dump $time -b $time2 -o out_$value.gro -pbc mol -center
	gmx_avx editconf -f out_$value.gro -o out_$value.pdb
	gmx_avx grompp -f ../cg-umbrella-dir.mdp -c out_$value.gro -p system.top -o PMF_${value}.tpr -n pull.ndx -maxwarn 2
done
