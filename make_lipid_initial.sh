#!/bin/bash

echo "[ moleculetype ]
PC12       1
"

for component in atoms 
do
	echo "
[ $component ]"
	for i in {0..11}
	do
		echo "; sugar $i"
		awk -v var=$i '{printf ($1+(var*8)"\t"$2"\t"var+1"\t"$4"\t"$5"\t"$6+(var*8)"\t"$7"\t"$8"\n")}' ${component}_in.itp
	done
done

for component in bonds constraints 
do
        echo "
[ $component ]"
        for i in {0..11}
        do
                awk -v var=$i '{printf ($1+(var*8)"\t"$2+(var*8)"\t"$3"\t"$4"\t"$5"\n")}' ${component}_in.itp
        done
done

for component in angles 
do
        echo "
[ $component ]"
        for i in {0..11}
        do
                awk -v var=$i '{printf ($1+(var*8)"\t"$2+(var*8)"\t"$3+(var*8)"\t"$4"\t"$5"\t"$6"\n")}' ${component}_in.itp
        done
done

for component in PE_angles mannan_angles
do
        echo "
;; extra angles for $component"
        for i in {0..11}
        do
                awk -v var=$i '{printf ($1+(var*8)"\t"$2+(var*8)"\t"$3+(var*8)"\t"$4"\t"$5"\t"$6"\n")}' ${component}_in.itp
        done
done

for component in dihedrals
do
        echo "
[ $component ]"
        for i in {0..11}
        do
                awk -v var=$i '{printf ($1+(var*8)"\t"$2+(var*8)"\t"$3+(var*8)"\t"$4+(var*8)"\t"$5"\t"$6"\t"$7"\t"$8"\n")}' ${component}_in.itp
        done
done
