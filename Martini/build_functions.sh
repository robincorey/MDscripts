#!/bin/bash

prot1=${1}
GRO=${2}
DIR=/sansom/s137/bioc1535/Desktop/CG_KIT
lipids=(CHOL POPE POPC POPS)
lipidratios=(67 23 10)
run=yes
CD=`pwd`

prot=${prot1//.pdb}

RebuildGro () {
sed '/^$/d' ${1} -i
head -n 2 ${1} > rebuild_${1}
sed '1s/^.*$/Membrane protein/g' rebuild_${1} -i
tail -n +2 ${1} | grep "BB\|SC" >> rebuild_${1}
for (( i=0; i<${#lipids[@]}; i++ )); do
        tail -n +2 ${1} | grep "${lipids[i]}" >> rebuild_${1}
done
tail -n +2 ${1} | grep "W" >> rebuild_${1}
tail -n +2 ${1} | grep "NA" >> rebuild_${1}
tail -n +2 ${1} | grep "CL-" >> rebuild_${1}
tail -n 1 ${1} >> rebuild_${1}
}

MakeNdx () {
rm -f memsol.ndx index.ndx
echo q | gmx_sse make_ndx -f ${1} -o index >& ndx
cat ndx | grep [0-9] | grep atoms | awk '{print $2}' | sort | uniq -d > ndx2
if [[ -s ndx2 ]] ; then
        while read i; do
                sed '0,/'$i'/! s/'$i'/Dud/' index.ndx -i
        done < ndx2
fi
groups=`grep -c '\[' index.ndx`
memno=$(expr $groups + 0)
solno=$(expr $groups + 1)
lipid_ndxs=`for (( i=0; i<${#lipids[@]}; i++ )); do echo -n '"'${lipids[i]}'"' '|'; done`
echo -e ${lipid_ndxs::-1} '\n' '"'W'"' '|' '"'ION'"' '\n' name $memno Mem '\n' name $solno Sol '\n'q  | gmx_sse make_ndx -f ${1} -o memsol.ndx -n index >& ndx3
}

ReportError () {
if [ ! -f ${1} ]; then
        echo "
        #################################
        Error: no ${1}
        #################################"
        exit 0
else
        echo "${1} done"
fi
}

UpdateToplogy () {
for (( i=0; i<${#lipids[@]}; i++ )); do
        sed "/${lipids[i]}/d" ${prot}.top -i
        lipidnum=`tail -n +2 ${1} | grep ${lipids[i]} | awk '{print $1}' | sort -u | wc -l`
        echo "${lipids[i]} $lipidnum" >> ${prot}.top
done
if grep -q W ${1} ;
        then sed "/^W/d" ${prot}.top -i
        solnum=`tail -n +2 ${1}| grep -c "W "`
        echo "W $solnum" >> ${prot}.top
fi
if grep -q NA ${1} ;
        then sed "/^NA/d" ${prot}.top -i
        NAnum=`tail -n +2 ${1} | grep -c "NA"`
        echo "NA $NAnum" >> ${prot}.top
fi
if grep -q CL- ${1} ;
        then sed "/^CL-/d" ${prot}.top -i
        CLnum=`tail -n +2 ${1} | grep -c "CL-"`
        echo "CL- $CLnum" >> ${prot}.top
fi
}

ChangeIons () {
sed 's/NA+ /ION /g' ${1} -i
sed 's/ NA+/  NA/g' ${1} -i
sed 's/CL- /ION /g' ${1} -i
sed 's/ CL/CL-/g' ${1} -i
}

MakeNdx $GRO
UpdateToplogy $GRO
