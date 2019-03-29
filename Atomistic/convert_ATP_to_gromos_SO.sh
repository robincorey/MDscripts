#!/bin/bash

NAME=$1

cp $1 ${1}_ATP.pdb

sed 's/HETATM/ATOM /g' ${1}_ATP.pdb -i
sed "s/ O1/O1P/" ${1}_ATP.pdb -i
sed "s/ O2/O2P/" ${1}_ATP.pdb -i
sed "s/ O3/O3P/" ${1}_ATP.pdb -i
sed "s/1H5*/H61/"  ${1}_ATP.pdb -i
sed "s/2H5*/H62/"  ${1}_ATP.pdb -i
sed "/ H[0-9]\*/d" ${1}_ATP.pdb -i
sed "s/O2P\`/ O2*/g" ${1}_ATP.pdb -i
sed "s/O3P\`/ O3*/g" ${1}_ATP.pdb -i
sed "s/\`/*/g" ${1}_ATP.pdb -i

# reorder

sed 's/*/@/g'  ${1}_ATP.pdb -i

ATOMS=(N9 C4 N3 C2 N1 C6 N6 H61 H62 C5 N7 C8 C1@ O4@ C4@ C2@ O2@ H2@ C3@ O3@ H3@ C5@ O5@ PA O1PA O2PA O3PA PB O1PB O2PB O3PB PG O1PG O2PG O3PG)

rm -f new_ATP.pdb

for i in ${ATOMS[@]}
do
	if grep -q $i ${1}_ATP.pdb
	then
		grep " $i " ${1}_ATP.pdb >> new_ATP.pdb
	else
		echo no $i >> new_ATP.pdb
	fi
done

sed 's/@/*/g' new_ATP.pdb -i
