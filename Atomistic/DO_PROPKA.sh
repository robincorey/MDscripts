#!/bin/bash
# prep pdb and run propka

SYS1=$1 # input pdb file

if [[ -z $SYS1 ]]; then
	echo "Add pdb as argument you fool"
	exit 0
fi 

SYS=${SYS1//.pdb}

sed 's/HSE/HIS/g' ${SYS}.pdb -i
sed 's/HSD/HIS/g' ${SYS}.pdb -i

rm -f chains

awk '{print $5}' ${SYS}.pdb | uniq | grep -v [0-9] > chains

if [ ! -s chains ]; then
	echo no chain identifiers??
	exit 0
fi

sed '/^$/d' chains -i

readarray ARRAY < chains

# countcolumns
cols=`awk '{print NF; exit}' ${SYS}.pdb`
if [[ $cols -gt 10 ]]; then
	echo there are $cols columns
else
	echo there are $cols columns
fi 

rm -f forpropka*pdb
rm -f preprotonation

for (( i=0; i<${#ARRAY[@]}; i++ )); do
	A=`echo ${ARRAY[i]} | tr -d '\n'`
	#echo "chain $A" > ${A}.pdb
	#while read line; do
#		id=`echo $line | awk '{print $5}'`
#		if [[ $id == ${A} ]]; then
#			echo $line >> ${A}.pdb
#		fi
#	done < ${SYS}.pdb
	convpdb.pl -chain ${A} ${SYS}.pdb > ${A}.pdb
	sed 's/HSE/HIS/g' ${A}.pdb -i
	sed 's/HSD/HIS/g' ${A}.pdb -i
	pdb2gmx -f ${A}.pdb -q forpropka_${A}.pdb -ignh -ff oplsaa -water spc >& chain1
	rm -f forpropka_${A}.pka
	propka31 -f forpropka_${A}.pdb >& prop1
	if [ -f forpropka_${A}.pka ]; then
		echo propka ${A} worked!
	fi
	grep "   LYS" forpropka_${A}.pka > lys
	grep "   ARG" forpropka_${A}.pka > arg
	grep "   ASP" forpropka_${A}.pka > asp
	grep "   GLU" forpropka_${A}.pka > glu
	grep "   HIS" forpropka_${A}.pka > his
	cat lys arg asp glu his > preprotonation_${A}
	cat preprotonation_${A} >> preprotonation
	cat forpropka_${A}.pdb >> forpropka.pdb
	echo hello
done
echo goodbye
count1=`grep -c "CA  LYS\|CA  ARG\|CA  ASP\|CA  GLU\|CA  HIS" forpropka.pdb`
count2=`grep -c "LYS\|ARG\|ASP\|GLU\|HIS" preprotonation`

if [[ $count1 = $count2 ]];
then
	echo atom count ok!
elif [[ $count1 != $count2 ]];
then
	echo not ok!
	exit 0
fi

rm -f protonation

while read line
do
	pka=`echo $line | awk '{print $4}'`
	#pka1=`echo $pka | cut -f1 -d"."`
	if [[ ${#pka} == 5 ]] # i.e. double digits 
	then
		protonationstate=1
	elif [[ "$pka" > "7" ]]
	then
		protonationstate=1
	elif [[ "$pka" < 7 ]]
	then
		protonationstate=0
	fi
	echo $protonationstate >> protonation
done < preprotonation

rm -f ${SYS//.pdb}.top

sed '/OXT/d' forpropka.pdb -i

pdb2gmx -f ${SYS} -o ${SYS//.pdb}_postpdb2gmx -p ${SYS//.pdb} -i ${SYS//.pdb} -ignh -ff oplsaa -water spc -lys -arg -asp -glu -his < protonation >& pdboutput

if [[ ! -f ${SYS//.pdb}_postpdb2gmx.gro ]]
then
	echo "pdb2gmx failed! Read pdboutput"
	exit 0
else
	cp preprotonation Titrated.txt
	echo propka assignment worked. Check values in Titrated.txt 
fi

