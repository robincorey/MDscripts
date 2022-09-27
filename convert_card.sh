#!/bin/bash

GRO=$1
#  script to convert CARD from one file format another
#  highly format specific

if [ -z $GRO ]; then echo no gro; exit 0; fi

line=`grep "CARD   GL0" $GRO | sed "s/GL0/Dum/g"`
x=`grep "CARD   GL0" $GRO | awk '{print $4}'`
y=`grep "CARD   GL0" $GRO | awk '{print $5}'`
z=`grep "CARD   GL0" $GRO | awk '{print $6}'`
num=`head -n 2 $GRO | tail -n 1`
newnum=`echo "$num + 1" | bc`

echo $x $y $z

xn1=`echo "$x + 0.5" | bc`
xn2=`echo "$x - 0.5" | bc`

cp $GRO backup_$GRO

sed "0,/CARD   GL0/{s/ $x / $xn1 /}" $GRO -i
sed "/CARD   C5B/a \ $line" $GRO -i
sed "s/\ $line/$line/g" $GRO -i
sed "0,/CARD   Dum/{s/ $x / $xn2 /}" $GRO -i
#sed "0,/ [^0-9]*([0-9]+).CARD   Dum/{s/CARD   Dum/}" $GRO -i
sed "s/^$num/$newnum/g" $GRO -i

echo q | gmx_avx make_ndx -f $GRO -o index > ndx
cat ndx | grep [0-9] | grep atoms | awk '{print $2}' | sort | uniq -d > ndx2
if [[ -s ndx2 ]] ; then
        while read i; do
                sed '0,/'$i'/! s/'$i'/Dud/' index.ndx -i
        done < ndx2
fi

groups=`grep -c '\[' index.ndx`

mem=$(expr $groups + 0)
sol=$(expr $groups + 1)

echo -e "!aW & !aNA+ & !aCL-" '\n' '"'W'"' '|' '"'ION'"' '\n' name $mem LIPID '\n' name $sol SOL_ION '\n'q  | gmx_avx make_ndx -f $GRO -o sys.ndx -n index.ndx
