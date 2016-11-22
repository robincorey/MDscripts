#!/bin/bash 
# Cards selected from pymol and saved 
## Very very specific to my system  
ARRAY=(817 818 819 821 822 823 824 826 829 831 834 835 837 840 842 843 845 846 847 850 852 853 854 855 856 1221 1222 1224 1225 1226 1230 1232 1233 1234 1235 1236 1237 1239 1241 1242 1244 1247 1249 1250 1251 1252 1253 1254 1255 1256 1257 1259 1261)  
rm delcard.gro  
cp prod_micro_10033.gro delcard.gro  
for i in `seq 0 52`; 	
	do sed "/${ARRAY[i]}CARD/d" delcard.gro -i 
done 
