# makin' mutated itps

ARRAY=(Sec*E.pdb)
ARRAYLEN=${#ARRAY[@]}
ARRAY2=`echo ${ARRAYLEN} - 1 | bc`
CURRENT_DIR=`pwd`

for i in `seq 0 ${ARRAY2}`; do
	cd ${CURRENT_DIR}
	if [[ ! -d ${ARRAY[i]//.pdb} ]]; then
		mkdir ${ARRAY[i]//.pdb}
	fi
	cd ${ARRAY[i]//.pdb}
	python ../martinize.py -f ../${ARRAY[i]} -o ${ARRAY[i]//.pdb}_orig.top -p ${ARRAY[i]//.pdb}-posres.itp -x ${ARRAY[i]//.pdb}-CG.pdb -dssp mkdssp -p backbone -ff martini22 -elastic -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0
	cp /home/birac/Desktop/CoarseGrainSims/3DIN/Elastic/Run2/Run3/Run4/eq4.gro ${ARRAY[i]//.pdb}.gro
	cp /home/birac/Desktop/CoarseGrainSims/3DIN/Elastic/Run2/Run3/Run4/3DIN_carddel3.top ${ARRAY[i]//.pdb}.top
	itppresent=(*.itp)
	itpavail=(/home/birac/Desktop/CoarseGrainSims/3DIN/Elastic/Run2/Run3/Run4/*.itp)
	itplen=`echo ${#itpavail[@]} - 1 | bc`
	for j in `seq 0 ${itplen}`; do 
		#echo ${itpavail[j]}
		itpavailname=${itpavail[j]##*/}
		#echo $itpavailname
		if [[ ! ${itpavailname} == ${itppresent[*]} ]]; then
			cp ${itpavail[j]} .
		fi
	done
done
