#!/bin/bash

for i in `seq 47500 22500 250000`
do
	echo -n "$i"
	for j in {1..5}
	do
		value=`awk -v var=$i '{if ($1 == var) print $2}' run_${j}_eq_f/dF_t.txt`
		echo -n ",$value"
	done
	printf "\n"
done

<<'END'
   47500.0     139.554 +- 0.196            25000.0     142.943 +- 0.050
   70000.0     143.655 +- 0.114            47500.0     143.090 +- 0.052
   92500.0     142.462 +- 0.097            70000.0     142.765 +- 0.056
  115000.0     141.673 +- 0.084            92500.0     143.025 +- 0.059
  137500.0     140.914 +- 0.077           115000.0     143.452 +- 0.063
  160000.0     140.478 +- 0.071           137500.0     144.158 +- 0.068
  182500.0     141.341 +- 0.063           160000.0     144.990 +- 0.075
  205000.0     142.537 +- 0.056           182500.0     145.445 +- 0.087
  227500.0     142.691 +- 0.053           205000.0     144.468 +- 0.113
  250000.0     142.943 +- 0.050
END
