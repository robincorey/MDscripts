#!/bin/bash

for i in {1..5}
do
	tail -n 1 run_${i}_eq_f/results.txt 
done
