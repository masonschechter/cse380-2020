#!/bin/bash

#SBATCH -J asymptotic_analysis
#SBATCH -o job.out
#SBATCH -e job.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A cse38018


for i in 10 25 50 100 250 500 1000 2500 5000 10000 25000 50000 100000 250000 500000 1000000 2500000 5000000 10000000
do	
	./integrate simpson $i >> simpson.out
	./integrate trapezoid $i >> trapezoid.out
done
