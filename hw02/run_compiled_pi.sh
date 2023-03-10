#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -o pi.compiled.log
#SBATCH -J bash_pi
#SBATCH -p normal
#SBATCH -A cse38018
#SBATCH -t 1:30:00

for i in 10 100 500 1000 5000 10000 50000
do
	./a.out $i
done
