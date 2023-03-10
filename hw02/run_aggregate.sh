#!/bin/bash

#SBATCH -N 1
#SBATCH -n 48
#SBATCH -o iter.log
#SBATCH -e iter.err
#SBATCH -J iter_run
#SBATCH -p skx-normal
#SBATCH -A cse38018
#SBATCH -t 03:59:99

module purge
module load TACC
module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=./iter_job_file.sh

START=$(date +"%s")
pi=`echo "scale=15; 4*a(1)" | bc -l`
num_iter=1
error=1
error_goal=`echo '.000005'`
echo "#iter" "num_samples" "num_i" "pi" "relative_error" "time_accum"

while (( $(echo "$error >= $error_goal" | bc -l) ))
do
	$LAUNCHER_DIR/paramrun >> garbage.dump
	agg=`awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}' iter.temp`
	num_samples=`echo $agg | awk '{ print $1 }'`
	num_i=`echo $agg | awk '{ print $2 }'`
	num_o=`echo $agg | awk '{ print $3 }'`
	pi_est=`echo "4 * $num_i / $num_samples" | bc -l`
	diff=`echo "$pi - $pi_est" | bc -l`
	diff=`echo ${diff#-}`
	error=`echo "$diff / $pi" | bc -l`
	STOP=$(date +"%s")
	time_diff=`echo "$STOP - $START" | bc -l`
	echo $num_iter $num_samples $num_i $pi_est $error $time_diff
	num_iter=$(( $num_iter + 1 ))
done
