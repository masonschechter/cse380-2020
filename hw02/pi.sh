#!/bin/bash
START=$(date +%s)

function RandomDraw {
	echo "scale=15; $RANDOM/32767" | bc -l
}

pi=`echo "4*a(1)" | bc -l`

num_loops=$1
ctr_inside=0
ctr_outside=0

for (( i=1; i<=$num_loops; i++ ))
do
	x=`RandomDraw`
	y=`RandomDraw`
	x_sqr=`echo $x*$x | bc -l`
	y_sqr=`echo $y*$y | bc -l`
	dist=`echo "$x_sqr + $y_sqr" | bc -l`
	if (( $(echo "$dist<=1" | bc -l) ))
	then
		((ctr_inside++))
	else
		((ctr_outside++))
	fi
done

pi_est=`echo "4 * $ctr_inside / $num_loops" | bc -l`
diff=`echo "$pi - $pi_est" | bc -l`
diff=`echo ${diff#-}`
error=`echo "$diff / $pi" | bc -l`
END=$(date +%s)
time_diff=$(( $END - $START ))
echo $num_loops $ctr_inside $ctr_outside $pi_est $error $time_diff
