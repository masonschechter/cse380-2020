#!/bin/bash

input=$1

if [[ $# -ne 1 ]]; then
	echo "This program will only run with exactly 1 integer arguments. Exiting..."
	exit 2
fi

factorial=1
current_num=$input
while [[ $current_num -gt 1 ]]
do
	factorial=`expr $factorial \* $current_num`
	current_num=`expr $current_num - 1`
done
echo $input ' factorial is ' $factorial
