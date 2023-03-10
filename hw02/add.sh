#!/bin/bash

first_int=$1
second_int=$2

if [[ $# -ne 2 ]]; then
	echo "This program will only run with exactly 2 integer arguments. Exiting..."
	exit 2
fi

sum=$(( first_int + $second_int ))
echo 'The sum of ' $first_int ' and ' $second_int ' is ' $sum
