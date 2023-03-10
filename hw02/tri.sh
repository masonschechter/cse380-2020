#!/bin/bash


total_loops=5

for i in {1..5}
do
	declare -a my_array=()
	for (( j=1; j <= $i; j++ ))
	do
		my_array+=($i)
	done
	total_spaces=$(( $total_loops - $j + 1))
	spaces=""
	for (( k=1; k<=$total_spaces; k++ ))
	do
		spaces+=" "
	done
	echo -n "$spaces"
	for var in "${my_array[@]}"
	do
		echo -n $var' '
	done
	echo ""
done

for i in {1..5}
do
        declare -a my_array=()
        for (( j=1; j <= $i; j++ ))
        do
                my_array+=('.')
        done
        total_spaces=$(( $total_loops - $j + 1))
        spaces=""
        for (( k=1; k<=$total_spaces; k++ ))
        do
                spaces+=" "
        done
        echo -n "$spaces"
        for var in "${my_array[@]}"
        do
                echo -n $var' '
        done
        echo ""
done

