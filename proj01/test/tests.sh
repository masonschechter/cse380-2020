#!/usr/bin/env bats

@test "check outputs of 1D linear system initialization" {
	# does it run?
	run ./1D_linear_system_check
	[ "$status" -eq 0 ]
	
	# is the output correct?
	rm -f .tmp_output
	./1D_linear_system_check > .tmp_output
	run diff .tmp_output 1D_check_ref_output
	[ "$status" -eq 0 ]
}

@test "check outputs of 2D linear system initialization" {
	# does it run?
	run ./2D_linear_system_check
	[ "$status" -eq 0 ]

	# is the output correct?
	rm -f .tmp_output
	./2D_linear_system_check > .tmp_output
	run diff .tmp_output 2D_check_ref_output
	[ "$status" -eq 0 ]
}

@test "Jacobi solver test 1D" {

	run ./jacobi_test_1D
	[ "$status" -eq 0 ]
}

@test "Jacobi solver test 2D" {

	run ./jacobi_test_2D
	[ "$status" -eq 0 ]
}

@test "Gauss-Seidel solver test 1D" {

	run ./gauss_seidel_test_1D
	[ "$status" -eq 0 ]
}

@test "Gauss-Seidel solver test 2D" {

	run ./gauss_seidel_test_2D
	[ "$status" -eq 0 ]
}



