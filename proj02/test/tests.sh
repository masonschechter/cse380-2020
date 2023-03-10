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

@test "Jacobi solver test 1D, 2nd order" {
	rm -f jacobi_1D_2nd_order_test_out.h5
	run ./jacobi_test_1D
	h5diff -d 1e-6 jacobi_1D_2nd_order_test_out.h5 jacobi_1D_2nd_order_ref_out.h5 /numerical_solution
	[ "$status" -eq 0 ]
}

@test "Jacobi solver test 2D, 2nd order" {
	rm -f jacobi_2D_2nd_order_test_out.h5
	run ./jacobi_test_2D
	h5diff -d 1e-6 jacobi_2D_2nd_order_test_out.h5 jacobi_2D_2nd_order_ref_out.h5 /numerical_solution
	[ "$status" -eq 0 ]
}

@test "Gauss-Seidel solver test 1D, 2nd order" {
	rm -f gauss_seidel_2D_2nd_order_test_out.h5
	run ./gauss_seidel_test_1D_2nd_order
	h5diff -d 1e-6 gauss_seidel_1D_2nd_order_test_out.h5 gauss_seidel_1D_2nd_order_ref_out.h5 /numerical_solution
	[ "$status" -eq 0 ]
}

@test "Gauss-Seidel solver test 1D, 4th order" {
	rm -f gauss_seidel_1D_4th_order_test_out.h5
	run ./gauss_seidel_test_1D_4th_order
	h5diff -d 1e-6 gauss_seidel_1D_4th_order_test_out.h5 gauss_seidel_1D_4th_order_ref_out.h5 /numerical_solution
	[ "$status" -eq 0 ]
}

@test "Gauss-Seidel solver test 2D, 2nd order" {
	rm -f gauss_seidel_2D_2nd_order_test_out.h5
	run ./gauss_seidel_test_2D_2nd_order
	h5diff -d 1e-6 gauss_seidel_2D_2nd_order_test_out.h5 gauss_seidel_2D_2nd_order_ref_out.h5 /numerical_solution
	[ "$status" -eq 0 ]
}

@test "Gauss-Seidel solver test 2D, 4th order" {
	rm -f gauss_seidel_2D_4th_order_test_out.h5
	run ./gauss_seidel_test_2D_4th_order
	h5diff -d 1e-6 gauss_seidel_2D_4th_order_test_out.h5 gauss_seidel_2D_4th_order_ref_out.h5 /numerical_solution
	[ "$status" -eq 0 ]
}


