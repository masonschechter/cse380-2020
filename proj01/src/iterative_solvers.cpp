#include <iostream>
#include <stdlib.h>

#include <helper_functions.cpp>

using namespace std;

double* Jacobi(int discretization_points, double** coefficient_matrix, double* T_old, double* RHS, double error_threshold, int max_iterations, string OUTPUT_MODE) {
	double* T_new = new double[discretization_points];
	if ( !T_new ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for temperature vector \n");
		exit(1);
	}
	double initValue = 0.0;
	fill_n(T_new, discretization_points, initValue); // Initialize all values to 0.0
	double row_sum;
	double* swap_ptr;
	int iter = 1; // counter
	while ( L2_norm(discretization_points, T_old, T_new) > error_threshold && iter < max_iterations) {
		if (OUTPUT_MODE == "DEBUG") {
			// Print iteration number, T_old and T_new
			grvy_printf(GRVY_DEBUG, "Iteration number: %d \n", iter);
			grvy_printf(GRVY_DEBUG, "T_old vector: \n");
			for (int i = 0; i < discretization_points; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", T_old[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
			grvy_printf(GRVY_DEBUG, "T_new vector: \n");
			for (int i = 0; i < discretization_points; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", T_new[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
		}
		if (iter++ > 1) { // Edge case for initial iteration. Each subsequent iteration we increment then swap the pointers of T_old and T_new

			swap_ptr = T_old;
			T_old = T_new;
			T_new = swap_ptr;
		}
		// Solves 1 iteration
		for (int i = 0; i < discretization_points; i++) {
			row_sum = 0.0; // reset value for each row
			for (int j = 0; j < discretization_points; j++) {
				if (i != j) { // all of the terms except for our current T_i
					row_sum += (coefficient_matrix[i][j] * T_old[j]);
				}
			} 
			T_new[i] = (RHS[i] - row_sum) / coefficient_matrix[i][i];
		}
	}
	grvy_printf(GRVY_INFO, "This solution took %d iterations \n", iter);
	return T_new;
}

double* Gauss_Seidel(int discretization_points, double** coefficient_matrix, double* T_old, double* RHS, double error_threshold, int max_iterations, string OUTPUT_MODE) {
	double* T_new = new double[discretization_points];
	if ( !T_new ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for temperature vector \n");
		exit(1);
	}
	double initValue = 0.0;
	fill_n(T_new, discretization_points, initValue); // Initialize all values to 0.0
	double row_sum;
	int iter = 1;
	while ( L2_norm(discretization_points, T_old, T_new) > error_threshold && iter < max_iterations) {
		if (OUTPUT_MODE == "DEBUG") {
			// Print iteration number, T_old and T_new
			grvy_printf(GRVY_DEBUG, "Iteration number: %d \n", iter);
			grvy_printf(GRVY_DEBUG, "T_old vector: \n");
			for (int i = 0; i < discretization_points; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", T_old[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
			grvy_printf(GRVY_DEBUG, "T_new vector: \n");
			for (int i = 0; i < discretization_points; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", T_new[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
		}
		if (iter++ > 1) { // Edge case for initial iteration. Each subsequent iteration we increment then copy T_new to T_old
			memcpy(T_old, T_new, sizeof(double) * discretization_points); // T_old = T_new, saving values for next iteration
		} else {
			memcpy(T_new, T_old, sizeof(double) * discretization_points); // T_new = T_old for first iteration
		}
		// Solves 1 iteration
		for (int i = 0; i < discretization_points; i++) {
			row_sum = 0.0; // reset value for each row
			for (int j = 0; j < discretization_points; j++) {
				if (i != j) { // all of the terms except for our current T_i
					row_sum += (coefficient_matrix[i][j] * T_new[j]);
				}
			}
			T_new[i] = (RHS[i] - row_sum) / coefficient_matrix[i][i];
		}
	}
	grvy_printf(GRVY_INFO, "This solution took %d iterations \n", iter);
	return T_new;
}