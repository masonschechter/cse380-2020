#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace MASA;


double** make_1D_coefficient_matrix(int discretization_points, int order) {
	double initValue = 0.0;
	double** matrix = new double*[discretization_points];
	if ( !matrix ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for coefficient matrix \n");
		exit(1);
	}
	// Initialize each element in matrix to 0
	for (int i = 0; i < discretization_points; i++) {
		matrix[i] = new double[discretization_points];
		if ( !matrix[i] ) { // check for NULL pointer
			grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for a row in coefficient matrix \n");
			exit(1);
		}
		fill_n(matrix[i], discretization_points, initValue); // Initialize all values to 0.0
	}
	if (order == 2){
		for (int i = 0; i < discretization_points; i++) {
			// Initializes rows for the boundaries
			if (i == 0 || i == discretization_points - 1) {
				matrix[i][i] = 1;
			// Initializes rows for interior points
			} else {
				matrix[i][i - 1] = 1;
				matrix[i][i] = -2;
				matrix[i][i + 1] = 1;
			}
		}
	} else if (order == 4) {
        for (int i = 0; i < discretization_points; i++) {
      	    // Initializes rows for the boundaries
            if (i == 0 || i == 1 || i == discretization_points - 2 || i == discretization_points - 1) {
                matrix[i][i] = 1;
			// Initializes rows for interior points
			} else {
				matrix[i][i - 2] = -1.0;
				matrix[i][i - 1] =  16.0;
                matrix[i][i]     = -30.0;
                matrix[i][i + 1] =  16.0;
				matrix[i][i + 2] = -1.0;
			}
		}
	} 
return matrix;
}

double* initialize_T_vector_1D(int discretization_points, int order, double X_MIN, double X_MAX) {
	double* T_vector = new double[discretization_points];
	if ( !T_vector ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for temperature vector \n");
		exit(1);
	}
	double initValue = 0.0;
	double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);
	double pos;
	fill_n(T_vector, discretization_points, initValue); // Initialize all values to 0.0
	if (order == 2) {
		T_vector[0] = masa_eval_1d_exact_t(X_MIN); // x = 0
		T_vector[discretization_points - 1] = masa_eval_1d_exact_t(X_MAX); // x = L
	} else if (order == 4) {
		T_vector[0] = masa_eval_1d_exact_t(X_MIN); // x = 0
		T_vector[discretization_points - 1] = masa_eval_1d_exact_t(X_MAX); // x = L
		pos = X_MIN + delta_x; // "interior boundary point", left side
		T_vector[1] = masa_eval_1d_exact_t(pos); 
		pos = (discretization_points - 2) * delta_x; // "interior boundary point", right side
		T_vector[discretization_points - 2] = masa_eval_1d_exact_t(pos); 
	}
	return T_vector;
}

double* initialize_RHS_1D(int discretization_points, int order, double X_MIN, double X_MAX, double k_0) {
	double* RHS = new double[discretization_points];
	if ( !RHS ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for right-hand side vector \n");
		exit(1);
	}
	double q;
	double pos;
	double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);

	if (order == 2) {
		for (int i = 0; i < discretization_points; i++) {
			if (i == 0) { // x = 0
				RHS[i] = masa_eval_1d_exact_t(X_MIN);
			} else if (i == discretization_points - 1) { // x = L
				RHS[i] = masa_eval_1d_exact_t(X_MAX);
			} else {
				pos = i * delta_x;
				q = masa_eval_1d_source_t(pos);
				RHS[i] = -(q * delta_x * delta_x) / k_0;	
			}
		}
	} else if (order == 4) {
		for (int i = 0; i < discretization_points; i++) {
			if (i == 0) { // x = 0
				RHS[i] = masa_eval_1d_exact_t(X_MIN);
			} else if (i == 1) { // "interior boundary point", left side
				pos = i * delta_x;
				RHS[i] = masa_eval_1d_exact_t(pos);
			} else if (i == discretization_points - 2) { // "interior boundary point", right side
				pos = i * delta_x;
				RHS[i] = masa_eval_1d_exact_t(pos);
			} else if (i == discretization_points - 1) { // x = L
				RHS[i] = masa_eval_1d_exact_t(X_MAX);
			} else {
				pos = i * delta_x;
				q = masa_eval_1d_source_t(pos);
				RHS[i] = -(q * delta_x * delta_x * 12.0) / k_0;
			}
		}
	}
	return RHS;
}