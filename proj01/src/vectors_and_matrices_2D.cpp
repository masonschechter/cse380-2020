#include <iostream>
#include <stdlib.h>


using namespace std;
using namespace MASA;

double** make_2D_coefficient_matrix(int discretization_points, int order) {
	double initValue = 0.0;
	int total_points = (discretization_points * discretization_points);
	double** matrix = new double*[total_points];
	if ( !matrix ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for coefficient matrix \n");
		exit(1);
	}
	// Initialize each element in matrix to 0
	for (int i = 0; i < total_points; i++) {
		matrix[i] = new double[total_points];
		if ( !matrix[i] ) { // check for NULL pointer
			grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for a row in coefficient matrix \n");
			exit(1);
		}
		fill_n(matrix[i], total_points, initValue);
	}
	if (order == 2) {
		for (int i = 0; i < total_points; i++) {
			if (i < discretization_points) { // top boundary
				matrix[i][i] = 1;
			} else if (i > total_points - discretization_points - 1) { // bottom boundary
				matrix[i][i] = 1;
			} else if (i % discretization_points == 0) { // left boundary
				matrix[i][i] = 1;
			} else if (i % discretization_points == discretization_points - 1) { // right boundary
				matrix[i][i] = 1;
			} else { // interior points
				matrix[i][i - discretization_points] = 1;
				matrix[i][i - 1] = 1;
				matrix[i][i]     = -4;
				matrix[i][i + 1] = 1;
				matrix[i][i + discretization_points] = 1;
			}
		}
	} else if (order == 4) {
		for (int i = 0; i < total_points; i++) {
			if (i < 2 * discretization_points) { // top boundary
				matrix[i][i] = 1; 
			} else if (i > total_points - (2 * discretization_points) - 1) { // bottom boundary
				matrix[i][i] = 1; 
			} else if (i % discretization_points < 2) { // left boundary
				matrix[i][i] = 1;
			} else if (i % discretization_points > discretization_points - 3) { // right boundary
				matrix[i][i] = 1; 
			} else { // interior points
				matrix[i][i - (2 * discretization_points)] = -1;
				matrix[i][i - discretization_points] = 16;
				matrix[i][i - 2] = -1.0;
				matrix[i][i - 1] =  16.0;
				matrix[i][i]     = -60.0;
				matrix[i][i + 1] =  16.0;
				matrix[i][i + 2] = -1.0;
				matrix[i][i + discretization_points] = 16;
				matrix[i][i + (2 * discretization_points)] = -1;
			}
		}
	}
	return matrix;
}

double* initialize_T_vector_2D(int discretization_points, int order, double X_MIN, double X_MAX, double Y_MIN, double Y_MAX){
	int total_points = discretization_points * discretization_points;
	double initValue = 0.0;
	double* T_vector = new double[total_points];
	if ( !T_vector ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for temperature vector \n");
		exit(1);
	}
	double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);
	double delta_y = (Y_MAX - Y_MIN) / (discretization_points - 1);
	double* pos = new double[2]; // x,y
	int i = 0;
	int j = 0;
	fill_n(T_vector, total_points, initValue); // Initialize all values to 0.0
	if (order == 2) {
		for (int h = 0; h < total_points; h++) {
			if (h > 0) { // initial point (0,0) treated differently
				i++; // x-position
				if (h % discretization_points == 0) { // beginning of a row in 2D space
					j++; // y-position
					i = 0; // reset x = 0 after each row of 2D space
				}
			}
			if (j == 0) { // bottom row of 2D space, y = 0
				if (i == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (i == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else {
					pos[0] = i * delta_x;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				}
			} else if (j == discretization_points - 1) { // top row of 2D space, y = W
				if (i == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else if (i == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[0] = i * delta_x;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				}
			} else if (i == 0) { // left column of 2D space, x = 0
				if (j == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (j == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				}
			} else if (i == discretization_points - 1) { // right column of 2D space, x = L
				if (j == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else if (j == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				}
			}
		}
	} else if (order == 4) {
		for (int h = 0; h < total_points; h++) {
			if (h > 0) { // initial point (0,0) treated differently
				i++; // x-position
				if (h % discretization_points == 0) { // beginning of a row in 2D space
					j++; // y-position
					i = 0; // reset x = 0 after each row of 2D space
				}
			}
			if (j == 0) { // bottom row of 2D space, y = 0
				if (i == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (i == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else {
					pos[0] = i * delta_x;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				}
			} else if (j == 1) {// bottom, "internal boundary", y = delta_y
				pos[1] = j * delta_y;
				if (i == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				} else if (i == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				} else {
					pos[0] = i * delta_x;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}
			} else if (j == discretization_points - 1) { // top row of 2D space, y = W
				if (i == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else if (i == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[0] = i * delta_x;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				}
			} else if (j == discretization_points - 2) { // top, "internal boundary", y = W - delta_y
				pos[1] = j * delta_y;
				if (i == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				} else if (i == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				} else {
					pos[0] = i * delta_x;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}	
			} else if (i == 0) { // left column of 2D space, x = 0
				if (j == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (j == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					T_vector[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				}
			} else if (i == 1) { // left, "internal boundary", x = delta_x
				pos[0] = i * delta_x;
				if (j == 0) {
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				} else if (j == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				} else {
					pos[1] = j * delta_y;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}
			} else if (i == discretization_points - 1) { // right column of 2D space, x = L
				if (j == 0) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else if (j == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					T_vector[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				}
			} else if (i == discretization_points - 2) { // right, "internal boundary", x = L - delta_x
				pos[0] = i * delta_x;
				if (j == 0) {
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				} else if (j == discretization_points - 1) {
					T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				} else {
					pos[1] = j * delta_y;
					T_vector[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}
			} 
		}
	}
	delete[] pos;
	return T_vector;
}

double* initialize_RHS_2D(int discretization_points, int order, double X_MIN, double X_MAX, double Y_MIN, double Y_MAX, double k_0) {
	int total_points = discretization_points * discretization_points;
	double* RHS = new double[total_points];
	if ( !RHS ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for right-hand side vector \n");
		exit(1);
	}
	double q;
	int i = 0;
	int j = 0;
	double* pos = new double[2]; // x,y
	double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);
	double delta_y = (Y_MAX - Y_MIN) / (discretization_points - 1);
	
	if (order == 2) {
		for (int h = 0; h < total_points; h++) {
			if (h > 0) { // initial point (0,0) treated differently
				i++; // x-position
				if (h % discretization_points == 0) { // beginning of a row in 2D space
					j++; // y-position
					i = 0; // reset x = 0 after each row of 2D space
				}
			}
			if (j == 0) { // bottom row of 2D space, y = 0
				if (i == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (i == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else {
					pos[0] = i * delta_x;
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				}
			} else if (j == discretization_points - 1) { // top row of 2D space, y = W
				if (i == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else if (i == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[0] = i * delta_x;
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				}
			} else if (i == 0) { // left column of 2D space, x = 0
				if (j == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (j == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					RHS[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				}
			} else if (i == discretization_points - 1) { // right column of 2D space, x = L
				if (j == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else if (j == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					RHS[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				}
			} else {
				pos[0] = i * delta_x;
				pos[1] = j * delta_y;
				q = masa_eval_2d_source_t(pos[0], pos[1]);
				RHS[h] = -(q * delta_x * delta_y) / k_0;
			}
		}
	} else if (order == 4) {
		for (int h = 0; h < total_points; h++) {
			if (h > 0) { // initial point (0,0) treated differently
				i++; // x-position
				if (h % discretization_points == 0) { // beginning of a row in 2D space
					j++; // y-position
					i = 0; // reset x = 0 after each row of 2D space
				}
			}
			if (j == 0) { // bottom row of 2D space, y = 0
				if (i == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (i == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else {
					pos[0] = i * delta_x;
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				}
			} else if (j == 1) {// bottom, "internal boundary", y = delta_y
				pos[1] = j * delta_y;
				if (i == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				} else if (i == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				} else {
					pos[0] = i * delta_x;
					RHS[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}
			} else if (j == discretization_points - 1) { // top row of 2D space, y = W
				if (i == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else if (i == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[0] = i * delta_x;
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				}
			} else if (j == discretization_points - 2) { // top, "internal boundary", y = W - delta_y
				pos[1] = j * delta_y;
				if (i == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				} else if (i == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				} else {
					pos[0] = i * delta_x;
					RHS[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}	
			} else if (i == 0) { // left column of 2D space, x = 0
				if (j == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MIN);
				} else if (j == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MIN, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					RHS[h] = masa_eval_2d_exact_t(X_MIN, pos[1]);
				}
			} else if (i == 1) { // left, "internal boundary", x = delta_x
				pos[0] = i * delta_x;
				if (j == 0) {
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				} else if (j == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				} else {
					pos[1] = j * delta_y;
					RHS[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}
			} else if (i == discretization_points - 1) { // right column of 2D space, x = L
				if (j == 0) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MIN);
				} else if (j == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
				} else {
					pos[1] = j * delta_y;
					RHS[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
				}
			} else if (i == discretization_points - 2) { // right, "internal boundary", x = L - delta_x
				pos[0] = i * delta_x;
				if (j == 0) {
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MIN);
				} else if (j == discretization_points - 1) {
					RHS[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
				} else {
					pos[1] = j * delta_y;
					RHS[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
				}
			} else {
				pos[0] = i * delta_x;
				pos[1] = j * delta_y;
				q = masa_eval_2d_source_t(pos[0], pos[1]);
				RHS[h] = (-q * delta_x * delta_y * 12) / k_0;
			}
		}
	}
	delete[] pos;
	return RHS;
}