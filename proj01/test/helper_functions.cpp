#include <iostream>
#include <stdlib.h>
#include <cmath>

using namespace std;
using namespace MASA;

double L2_norm(int discretization_points, double* T_old, double* T_new) {
	double distance = 0.0;
	discretization_points = (double)discretization_points;
	for (int i = 0; i < discretization_points; i++) {
		distance += ((T_old[i] - T_new[i]) * (T_old[i] - T_new[i]));
	}
	distance = distance / discretization_points;
	distance = sqrt(distance);
	return distance;
}

double* analytical_solution_1D(int discretization_points, double X_MIN, double X_MAX) {
	double* T_vector = new double[discretization_points];
	if ( !T_vector ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for analytical temperature vector \n");
		exit(1);
	}
	double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);
	double pos;
	T_vector[discretization_points - 1] = masa_eval_1d_exact_t(X_MAX); // x = L
	for (int i = 0; i < discretization_points - 1; i++) {
		pos = i * delta_x;
		T_vector[i] = masa_eval_1d_exact_t(pos);
	}
	return T_vector;
}

double* analytical_solution_2D(int discretization_points, double X_MIN, double X_MAX, double Y_MIN, double Y_MAX) {
	int total_points = discretization_points * discretization_points;
	double* T_vector = new double[total_points];
	if ( !T_vector ) { // check for NULL pointer
		grvy_printf(GRVY_ERROR, "ERROR: Could not allocate memory for analytical temperature vector \n");
		exit(1);
	}
	double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);
	double delta_y = (Y_MAX - Y_MIN) / (discretization_points - 1);
	double* pos = new double[2]; // x,y
	int i;
	int j;
	for (int h = 0; h < total_points; h++) {
		if (h > 0) { // initial point (0,0) treated differently
			i++; // x-position
			if (h % discretization_points == 0) { // beginning of a row in 2D space
				j++; // y-position
				i = 0; // reset x = 0 after each row of 2D space
			}
		}
		if (i == discretization_points - 1 && j == discretization_points - 1) {
			T_vector[h] = masa_eval_2d_exact_t(X_MAX, Y_MAX);
		} else {
			pos[0] = i * delta_x;
			pos[1] = j * delta_y;
			if (i == discretization_points - 1) {
				T_vector[h] = masa_eval_2d_exact_t(X_MAX, pos[1]);
			} else if (j == discretization_points - 1) {
				T_vector[h] = masa_eval_2d_exact_t(pos[0], Y_MAX);
			} else {
				T_vector[h] = masa_eval_2d_exact_t(pos[0], pos[1]);
			}
		}
	}
	return T_vector;
}