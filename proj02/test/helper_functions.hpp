#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <hdf5.h>

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
	int i = 0;
	int j = 0;
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
	delete[] pos;
	return T_vector;
}

void output_to_file(int discretization_points, int dimensions, string RUN_MODE, double X_MIN, double X_MAX, double* numerical_solution, double* analytical_solution, char* filename) {
	hid_t file, dataset, dataspace;
	herr_t status;
	file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	double delta_x = (double)(X_MAX - X_MIN) / (discretization_points - 1);

	if (dimensions == 1) {
		hsize_t coord_dims[2], temp_dims[1];
		coord_dims[0] = discretization_points;
		coord_dims[1] = 2; // x,y
		temp_dims[0] = discretization_points;
		double coordinates[discretization_points][2];
		// Populate 2D coordinate array, N x (x,y)
		for (int i = 0; i < discretization_points; i++) {
			coordinates[i][0] = i * delta_x;
			coordinates[i][1] = 0.0;
		}
		// coordinate grid
		dataspace = H5Screate_simple(2, coord_dims, NULL);
		dataset = H5Dcreate2(file, "coordinates", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coordinates);

		// numerical solution grid
		dataspace = H5Screate_simple(1, temp_dims, NULL);
		dataset = H5Dcreate2(file, "numerical solution", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, numerical_solution);

		if (RUN_MODE == "VERIFICATION") {
			// analytical solution grid
			dataspace = H5Screate_simple(1, temp_dims, NULL);
			dataset = H5Dcreate2(file, "analytical solution", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, analytical_solution);
		}
		// Clean up
		H5Sclose(dataspace);
		H5Dclose(dataset);
		H5Fclose(file);

	} else if (dimensions == 2) {
		hsize_t coord_dims[3], temp_dims[2];
		coord_dims[0] = discretization_points;
		coord_dims[1] = discretization_points;
		coord_dims[2] = 2; // x,y
		temp_dims[0] = discretization_points;
		temp_dims[1] = discretization_points;
		double coordinates[discretization_points][discretization_points][2];
		double numerical_temps[discretization_points][discretization_points];
		double analytical_temps[discretization_points][discretization_points];
		// Populate 3D coordinate array, N x N x (x,y)
		for (int i = 0; i < discretization_points; i++) {
			for (int j = 0; j < discretization_points; j++) {
				coordinates[i][j][0] = j * delta_x;
				coordinates[i][j][1] = i * delta_x;
			}
		}
		// Map 1D solution vectors to 2D space to align with coordinate grid ^
		int h = 0; 
		for (int i = 0; i < discretization_points; i++) {
			for (int j = 0; j < discretization_points; j++) {
				numerical_temps[i][j] = numerical_solution[h];
				if (RUN_MODE == "VERIFICATION") {
					analytical_temps[i][j] = analytical_solution[h];
				}
				h++;
			}
		}
		// coordinate grid
		dataspace = H5Screate_simple(3, coord_dims, NULL);
		dataset = H5Dcreate2(file, "coordinates", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coordinates);

		// numerical solution grid
		dataspace = H5Screate_simple(2, temp_dims, NULL);
		dataset = H5Dcreate2(file, "numerical solution", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, numerical_temps);

		if (RUN_MODE == "VERIFICATION") {
			// analytical solution grid
			dataspace = H5Screate_simple(2, temp_dims, NULL);
			dataset = H5Dcreate2(file, "analytical solution", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, analytical_temps);
		}
		// Clean up
		H5Sclose(dataspace);
		H5Dclose(dataset);
		H5Fclose(file);
	}
}