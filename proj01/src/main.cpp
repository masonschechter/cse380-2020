#include <string.h>
#include <iomanip>
#include<sys/time.h>
#include<time.h>
#include<unistd.h>
#include <masa.h>


#include <input_parsing.cpp>
#include <iterative_solvers.cpp>
#include <vectors_and_matrices_1D.cpp>
#include <vectors_and_matrices_2D.cpp>

using namespace std;
using namespace MASA;
using namespace GRVY;

double* run_1D_pipepine(settings inputs, GRVY_Timer_Class &timer) {
	double* solution;

	// Initialize analytical solution
	timer.BeginTimer("Initializing analytical solution");
	masa_init("1D_case", "heateq_1d_steady_const");
	masa_set_param("A_x", inputs.A_X);
	masa_set_param("k_0", inputs.K_0);
	timer.EndTimer("Initializing analytical solution");
	

	// Initialize linear system
	timer.BeginTimer("Initializing linear system");
	double** mat = make_1D_coefficient_matrix(inputs.DISCRETIZATION_POINTS, inputs.ORDER);
	double* t0 = initialize_T_vector_1D(inputs.DISCRETIZATION_POINTS, inputs.ORDER, inputs.X_MIN, inputs.X_MAX);
	double* RHS = initialize_RHS_1D(inputs.DISCRETIZATION_POINTS, inputs.ORDER, inputs.X_MIN, inputs.X_MAX, inputs.K_0);
	if (inputs.OUTPUT_MODE == "DEBUG") {
		// Print the linear system
		grvy_printf(GRVY_DEBUG, "Coefficient Matrix: \n");
		for (int i = 0; i < inputs.DISCRETIZATION_POINTS; i++) {
			grvy_printf(GRVY_DEBUG, "Row %d: \n", i+1);
			for (int j = 0; j < inputs.DISCRETIZATION_POINTS; j++){
				grvy_printf(GRVY_DEBUG, "%.2f ", mat[i][j]);
				if (j == inputs.DISCRETIZATION_POINTS - 1) {
					grvy_printf(GRVY_DEBUG, "\n");
				}
			}
		}
		grvy_printf(GRVY_DEBUG, "\n");
		grvy_printf(GRVY_DEBUG, "Initial Temperature Vector: \n");
		for (int i = 0; i < inputs.DISCRETIZATION_POINTS; i++) {
			grvy_printf(GRVY_DEBUG, "%.15f ", t0[i]);
		}
		grvy_printf(GRVY_DEBUG, "\n\n");
		grvy_printf(GRVY_DEBUG, "Initial Right-Hand Side Vector: \n");
		for (int i = 0; i < inputs.DISCRETIZATION_POINTS; i++) {
			grvy_printf(GRVY_DEBUG, "%.15f ", RHS[i]);
		}
		grvy_printf(GRVY_DEBUG, "\n\n");
	}
	timer.EndTimer("Initializing linear system");

	// Solve the linear system
	timer.BeginTimer("Iterative solution");
	if (inputs.SOLVER_TYPE == "Jacobi") {
		solution = Jacobi(inputs.DISCRETIZATION_POINTS, mat, t0, RHS, inputs.ERROR_THRESHOLD, inputs.MAX_ITERATIONS, inputs.OUTPUT_MODE);
	} else if (inputs.SOLVER_TYPE == "Gauss-Seidel") {
		solution = Gauss_Seidel(inputs.DISCRETIZATION_POINTS, mat, t0, RHS, inputs.ERROR_THRESHOLD, inputs.MAX_ITERATIONS, inputs.OUTPUT_MODE);
	}
	timer.EndTimer("Iterative solution");
	delete[] t0;
	delete[] RHS;
	for (int i = 0; i < inputs.DISCRETIZATION_POINTS; i++) {
		delete[] mat[i];
	}
	delete[] mat;
	return solution;
}

double* run_2D_pipepine(settings inputs, GRVY_Timer_Class &timer) {
	double* solution;

	// Initialize analytical solution
	timer.BeginTimer("Initializing analytical solution");
	masa_init("2D_case", "heateq_2d_steady_const");
	masa_set_param("A_x", inputs.A_X);
	masa_set_param("B_y", inputs.B_Y);
	masa_set_param("k_0", inputs.K_0);
	timer.EndTimer("Initializing analytical solution");

	// Initialize linear system
	timer.BeginTimer("Initializing linear system");
	double** mat = make_2D_coefficient_matrix(inputs.DISCRETIZATION_POINTS, inputs.ORDER);
	double* t0 = initialize_T_vector_2D(inputs.DISCRETIZATION_POINTS, inputs.ORDER, inputs.X_MIN, inputs.X_MAX, inputs.Y_MIN, inputs.Y_MAX);
	double* RHS = initialize_RHS_2D(inputs.DISCRETIZATION_POINTS, inputs.ORDER, inputs.X_MIN, inputs.X_MAX, inputs.Y_MIN, inputs.Y_MAX, inputs.K_0);
	if (inputs.OUTPUT_MODE == "DEBUG") {
		// Print the linear system
		grvy_printf(GRVY_DEBUG, "Coefficient Matrix: \n");
		for (int i = 0; i < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; i++) {
			grvy_printf(GRVY_DEBUG, "Row %d: \n", i+1);
			for (int j = 0; j < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; j++){
				grvy_printf(GRVY_DEBUG, "%.2f ", mat[i][j]);
				if (j == inputs.DISCRETIZATION_POINTS - 1) {
					grvy_printf(GRVY_DEBUG, "\n");
				}
			}
		}
		grvy_printf(GRVY_DEBUG, "\n");
		grvy_printf(GRVY_DEBUG, "Initial Temperature Vector: \n");
		for (int i = 0; i < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; i++) {
			grvy_printf(GRVY_DEBUG, "%.15f ", t0[i]);
		}
		grvy_printf(GRVY_DEBUG, "\n\n");
		grvy_printf(GRVY_DEBUG, "Initial Right-Hand Side Vector: \n");
		for (int i = 0; i < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; i++) {
			grvy_printf(GRVY_DEBUG, "%.15f ", RHS[i]);
		}
		grvy_printf(GRVY_DEBUG, "\n\n");
	}
	timer.EndTimer("Initializing linear system");
	// Solve the linear system
	timer.BeginTimer("Iterative solution");
	int total_points = inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS;
	if (inputs.SOLVER_TYPE == "Jacobi") {
		solution = Jacobi(total_points, mat, t0, RHS, inputs.ERROR_THRESHOLD, inputs.MAX_ITERATIONS, inputs.OUTPUT_MODE);
	} else if (inputs.SOLVER_TYPE == "Gauss-Seidel") {
		solution = Gauss_Seidel(total_points, mat, t0, RHS, inputs.ERROR_THRESHOLD, inputs.MAX_ITERATIONS, inputs.OUTPUT_MODE);
	}
	timer.EndTimer("Iterative solution");
	delete[] t0;
	delete[] RHS;
	for (int i = 0; i < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; i++) {
		delete[] mat[i];
	}
	delete[] mat;
	return solution;
}

int main() {
	GRVY_Timer_Class timer;
	timer.Init("My Timer");
	cout << setprecision(15);

	settings inputs = parse_inputs();

	if (inputs.OUTPUT_MODE == "INFO") {
		grvy_log_setlevel(GRVY_INFO);
	} else if (inputs.OUTPUT_MODE == "DEBUG") {
		grvy_log_setlevel(GRVY_DEBUG);
	}

	double* numerical_solution;
	double* analytical_solution;
	double error;

	if (inputs.DIMENSIONS == 1) {
		numerical_solution = run_1D_pipepine(inputs, timer);
		analytical_solution = analytical_solution_1D(inputs.DISCRETIZATION_POINTS, inputs.X_MIN, inputs.X_MAX);
		if (inputs.OUTPUT_MODE == "DEBUG") {
			// Print the numerical and analytical solution vectors
			grvy_printf(GRVY_DEBUG, "Numerical Solution Temperature Vector: \n");
			for (int i = 0; i < inputs.DISCRETIZATION_POINTS; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", numerical_solution[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
			grvy_printf(GRVY_DEBUG, "Analytical Solution Temperature Vector: \n");
			for (int i = 0; i < inputs.DISCRETIZATION_POINTS; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", analytical_solution[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
		}
		for (int i = 0; i < inputs.DISCRETIZATION_POINTS; i++) {
			// Print solution vector
			grvy_printf(GRVY_INFO, "%.15f ", numerical_solution[i]);
		}
		grvy_printf(GRVY_INFO, "\n");
		error = L2_norm(inputs.DISCRETIZATION_POINTS, numerical_solution, analytical_solution);
	} else if (inputs.DIMENSIONS == 2) {
		numerical_solution = run_2D_pipepine(inputs, timer);
		analytical_solution = analytical_solution_2D(inputs.DISCRETIZATION_POINTS, inputs.X_MIN, inputs.X_MAX, inputs.Y_MIN, inputs.Y_MAX);
		if (inputs.OUTPUT_MODE == "DEBUG") {
			// Print the numerical and analytical solution vectors
			grvy_printf(GRVY_DEBUG, "Numerical Solution Temperature Vector: \n");
			for (int i = 0; i < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", numerical_solution[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
			grvy_printf(GRVY_DEBUG, "Analytical Solution Temperature Vector: \n");
			for (int i = 0; i < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; i++) {
				grvy_printf(GRVY_DEBUG, "%.15f ", analytical_solution[i]);
			}
			grvy_printf(GRVY_DEBUG, "\n\n");
		} else {
			for (int i = 0; i < inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS; i++) {
				// Print solution vector
				grvy_printf(GRVY_INFO, "%.15f ", numerical_solution[i]);
			}
		}
		grvy_printf(GRVY_INFO, "\n");
		error = L2_norm(inputs.DISCRETIZATION_POINTS * inputs.DISCRETIZATION_POINTS, numerical_solution, analytical_solution);
	}
	if (inputs.RUN_MODE == "VERIFICATION"){
		grvy_printf(GRVY_INFO, "L2 Norm of numerical vs analytical solutions: %.15f", error);
	}
	timer.Finalize();
	timer.Summarize();
	delete[] numerical_solution;
	delete[] analytical_solution;
	return 0;
}
