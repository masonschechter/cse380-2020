#include <masa.h>
#include <grvy.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>

#include <iterative_solvers.hpp>
#include <vectors_and_matrices_2D.hpp>

using namespace std;
using namespace MASA;
using namespace GRVY;

int main() {
    double X_MIN = 0.0;
    double X_MAX = 1.0;
    double Y_MIN = 0.0;
    double Y_MAX = 1.0;
    int DISCRETIZATION_POINTS = 10;
    int ORDER = 2;
    double** matrix;
    double* t0;
    double* RHS;
    double* solution;
    double* analytical_sol;

    masa_init("2D_tests", "heateq_2d_steady_const");
	masa_set_param("A_x", 10.0);
    masa_set_param("B_y", 10.0);
	masa_set_param("k_0", 1.0);

    matrix = make_2D_coefficient_matrix(DISCRETIZATION_POINTS, ORDER);
    t0 = initialize_T_vector_2D(DISCRETIZATION_POINTS, ORDER, X_MIN, X_MAX, Y_MIN, Y_MAX);
    RHS = initialize_RHS_2D(DISCRETIZATION_POINTS, ORDER, X_MIN, X_MAX, Y_MIN, Y_MAX, 1.0);

    solution = Jacobi(DISCRETIZATION_POINTS * DISCRETIZATION_POINTS, matrix, t0, RHS, 1e-10, 250000, "INFO");
    analytical_sol = analytical_solution_2D(DISCRETIZATION_POINTS, X_MIN, X_MAX, X_MIN, X_MAX);


    output_to_file(DISCRETIZATION_POINTS, 2, "VERIFICATION", X_MIN, X_MAX, solution, analytical_sol, "jacobi_2D_2nd_order_test_out.h5");

    return 0;
}