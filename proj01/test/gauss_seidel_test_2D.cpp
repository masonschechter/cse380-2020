#include <masa.h>
#include <grvy.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>

#include <iterative_solvers.cpp>
#include <vectors_and_matrices_2D.cpp>

using namespace std;
using namespace MASA;
using namespace GRVY;

int main() {
    double X_MIN = 0.0;
    double X_MAX = 1.0;
    double Y_MIN = 0.0;
    double Y_MAX = 1.0;
    int DISCRETIZATION_POINTS = 10;
    double** matrix;
    double* t0;
    double* RHS;
    double* solution;

    masa_init("2D_tests", "heateq_2d_steady_const");
	masa_set_param("A_x", 10.0);
    masa_set_param("B_y", 10.0);
	masa_set_param("k_0", 1.0);

    for (int order = 2; order < 5; order += 2) {
        matrix = make_2D_coefficient_matrix(DISCRETIZATION_POINTS, order);
        t0 = initialize_T_vector_2D(DISCRETIZATION_POINTS, order, X_MIN, X_MAX, Y_MIN, Y_MAX);
        RHS = initialize_RHS_2D(DISCRETIZATION_POINTS, order, X_MIN, X_MAX, Y_MIN, Y_MAX, 1.0);

        solution = Gauss_Seidel(DISCRETIZATION_POINTS * DISCRETIZATION_POINTS, matrix, t0, RHS, 1e-10, 250000, "INFO");
    }

    return 0;
}