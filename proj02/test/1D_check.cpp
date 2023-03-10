#include <masa.h>
#include <grvy.h>

#include <vectors_and_matrices_1D.hpp>

using namespace std;
using namespace MASA;

int main() {
    double X_MIN = 0.0;
    double X_MAX = 1.0;
    int DISCRETIZATION_POINTS = 10;
    double** matrix;
    double* t0;
    double* RHS;

    masa_init("1D_tests", "heateq_1d_steady_const");
	masa_set_param("A_x", 10.0);
	masa_set_param("k_0", 1.0);

    for (int order = 2; order < 5; order += 2) {
        matrix = make_1D_coefficient_matrix(DISCRETIZATION_POINTS, order);
        for (int i = 0; i < (DISCRETIZATION_POINTS); i++) {
	    	cout << endl;
	    	for (int j = 0; j < DISCRETIZATION_POINTS; j++){
	    		cout << matrix[i][j] << ' ';
	    	}
        }
        cout << '\n' << endl;
        t0 = initialize_T_vector_1D(DISCRETIZATION_POINTS, order, X_MIN, X_MAX);
        for (int i = 0; i < DISCRETIZATION_POINTS; i++) {
            cout << t0[i] << " ";
        }
        cout << '\n' << endl;
        RHS = initialize_RHS_1D(DISCRETIZATION_POINTS, order, X_MIN, X_MAX, 1.0);
        for (int i = 0; i < DISCRETIZATION_POINTS; i++) {
            cout << RHS[i] << " ";
        }
        cout << '\n' << endl;
    }
    return 0;
}
