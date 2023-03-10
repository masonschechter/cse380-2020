#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <petsc.h>

using namespace std;
using namespace MASA;


PetscErrorCode initialize_coefficient_matrix(Mat* mat, PetscErrorCode ierr, int discretization_points, int order, int dimensions) {
    
    PetscInt j;
    if (dimensions == 1) {
        PetscInt disc_pts = discretization_points;
        MatCreate(PETSC_COMM_WORLD, mat);
        MatSetSizes(*mat, PETSC_DECIDE, PETSC_DECIDE, disc_pts, disc_pts);
        MatSetUp(*mat);

        if (order == 2) {
            PetscInt cols[3];
            PetscReal values[3];
            for (int i = 0; i < discretization_points; i++) {
                j = i;
                if (i == 0 || i == discretization_points - 1) {
                    cols[0] = i; values[0] = 1.0;
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else {
                    cols[0] = i - 1; cols[1] = i; cols[2] = i + 1;
                    values[0] = 1.0; values[1] = -2.0; values[2] = 1.0;
                    ierr = MatSetValues(*mat, 1, &j, 3, cols, values, INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        } else if (order == 4) {
            PetscInt cols[5];
            PetscReal values[5];
            for (int i = 0; i < discretization_points; i++) {
                j = i;
                if (i == 0 || i == 1 || i == discretization_points - 2 || i == discretization_points - 1) {
                    cols[0] = i; values[0] = 1.0;
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else {
                    cols[0] = i - 2; cols[1] = i - 1; cols[2] = i; cols[3] = i + 1; cols[4] = i +2;
                    values[0] = -1.0; values[1] = 16.0; values[2] = -30.0; values[3] = 16.0; values[4] = -1.0;
                    ierr = MatSetValues(*mat, 1, &j, 5, cols, values, INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }
    } else if (dimensions == 2) {
        PetscInt total_points = discretization_points * discretization_points;
        MatCreate(PETSC_COMM_WORLD, mat);
        MatSetSizes(*mat, PETSC_DECIDE, PETSC_DECIDE, total_points, total_points);
        MatSetUp(*mat);

        if (order == 2) {
            PetscInt cols[5];
            PetscReal values[5];
            for (int i = 0; i < total_points; i++) {
                j = i;
                cols[0] = i; values[0] = 1.0;
                if (i < discretization_points) { // top boundary
				    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else if (i > total_points - discretization_points - 1) { // bottom boundary
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else if (i % discretization_points == 0) { // left boundary
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else if (i % discretization_points == discretization_points - 1) { // right boundary
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else { // interior points
                    cols[0] = i - discretization_points; cols[1] = i - 1; cols[2] = i; cols[3] = i + 1; cols[4] = i + discretization_points;
                    values[0] = 1.0; values[1] = 1.0; values[2] = -4.0; values[3] = 1.0; values[4] = 1.0;
                    ierr = MatSetValues(*mat, 1, &j, 5, cols, values, INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        } if (order == 4) {
            PetscInt cols[9];
            PetscReal values[9];
            for (int i = 0; i < discretization_points*discretization_points; i++) {
                j = i;
                cols[0] = i; values[0] = 1.0;
                if (i < 2 * discretization_points) { // top boundary
				    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr); 
                } else if (i > total_points - (2 * discretization_points) - 1) { // bottom boundary
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else if (i % discretization_points < 2) { // left boundary
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else if (i % discretization_points > discretization_points - 3) { // right boundary
                    ierr = MatSetValues(*mat, 1, &j, 1, &cols[0], &values[0], INSERT_VALUES); CHKERRQ(ierr);
                } else { // interior points
                    cols[0] = i - (2 * discretization_points); cols[1] = i - discretization_points; cols[2] = i - 2; cols[3] = i - 1; cols[4] = i;
                    cols[5] = i + 1; cols[6] = i + 2; cols[7] = i + discretization_points; cols[8] = i + (2 * discretization_points);
                    values[0] = -1.0; values[1] = 16; values[2] = -1.0; values[3] = 16.0; values[4] = -60.0;
                    values[5] = 16.0; values[6] = -1.0; values[7] = 16.0; values[8] = -1.0;
                    ierr = MatSetValues(*mat, 1, &j, 9, cols, values, INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }
    }
    ierr = MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
}

PetscErrorCode initialize_RHS_and_Sol_vectors(Vec* RHS, Vec* Sol, PetscErrorCode ierr, int discretization_points, int order, int dimensions, double X_MIN, double X_MAX, double Y_MIN, double Y_MAX, double k_0) {
    cout << setprecision(15);
    if (dimensions == 1) {
        PetscInt disc_pts = discretization_points;
        VecCreate(PETSC_COMM_WORLD, RHS);
        VecCreate(PETSC_COMM_WORLD, Sol);
        VecSetSizes(*RHS, PETSC_DECIDE, disc_pts);
        VecSetSizes(*Sol, PETSC_DECIDE, disc_pts);
        VecSetUp(*RHS);
        VecSetUp(*Sol);
        double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);
        double pos;
        double q;
        PetscReal temp;
        PetscInt j;

        if (order == 2) {
            for (int i = 0; i < discretization_points; i++) {
                j = i;
                if (i == 0) {
                    temp = masa_eval_1d_exact_t(X_MIN);
                    ierr = VecSetValues(*RHS, 1, &j, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    ierr = VecSetValues(*Sol, 1, &j, &temp, INSERT_VALUES); CHKERRQ(ierr);
                } else if (i == discretization_points - 1) {
                    temp = masa_eval_1d_exact_t(X_MAX);
                    ierr = VecSetValues(*RHS, 1, &j, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    ierr = VecSetValues(*Sol, 1, &j, &temp, INSERT_VALUES); CHKERRQ(ierr);
                } else {
                    pos = i * delta_x;
                    q = masa_eval_1d_source_t(pos);
                    temp = -(q * delta_x * delta_x) / k_0;
                    PetscScalar zero = 0.0;
                    ierr = VecSetValues(*RHS, 1, &j, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    ierr = VecSetValues(*Sol, 1, &j, &zero, INSERT_VALUES); CHKERRQ(ierr);
                }

            }
        }
    } else {
        PetscInt disc_pts = discretization_points * discretization_points;
        VecCreate(PETSC_COMM_WORLD, RHS);
        VecCreate(PETSC_COMM_WORLD, Sol);
        VecSetSizes(*RHS, PETSC_DECIDE, disc_pts);
        VecSetSizes(*Sol, PETSC_DECIDE, disc_pts);
        VecSetUp(*RHS);
        VecSetUp(*Sol);
        double delta_x = (X_MAX - X_MIN) / (discretization_points - 1);
        double delta_y = (Y_MAX - Y_MIN) / (discretization_points - 1);
        double pos[2];
        double q;
        int i = 0;
        int j = 0;
        PetscReal temp;
        PetscInt k;
        if (order == 2) {
            for (int h = 0; h < disc_pts; h++) {
                k = h;
                if (h > 0) { // initial point (0,0) treated differently
                    i++; // x-position
                    if (h % discretization_points == 0) { // beginning of a row in 2D space
                        j++; // y-position
                        i = 0; // reset x = 0 after each row of 2D space
                    }
                }
                if (j == 0) { // bottom row of 2D space, y = 0
                    if (i == 0) {
                        temp = masa_eval_2d_exact_t(X_MIN, Y_MIN);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else if (i == discretization_points - 1) {
                        temp = masa_eval_2d_exact_t(X_MAX, Y_MIN);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else {
                        pos[0] = i * delta_x;
                        temp = masa_eval_2d_exact_t(pos[0], Y_MIN);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    }
                } else if (j == discretization_points - 1) { // top row of 2D space, y = W
                    if (i == 0) {
                        temp = masa_eval_2d_exact_t(X_MIN, Y_MAX);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else if (i == discretization_points - 1) {
                        temp = masa_eval_2d_exact_t(X_MAX, Y_MAX);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else {
                        pos[0] = i * delta_x;
                        temp = masa_eval_2d_exact_t(pos[0], Y_MAX);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    }
                } else if (i == 0) { // left column of 2D space, x = 0
                    if (j == 0) {
                        temp = masa_eval_2d_exact_t(X_MIN, Y_MIN);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else if (j == discretization_points - 1) {
                        temp = masa_eval_2d_exact_t(X_MIN, Y_MAX);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else {
                        pos[1] = j * delta_y;
                        temp = masa_eval_2d_exact_t(X_MIN, pos[1]);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    }
                } else if (i == discretization_points - 1) { // right column of 2D space, x = L
                    if (j == 0) {
                        temp = masa_eval_2d_exact_t(X_MAX, Y_MIN);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else if (j == discretization_points - 1) {
                        temp = masa_eval_2d_exact_t(X_MAX, Y_MAX);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    } else {
                        pos[1] = j * delta_y;
                        temp = masa_eval_2d_exact_t(X_MAX, pos[1]);
                        ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = VecSetValues(*Sol, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    }
                } else {
                    pos[0] = i * delta_x;
				    pos[1] = j * delta_y;
                    q = masa_eval_2d_source_t(pos[0], pos[1]);
                    temp = -(q * delta_x * delta_y) / k_0;
                    PetscScalar zero = 0.0;
                    ierr = VecSetValues(*RHS, 1, &k, &temp, INSERT_VALUES); CHKERRQ(ierr);
                    ierr = VecSetValues(*Sol, 1, &k, &zero, INSERT_VALUES); CHKERRQ(ierr);
                } 
            }
        }
    }
    ierr = VecAssemblyBegin(*RHS); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*RHS); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(*Sol); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*Sol); CHKERRQ(ierr);
}


PetscErrorCode petsc_pipeline(settings inputs, GRVY_Timer_Class &timer) {

    PetscErrorCode ierr;
    KSP Solver;
    PC Prec;
    Mat mat;
    Vec RHS, Sol;
    PetscReal ERROR_THRESHOLD = inputs.ERROR_THRESHOLD;
    PetscReal MAX_ITERATIONS = inputs.MAX_ITERATIONS;

    ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
    if (inputs.DIMENSIONS == 1) {
        // Initialize analytical solution
	    timer.BeginTimer("Initializing analytical solution");
	    masa_init("1D_case", "heateq_1d_steady_const");
	    masa_set_param("A_x", inputs.A_X);
	    masa_set_param("k_0", inputs.K_0);
	    timer.EndTimer("Initializing analytical solution");
    } else if (inputs.DIMENSIONS == 2) {
        // Initialize analytical solution
        timer.BeginTimer("Initializing analytical solution");
        masa_init("2D_case", "heateq_2d_steady_const");
        masa_set_param("A_x", inputs.A_X);
        masa_set_param("B_y", inputs.B_Y);
        masa_set_param("k_0", inputs.K_0);
        timer.EndTimer("Initializing analytical solution");
    }

    // Initialize linear system
	timer.BeginTimer("Initializing linear system");
    ierr = initialize_coefficient_matrix(&mat, ierr, inputs.DISCRETIZATION_POINTS, inputs.ORDER, inputs.DIMENSIONS);
    CHKERRQ(ierr);
    ierr = initialize_RHS_and_Sol_vectors(&RHS, &Sol, ierr, inputs.DISCRETIZATION_POINTS, inputs.ORDER, inputs.DIMENSIONS, inputs.X_MIN, inputs.X_MAX, inputs.Y_MIN, inputs.Y_MAX, inputs.K_0);
    CHKERRQ(ierr);
    timer.EndTimer("Initializing linear system");

    timer.BeginTimer("Iterative solution");
    ierr = KSPCreate(PETSC_COMM_WORLD, &Solver); CHKERRQ(ierr);
    ierr = KSPSetOperators(Solver, mat, mat); CHKERRQ(ierr);
    // ierr = KSPSetType(Solver, KSPGMRES); CHKERRQ(ierr);
    // ierr = KSPGetPC(Solver, &Prec); CHKERRQ(ierr);
    // ierr = PCSetType(Prec, PCJACOBI); CHKERRQ(ierr);
    ierr = KSPSetTolerances(Solver, PETSC_DEFAULT, ERROR_THRESHOLD, PETSC_DEFAULT, MAX_ITERATIONS); CHKERRQ(ierr);
    ierr = KSPSetUp(Solver); CHKERRQ(ierr);
    VecView(Sol, 0);
    ierr = KSPSolve(Solver, RHS, Sol); CHKERRQ(ierr);
    VecView(Sol, 0);
    timer.EndTimer("Iterative solution");
    












    ierr = PetscFinalize(); CHKERRQ(ierr);
    return ierr;
    // double* fake_output;
    // return fake_output;
} 