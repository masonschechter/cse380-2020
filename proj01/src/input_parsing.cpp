#include <grvy.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
using namespace GRVY;



struct settings
{
	int ORDER, DISCRETIZATION_POINTS, DIMENSIONS, MAX_ITERATIONS;
	double K_0, A_X, B_Y, X_MIN, X_MAX, Y_MIN, Y_MAX, ERROR_THRESHOLD;
	string SOLVER_TYPE, OUTPUT_MODE, RUN_MODE;
};

settings parse_inputs() {
	GRVY_Input_Class iparse; // Input parsing object
	settings inputs;
	
	if (! iparse.Open("./input.dat")) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not open input file \n");
		exit(1);
	}

	if(! iparse.Read_Var("ORDER",&inputs.ORDER) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: ORDER \n");
		exit(1);
	} else {
		if (inputs.ORDER != 2 && inputs.ORDER != 4) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: ORDER. Use 2 or 4 \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("DISCRETIZATION_POINTS",&inputs.DISCRETIZATION_POINTS) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: DISCRETIZATION_POINTS \n");
		exit(1);
	} else {
		if (inputs.DISCRETIZATION_POINTS < 1) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: DISCRETIZATION_POINTS. Value must be non-zero and positive \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("DIMENSIONS",&inputs.DIMENSIONS) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: DIMENSIONS \n");
		exit(1);
	} else {
		if (inputs.DIMENSIONS != 1 && inputs.DIMENSIONS != 2) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: DIMENSIONS. Use 1 or 2 \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("MAX_ITERATIONS",&inputs.MAX_ITERATIONS) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: MAX_ITERATIONS \n");
		exit(1);
	} else {
		if (inputs.MAX_ITERATIONS <= 0) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: MAX_ITERATIONS. Value must be non-zero and positive \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("K_0",&inputs.K_0) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: K_0 \n");
		exit(1);
	} else {
		if (inputs.K_0 <= 0) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: K_0. Value must be non-zero and positive \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("A_X",&inputs.A_X) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: A_X \n");
		exit(1);
	} else {
		if (inputs.A_X <= 0) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: A_X. Value must be non-zero and positive \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("B_Y",&inputs.B_Y) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: B_Y \n");
		exit(1);
	} else {
		if (inputs.B_Y <= 0) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: Value must be non-zero and positive \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("X_MAX",&inputs.X_MAX) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: X_MAX \n");
		exit(1);
	} 

	if(! iparse.Read_Var("X_MIN",&inputs.X_MIN) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: X_MIN \n");
		exit(1);
	} else {
		if (inputs.X_MIN > inputs.X_MAX) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported values for variables: X_MIN and X_MAX. Ensure that X_MIN < X_MAX \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("Y_MAX",&inputs.Y_MAX) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: Y_MAX \n");
		exit(1);
	} 

	if(! iparse.Read_Var("Y_MIN",&inputs.Y_MIN) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: Y_MIN \n");
		exit(1);
	} else {
		if (inputs.Y_MIN > inputs.Y_MAX) {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported values for variables: Y_MIN and Y_MAX. Ensure that Y_MIN < Y_MAX \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("ERROR_THRESHOLD",&inputs.ERROR_THRESHOLD) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: ERROR_THRESHOLD \n");
		exit(1);
	}

	if(! iparse.Read_Var("SOLVER_TYPE",&inputs.SOLVER_TYPE) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: SOLVER_TYPE \n");
		exit(1);
	} else {
		if (inputs.SOLVER_TYPE != "Jacobi" && inputs.SOLVER_TYPE != "Gauss-Seidel") {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: SOLVER_TYPE. Use 'Jacobi' or 'Gauss-Seidel' \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("OUTPUT_MODE",&inputs.OUTPUT_MODE) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: OUTPUT_MODE \n");
		exit(1);
	} else {
		if (inputs.OUTPUT_MODE != "INFO" && inputs.OUTPUT_MODE != "DEBUG") {
			grvy_printf(GRVY_ERROR, "ERROR: Unsupported value for variable: OUTPUT_MODE. Use 'INFO' or 'DEBUG' \n");
			exit(1);
		}
	}

	if(! iparse.Read_Var("RUN_MODE",&inputs.RUN_MODE) ) {
		grvy_printf(GRVY_ERROR, "ERROR: Could not find variable: RUN_MODE \n");
		exit(1);
	}

	iparse.Fdump("# ");
	iparse.Close();
	return inputs;
}