#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <gsl/gsl_math.h>
using namespace std;

double f(double x){
	// Derivative of erf(x), to find instantaneous y-values
	return M_2_SQRTPI * pow(M_E, -(x*x));
}

double trapezoid(double LOWER_LIM, double UPPER_LIM, int NUM_POINTS){
	// Implementation of the trapezoid rule to estimate an integral. 
	// Arguments: lower limit, upper limit, number of discretization points
	double area;	
	double x_coords[NUM_POINTS];
        double y_coords[NUM_POINTS];
	// Populate arrays of the x- and corresponding y-coordinates for a given number of discretization points
        for (int i = 0; i < NUM_POINTS; i++){
                x_coords[i] = ((UPPER_LIM - LOWER_LIM) / (NUM_POINTS - 1)) * i;
        }
        for (int i = 0; i < NUM_POINTS; i++){
                y_coords[i] = f(x_coords[i]);
        }
	// Calculate and sum all of the trapezoids created in previous step
	for (int i = 0; i < NUM_POINTS - 1; i++){
                double x_1 = x_coords[i];
                double x_2 = x_coords[i + 1];
                double y_1 = y_coords[i];
                double y_2 = y_coords[i + 1];
                area += (y_1 + y_2) / 2 * (x_2 - x_1);
        }
	return area;
}

double simpson(double LOWER_LIM, double UPPER_LIM, int NUM_POINTS){
	// Implementation of Simpson's rule to estimate an integral.
	// Arguments: lower limit, upper limit, number of discretization points
	double area;
	double delta_x = (UPPER_LIM - LOWER_LIM) / (NUM_POINTS - 1);
	double x_coords[NUM_POINTS];
        double y_coords[NUM_POINTS];
	// Populate arrays of the x- and corresponding y-coordinates for a given number of discretization points
        for (int i = 0; i < NUM_POINTS; i++){
                x_coords[i] = delta_x * i;
        }
        for (int i = 0; i < NUM_POINTS; i++){
                y_coords[i] = f(x_coords[i]);
        }
	// Calculate the sum of the odd and even y-values
	double sum_of_odds = 0;
	double sum_of_evens = 0;
	for (int i = 1; i < NUM_POINTS; i += 2){
		sum_of_odds += y_coords[i];
	}
	for (int i = 2; i < NUM_POINTS; i += 2){
                sum_of_evens += y_coords[i];
	}
	// Calculate area using a common re-writing of Simpson's rule
	area = delta_x / 3 * (y_coords[0] + (4 * sum_of_odds) + (2 * sum_of_evens) + y_coords[NUM_POINTS]);
	return area;
}

int main(int argc, char *argv[]){
	// Integration method
	string int_method = argv[1];
	// Number of discretization points
	int NUM_POINTS = atoi(argv[2]);
	double LOWER_LIM = 0;
	double UPPER_LIM = 1;
	double area;
	double error = 0.84270079295;
	if (int_method == "simpson"){
		area = simpson(LOWER_LIM, UPPER_LIM, NUM_POINTS);
	} else if (int_method == "trapezoid"){
		area = trapezoid(LOWER_LIM, UPPER_LIM, NUM_POINTS);
	}
	cout << setprecision(11);
	double std_error = fabs(area - error) / error;
	// Output formatted for downstream applications
	cout << NUM_POINTS  << ' ' << std_error << '\n';
	return 0;
}
