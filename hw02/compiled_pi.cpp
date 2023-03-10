#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<cmath>
#include<chrono>

double RandomDraw() {
	double num = ((double)rand() / (double)(RAND_MAX));
	return num;
}

int main(int argc, char *argv[]) {
	auto start = std::chrono::high_resolution_clock::now();
	srand(time(0));
	const double pi = 3.141592653;
	int num_loops = atoi(argv[1]);
	int ctr_inside = 0;
	int ctr_outside = 0;
	for (int i = 0; i < num_loops; i++){
		double x = RandomDraw();
		double y = RandomDraw();
		double dist = x*x + y*y;
		if (dist <= 1){
			ctr_inside++;
		} else {
			ctr_outside++;
		}
	}
	double pi_est = 4.0000000 * ctr_inside / num_loops;
	double error = pi_est - pi;
	error = fabs(error) / pi;
	std::cout.precision(12);
	auto end = std::chrono::high_resolution_clock::now();
	auto time = std::chrono::duration_cast<std::chrono::seconds>( end - start ).count();
	std::cout << num_loops << ' ' << ctr_inside << ' ' << ctr_outside << ' ' << pi_est << ' ' << error << ' ' << time << '\n';
	return 0;
}
	
