/**
 * Assignment 2: Exercise 5:
 * OpenMP Critical Pi Calculation
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "timer.h"

int main(int argc, char **argv)
{
	double pi = 0.0;
	double x = 0.0;
	unsigned long steps = 1000000;
	double dx= 0.0;
	int i;
	double val;	
	double time_taken;
	time_marker_t time;

	if(argc > 1)
	{
		steps = strtoul(argv[1], NULL, 0);
	}
	dx = 1.0 / steps;

	printf("Step Size = %e\n", dx);	

	time = get_time();

	#pragma omp parallel for schedule(static) private(val) private(x)
	for(i = 0; i < steps; i++)
	{
		//printf("Num Threads = %d\n", omp_get_num_threads());
		x = (i * dx) + (dx / 2);
		val = (1.0 / (1.0 + x * x));
		#pragma omp critical
		{
			pi += val;
		} 

	}

	pi = pi * 4 * dx;

	time_taken = get_ToD_diff_time(time);

	printf("Value of pi = %e\n", pi);
	printf("Time Taken = %f\n", time_taken);

	return 0;
}

