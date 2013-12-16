/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the training material of the master's course           *
* Scientific Computing                                                        *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <immintrin.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <cilk/cilk.h>

/// store number of grid points in one dimension
std::size_t grid_points_1d = 0;

/// store begin timestep
struct timeval begin;
/// store end timestep
struct timeval end;

/**
 * initialize and start timer
 */
void timer_start()
{
	gettimeofday(&begin,(struct timezone *)0);
}

/**
 * stop timer and return measured time
 *
 * @return measured time
 */
double timer_stop()
{
	gettimeofday(&end,(struct timezone *)0);
	double seconds, useconds;
	double ret, tmp;

	if (end.tv_usec >= begin.tv_usec)
	{
		seconds = (double)end.tv_sec - (double)begin.tv_sec;
		useconds = (double)end.tv_usec - (double)begin.tv_usec;
	}
	else
	{
		seconds = (double)end.tv_sec - (double)begin.tv_sec;
		seconds -= 1;					// Correction
		useconds = (double)end.tv_usec - (double)begin.tv_usec;
		useconds += 1000000;			// Correction
	}

	// get time in seconds
	tmp = (double)useconds;
	ret = (double)seconds;
	tmp /= 1000000;
	ret += tmp;

	return ret;
}

/**
 * stores a given grid into a file
 * 
 * @param grid the grid that should be stored
 * @param filename the filename
 */
void store_grid(double** grid, std::string filename)
{
	std::fstream filestr;
	filestr.open (filename.c_str(), std::fstream::out);
	
	// calculate mesh width 
	double mesh_width = 1.0/((double)(grid_points_1d-1));

	// store grid incl. boundary points
	for (int i = 0; i < grid_points_1d; i++)
	{
		for (int j = 0; j < grid_points_1d; j++)
		{
			filestr << mesh_width*i << " " << mesh_width*j << " " << grid[i][j] << std::endl;
		}
		
		filestr << std::endl;
	}

	filestr.close();
}

/**
 * calculate the grid's initial values for given grid points
 *
 * @param x the x-coordinate of a given grid point
 * @param y the y-coordinate of a given grid point
 *
 * @return the initial value at position (x,y)
 */
__declspec(vector) double eval_init_func(double x, double y)
{
	return (x*x)*(y*y);
}

/**
 * initializes a given grid: inner points are set to zero
 * boundary points are initialized by calling eval_init_func
 *
 * @param grid the grid to be initialized
 */
 // Test
void init_grid_2(double* grid)
{
	double* x = (double*)_mm_malloc(grid_points_1d*sizeof(double), 64);
	double* y = (double*)_mm_malloc(grid_points_1d*sizeof(double), 64);
	double* a = (double*)_mm_malloc(grid_points_1d*sizeof(double), 64);

	grid[0:grid_points_1d*grid_points_1d] = 0.0;

	double mesh_width = 1.0/((double)(grid_points_1d-1));

	a[0] = 0.0;
        a[1:(grid_points_1d-1)] = a[0:(grid_points_1d-1)] + 1.0;

	x[0:grid_points_1d] = 0.0;
	y[0:grid_points_1d] = a[0:grid_points_1d] * mesh_width;

	//x-boundaries upper        
	grid[0:grid_points_1d] = eval_init_func(x[0:grid_points_1d], y[0:grid_points_1d]);        
	//y-boundaries left        
	grid[0:grid_points_1d:grid_points_1d] = eval_init_func(y[0:grid_points_1d], x[0:grid_points_1d]);

	x[0:grid_points_1d] = 1.0;

	//x-boundaries lower        
	grid[(grid_points_1d-1)*(grid_points_1d):grid_points_1d] = eval_init_func(x[0:grid_points_1d], y[0:grid_points_1d]);        
	//y-boundaries right        
	grid[(grid_points_1d-1):grid_points_1d:grid_points_1d] = eval_init_func(y[0:grid_points_1d], x[0:grid_points_1d]);

}

void init_grid(double** grid)
{
	
	double* x = (double*)_mm_malloc(grid_points_1d*sizeof(double), 64); 
	double* y = (double*)_mm_malloc(grid_points_1d*sizeof(double), 64);
	double* a = (double*)_mm_malloc(grid_points_1d*sizeof(double), 64);


	// set all points to zero
	/*
	for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
	{
		grid[i] = 0.0;
	}
	*/
	//grid[0:grid_points_1d*grid_points_1d] = 0.0;
	grid[0:grid_points_1d][0:grid_points_1d] = 0.0;	

	double mesh_width = 1.0/((double)(grid_points_1d-1));
	
	/*	
	for (int i = 0; i < grid_points_1d; i++)
	{
		// x-boundaries
		grid[i] = eval_init_func(0.0, ((double)i)*mesh_width);
		grid[i + ((grid_points_1d)*(grid_points_1d-1))] = eval_init_func(1.0, ((double)i)*mesh_width);
		// y-boundaries
		grid[i*grid_points_1d] = eval_init_func(((double)i)*mesh_width, 0.0);
		grid[(i*grid_points_1d) + (grid_points_1d-1)] = eval_init_func(((double)i)*mesh_width, 1.0);
	}
	*/

	a[0] = 0.0;        
	a[1:(grid_points_1d-1)] = a[0:(grid_points_1d-1)] + 1.0;
	
	/*
	for(int i = 0; i < grid_points_1d; i++)
	{
		printf("a[%d] = %e\n", i, a[i]);
	}
	*/
	
	x[0:grid_points_1d] = 0.0;
	y[0:grid_points_1d] = a[0:grid_points_1d] * mesh_width;

	//x-boundaries upper
	grid[0][0:grid_points_1d] = eval_init_func(x[0:grid_points_1d], y[0:grid_points_1d]);
	//y-boundaries left
	grid[0:grid_points_1d][0] = eval_init_func(y[0:grid_points_1d], x[0:grid_points_1d]);

	x[0:grid_points_1d] = 1.0;

	//x-boundaries lower
	grid[(grid_points_1d)][0:grid_points_1d] = eval_init_func(x[0:grid_points_1d], y[0:grid_points_1d]);
	//y-boundaries right
	grid[0:grid_points_1d][grid_points_1d] = eval_init_func(y[0:grid_points_1d], x[0:grid_points_1d]);
			
}

/**
 * initializes the right hand side, we want to keep it simple and
 * solve the Laplace equation instead of Poisson (-> b=0)
 *
 * @param b the right hand side
 */
 // test
void init_b_2(double* b)
{
	b[0:grid_points_1d*grid_points_1d] = 0.0;
}

void init_b(double** b)
{
	// set all points to zero
	/*
	for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
	{
		b[i] = 0.0;
	}
	*/
	//b[0:grid_points_1d*grid_points_1d] = 0.0;
	b[0:grid_points_1d][0:grid_points_1d] = 0.0;

}

/**
 * copies data from one grid to another
 *
 * @param dest destination grid
 * @param src source grid
 */
 // Test
void g_copy_2(double* dest, double* src)
{
	dest[0:grid_points_1d*grid_points_1d] = src[0:grid_points_1d*grid_points_1d];
}

void g_copy(double** dest, double** src)
{
	/*	
	for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
	{
		dest[i] = src[i];
	}
	*/

	//dest[0:grid_points_1d*grid_points_1d] = src[0:grid_points_1d*grid_points_1d];
	dest[0:grid_points_1d][0:grid_points_1d] = src[0:grid_points_1d][0:grid_points_1d];

}


/**
 * calculates the dot product of the two grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid1 first grid
 * @param grid2 second grid
 */
 // Test
double g_dot_product(double** grid1, double** grid2)
{
	double tmp = 0.0;
	double* mult = (double*)_mm_malloc((grid_points_1d-2)*(grid_points_1d-2)*sizeof(double), 64);
	
	double** mult_2d = (double**)malloc((grid_points_1d-2)*sizeof(double));

	mult_2d[0:grid_points_1d-2] = &(mult[0:grid_points_1d-2:grid_points_1d-2]);

	/*	
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		for (int j = 1; j < grid_points_1d-1; j++)
		{
			tmp += (grid1[(i*grid_points_1d)+j] * grid2[(i*grid_points_1d)+j]);
		}
	}
	*/
	
	/*	
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		mult[(i-1)*(grid_points_1d-2):(grid_points_1d-2)] = 
				grid1[((i*grid_points_1d)+1):(grid_points_1d-2)] * grid2[((i*grid_points_1d)+1):(grid_points_1d-2)];
		
		tmp = __sec_reduce_add(mult[0:((grid_points_1d-2)*(grid_points_1d-2))]);
	}*/

	mult_2d[0:grid_points_1d-2][0:grid_points_1d-2] = grid1[1:grid_points_1d-2][1:grid_points_1d-2] * grid2[1:grid_points_1d-2][1:grid_points_1d-2];

	tmp = __sec_reduce_add(mult_2d[0:grid_points_1d-2][0:grid_points_1d-2]);
	
		
	return tmp;
}

/**
 * scales a grid by a given scalar (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid grid to be scaled
 * @param scalar scalar which is used to scale to grid
 */
 // Test
void g_scale(double** grid, double scalar)
{
	
	int start;
	/*	
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		for (int j = 1; j < grid_points_1d-1; j++)
		{
			grid[(i*grid_points_1d)+j] *= scalar;
		}
	}
	*/
	
	/*
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		start = ((i*grid_points_1d) + 1);
		grid[start:(grid_points_1d-2)] = grid[start:(grid_points_1d-2)] * scalar;
	}
	*/
	grid[1:grid_points_1d-2][1:grid_points_1d-2] = grid[1:grid_points_1d-2][1:grid_points_1d-2] * scalar;

	//grid[(grid_points_1d+1):(grid_points_1d-2)*(grid_points_1d-2)] = grid[(grid_points_1d+1):(grid_points_1d-2)*(grid_points_1d-2)] *  scalar;
}

/**
 * implements BLAS's Xaxpy operation for grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param dest destination grid
 * @param src source grid
 * @param scalar scalar to scale to source grid
 */
 // Test
void g_scale_add(double** dest, double** src, double scalar)
{
	
	int start;
	/*
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		for (int j = 1; j < grid_points_1d-1; j++)
		{
			dest[(i*grid_points_1d)+j] += (scalar*src[(i*grid_points_1d)+j]);
		}
	}
	*/
	
	/*
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		start = ((i*grid_points_1d) + 1);
		dest[start:(grid_points_1d-2)] = dest[start:(grid_points_1d-2)] + (scalar * src[start:(grid_points_1d-2)]); 	
	}
	*/

	dest[1:grid_points_1d-2][1:grid_points_1d-2] = dest[1:grid_points_1d-2][1:grid_points_1d-2] + (scalar * src[1:grid_points_1d-2][1:grid_points_1d-2]);	

	//dest[(grid_points_1d+1):(grid_points_1d-2)*(grid_points_1d-2)] += (scalar*src[(grid_points_1d+1):(grid_points_1d-2)*(grid_points_1d-2)]);
	
}

/**
 * implements the the 5-point finite differences stencil (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 * 
 * @param grid grid for which the stencil should be evaluated
 * @param result grid where the stencil's evaluation should be stored
 */
 // Test
void g_product_operator(double** grid, double** result)
{
	double mesh_width = 1.0/((double)(grid_points_1d-1));

	/*
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		for (int j = 1; j < grid_points_1d-1; j++)
		{
			result[(i*grid_points_1d)+j] =  (
							(4.0*grid[(i*grid_points_1d)+j]) 
							- grid[((i+1)*grid_points_1d)+j]
							- grid[((i-1)*grid_points_1d)+j]
							- grid[(i*grid_points_1d)+j+1]
							- grid[(i*grid_points_1d)+j-1]
							) * (mesh_width*mesh_width);
		}
	}
	*/
	
	result[1:grid_points_1d-2][1:grid_points_1d-2] = (
							(4.0*grid[1:grid_points_1d-2][1:grid_points_1d-2])
							- grid[2:grid_points_1d-2][1:grid_points_1d-2]
							- grid[0:grid_points_1d-2][1:grid_points_1d-2]
							- grid[1:grid_points_1d-2][2:grid_points_1d-2]
							- grid[1:grid_points_1d-2][0:grid_points_1d-2]
							) * (mesh_width*mesh_width);

	/*
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		result[((i*grid_points_1d)+1):(grid_points_1d-2)] = (
								(4.0*grid[((i*grid_points_1d)+1):(grid_points_1d-2)])
								- grid[(((i+1)*grid_points_1d)+1):(grid_points_1d-2)]
								- grid[(((i-1)*grid_points_1d)+1):(grid_points_1d-2)]
								- grid[((i*grid_points_1d)+2):(grid_points_1d-2)]
								- grid[(i*grid_points_1d):(grid_points_1d-2)]
								) * (mesh_width*mesh_width);
	}
	*/
		
}

/**
 * The CG Solver (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * For details please see :
 * http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
 *
 * @param grid the grid containing the initial condition
 * @param b the right hand side
 * @param cg_max_iterations max. number of CG iterations 
 * @param cg_eps the CG's epsilon
 */
void solve(double** grid_2d, double** b_2d, std::size_t cg_max_iterations, double cg_eps)
{
	std::cout << "Starting Conjugated Gradients" << std::endl;

	double eps_squared = cg_eps*cg_eps;
	std::size_t needed_iters = 0;

	// define temporal vectors
	double* q = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
	double* r = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
	double* d = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
	double* b_save = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
	
	// define 2d arrays
	double** q_2d = (double**)malloc(grid_points_1d*sizeof(double*));
	double** r_2d = (double**)malloc(grid_points_1d*sizeof(double*));
	double** d_2d = (double**)malloc(grid_points_1d*sizeof(double*));
	double** b_save_2d = (double**)malloc(grid_points_1d*sizeof(double*));	
	
	q_2d[0:grid_points_1d] = &(q[0:grid_points_1d:grid_points_1d]);
	r_2d[0:grid_points_1d] = &(r[0:grid_points_1d:grid_points_1d]);
	d_2d[0:grid_points_1d] = &(d[0:grid_points_1d:grid_points_1d]);
	b_save_2d[0:grid_points_1d] = &(b_save[0:grid_points_1d:grid_points_1d]);
	
	

	g_copy(q_2d, grid_2d);
	g_copy(r_2d, grid_2d);
	g_copy(d_2d, grid_2d);
	g_copy(b_save_2d, b_2d);
	
	//g_copy2(q_2d, r_2d);

	double delta_0 = 0.0;
	double delta_old = 0.0;
	double delta_new = 0.0;
	double beta = 0.0;
	double a = 0.0;
	double residuum = 0.0;
	
	g_product_operator(grid_2d, d_2d);
	g_scale_add(b_2d, d_2d, -1.0);
	g_copy(r_2d, b_2d);
	g_copy(d_2d, r_2d);

	// calculate starting norm
	delta_new = g_dot_product(r_2d, r_2d);
	delta_0 = delta_new*eps_squared;
	residuum = (delta_0/eps_squared);
	
	//std::cout << "Starting norm of residuum: " << (delta_0/eps_squared) << std::endl;
	//std::cout << "Target norm:               " << (delta_0) << std::endl;

	while ((needed_iters < cg_max_iterations) && (delta_new > delta_0))
	{
		// q = A*d
		g_product_operator(d_2d, q_2d);

		// a = d_new / d.q
		a = delta_new/g_dot_product(d_2d, q_2d);
		
		// x = x + a*d
		g_scale_add(grid_2d, d_2d, a);
		
		if ((needed_iters % 50) == 0)
		{
			g_copy(b_2d, b_save_2d);
			g_product_operator(grid_2d, q_2d);
			g_scale_add(b_2d, q_2d, -1.0);
			g_copy(r_2d, b_2d);
		}
		else
		{
			// r = r - a*q
			g_scale_add(r_2d, q_2d, -a);
		}
		
		// calculate new deltas and determine beta
		delta_old = delta_new;
		delta_new = g_dot_product(r_2d, r_2d);
		beta = delta_new/delta_old;

		// adjust d
		g_scale(d_2d, beta);
		g_scale_add(d_2d, r_2d, 1.0);
		
		residuum = delta_new;
		needed_iters++;
		//std::cout << "(iter: " << needed_iters << ")delta: " << delta_new << std::endl;
	}

	std::cout << "Number of iterations: " << needed_iters << " (max. " << cg_max_iterations << ")" << std::endl;
	std::cout << "Final norm of residuum: " << delta_new << std::endl;
	
	_mm_free(d);
	_mm_free(q);
	_mm_free(r);
	_mm_free(b_save);

	
}

/**
 * main application
 *
 * @param argc number of cli arguments
 * @param argv values of cli arguments
 */
int main(int argc, char* argv[])
{
	// check if all parameters are specified
	if (argc != 4)
	{
		std::cout << std::endl;
		std::cout << "meshwidth" << std::endl;
		std::cout << "cg_max_iterations" << std::endl;
		std::cout << "cg_eps" << std::endl;
		std::cout << std::endl;
		std::cout << "example:" << std::endl;
		std::cout << "./app 0.125 100 0.0001" << std::endl;
		std::cout << std::endl;
		
		return -1;
	}
	
	// read cli arguments
	double mesh_width = atof(argv[1]);
	size_t cg_max_iterations = atoi(argv[2]);
	double cg_eps = atof(argv[3]);

	// calculate grid points per dimension
	grid_points_1d = (std::size_t)(1.0/mesh_width)+1;
	
	// initialize the gird and rights hand side
	double* grid = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
	double* b = (double*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(double), 64);
	
	init_grid_2(grid);

	//initialize 2d grid and array b
	double** grid_2d = (double**)malloc(grid_points_1d*sizeof(double*));
	double** b_2d = (double**)malloc(grid_points_1d*sizeof(double*));

	grid_2d[0:grid_points_1d] = &(grid[0:grid_points_1d:grid_points_1d]);
	/*
	for(int i = 0; i<grid_points_1d; i++){
		for(int j = 0; j<grid_points_1d; j++){
			printf("grid_2d[%d][%d] = %e\n", i, j, grid_2d[i][j]);
		}
	}
	*/
	
	store_grid(grid_2d, "initial_condition.gnuplot");

	init_b_2(b);

	b_2d[0:grid_points_1d] = &(b[0:grid_points_1d:grid_points_1d]);
	
	store_grid(b_2d, "b.gnuplot");
	
	// solve Poisson equation using CG method
	timer_start();
	solve(grid_2d, b_2d, cg_max_iterations, cg_eps);
	double time = timer_stop();
	store_grid(grid_2d, "solution.gnuplot");
	
	std::cout << std::endl << "Needed time: " << time << " s" << std::endl << std::endl;
	
	_mm_free(grid);
	_mm_free(b);

	return 0;
}

