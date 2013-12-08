/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the training material of the master's course           *
* Scientific Computing                                                        *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <x86intrin.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <mpi.h>

#define GRID_DIM_X 4
#define GRID_DIM_Y 4

#define HERE0() if(rank==0) printf("%s: %s: %d\n", __FILE__,__func__,__LINE__)
#define THD_0() if(rank == 0)

/// store number of grid points in one dimension
std::size_t grid_points_1d = 0;
std::size_t grid_points_1d_thread = 0;

/// store begin timestep
struct timeval begin;
/// store end timestep
struct timeval end;

// For creating MPI Topology
MPI_Comm grid_comm;
int dim_sizes[2] = {GRID_DIM_X, GRID_DIM_Y};
int wrap_around[2] = {0, 0};
int reorder = 1;
int size;
int rank;
int coordinates[2];
MPI_Status status;

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
void store_grid(double* grid, std::string filename)
{
	std::fstream filestr;
	int res_i;
	int res_j;
	if(rank == 0)
	{
		filestr.open (filename.c_str(), std::fstream::out);
		filestr.close();
	}
	
	// calculate mesh width 
	double mesh_width = 1.0/((double)(grid_points_1d-1));

	int res_rank;

	// store grid incl. boundary points
	for (int i = 0; i < grid_points_1d; i++)
	{
		for (int j = 0; j < grid_points_1d; j++)
		{
			res_rank = (((i-1) / GRID_DIM_X) * GRID_DIM_X) + (j-1) / GRID_DIM_Y;
			THD_0() printf("i=%d, j=%d, res_rank=%d\n", i, j, res_rank);
			if(rank == res_rank) 
			{
				printf("rank = %d\n", rank);
				res_i = i % GRID_DIM_X;
				res_j = j % GRID_DIM_Y;
				printf("res_i = %d, res_j = %d\n", res_i, res_j);
				filestr.open (filename.c_str(), std::fstream::app | std::fstream::out);
				filestr << mesh_width*i << " " << mesh_width*j << " " << grid[(res_i*grid_points_1d_thread)+res_j] << std::endl;
				filestr.close();
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		
		if(rank == 0)
		{
			filestr.open (filename.c_str(), std::fstream::app | std::fstream::out);
			filestr << std::endl;
			filestr.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
/**
 * calculate the grid's initial values for given grid points
 *
 * @param x the x-coordinate of a given grid point
 * @param y the y-coordinate of a given grid point
 *
 * @return the initial value at position (x,y)
 */
double eval_init_func(double x, double y)
{
	return (x*x)*(y*y);
}

/**
 * initializes a given grid: inner points are set to zero
 * boundary points are initialized by calling eval_init_func
 *
 * @param grid the grid to be initialized
 */
void init_grid(double* grid)
{
	// set all points to zero
	for (int i = 1; i < grid_points_1d_thread*grid_points_1d_thread; i++)
	{
		grid[i] = 0.0;
	}

	double mesh_width = 1.0/((double)(grid_points_1d-1));
	
	// x-boundaries
	if(coordinates[0] == 0)
	{
		for (int i = 1; i < grid_points_1d_thread; i++)
		{
			grid[i*grid_points_1d_thread] = eval_init_func(0.0, ((double)(coordinates[1] * (grid_points_1d_thread - 1) + i))*mesh_width);
		}
	}
	else if(coordinates[0] == GRID_DIM_X-1)
	{
		for (int i = 1; i < grid_points_1d_thread; i++)
		{
			grid[grid_points_1d_thread*(i+1)-1] = eval_init_func(1.0, ((double)(coordinates[1]*(grid_points_1d_thread-1) + i)*mesh_width));
		}
	}
	// y-boundaries	
	else if(coordinates[1] == 0)
	{
		for (int i = 1; i < grid_points_1d_thread; i++)
		{
			grid[i] = eval_init_func(((double)(coordinates[1] * (grid_points_1d_thread - 1) + i))*mesh_width, 0.0);
		}
	}
	else if(coordinates[1] == GRID_DIM_Y-1)
	{
		for (int i = 1; i < grid_points_1d_thread; i++)
		{
			grid[(grid_points_1d_thread)*(grid_points_1d_thread-1)+i] = eval_init_func(((double)(coordinates[1] * (grid_points_1d_thread - 1) + i))*mesh_width, 1.0);
			if(rank == 3) printf("| %f |", grid[(grid_points_1d_thread)*(grid_points_1d_thread-1)+i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * initializes the right hand side, we want to keep it simple and
 * solve the Laplace equation instead of Poisson (-> b=0)
 *
 * @param b the right hand side
 */
void init_b(double* b)
{
	// set all points to zero
	for (int i = 0; i < grid_points_1d_thread*grid_points_1d_thread; i++)
	{
		b[i] = 0.0;
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * copies data from one grid to another
 *
 * @param dest destination grid
 * @param src source grid
 */
void g_copy(double* dest, double* src)
{
	for (int i = 0; i < grid_points_1d_thread*grid_points_1d_thread; i++)
	{
		dest[i] = src[i];
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * calculates the dot product of the two grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid1 first grid
 * @param grid2 second grid
 */
double g_dot_product(double* grid1, double* grid2)
{
	double tmp = 0.0;
	double dot_prod = 0.0;

	for (int i = 1; i < grid_points_1d_thread-1; i++)
	{
		for (int j = 1; j < grid_points_1d_thread-1; j++)
		{
			tmp += (grid1[(i*grid_points_1d_thread)+j] * grid2[(i*grid_points_1d_thread)+j]);
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Reduce the sum of tmp in each thread and return it.
	MPI_Allreduce((void*)&tmp, (void*)&dot_prod, 1, MPI_DOUBLE_PRECISION, MPI_SUM, grid_comm);

	return dot_prod;
}

/**
 * scales a grid by a given scalar (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid grid to be scaled
 * @param scalar scalar which is used to scale to grid
 */
void g_scale(double* grid, double scalar)
{
	for (int i = 1; i < grid_points_1d_thread-1; i++)
	{
		for (int j = 1; j < grid_points_1d_thread-1; j++)
		{
			grid[(i*grid_points_1d_thread)+j] *= scalar;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);	
}

/**
 * implements BLAS's Xaxpy operation for grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param dest destination grid
 * @param src source grid
 * @param scalar scalar to scale to source grid
 */
void g_scale_add(double* dest, double* src, double scalar)
{
	for (int i = 1; i < grid_points_1d_thread-1; i++)
	{
		for (int j = 1; j < grid_points_1d_thread-1; j++)
		{
			dest[(i*grid_points_1d_thread)+j] += (scalar*src[(i*grid_points_1d_thread)+j]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * implements the the 5-point finite differences stencil (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 * 
 * @param grid grid for which the stencil should be evaluated
 * @param result grid where the stencil's evaluation should be stored
 */
void g_product_operator(double* grid, double* result)
{
	double mesh_width = 1.0/((double)(grid_points_1d-1));
	MPI_Datatype vector_lr;

	// All threads Send and receive the required border values.
	// int MPI_Type_vector(int count, int blocklength, int stride,
 	// 		MPI_Datatype oldtype, MPI_Datatype *newtype)
	MPI_Type_vector(grid_points_1d_thread-2, 1, grid_points_1d_thread, MPI_DOUBLE_PRECISION, &vector_lr);
	MPI_Type_commit(&vector_lr);

	if(coordinates[0] != GRID_DIM_X-1)
	{
		// int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
  		//            MPI_Comm comm)
		MPI_Send(&grid[2*grid_points_1d_thread-2], 1, vector_lr, (coordinates[0]+1)*GRID_DIM_X+coordinates[1], 1, MPI_COMM_WORLD);
		MPI_Recv(&grid[2*grid_points_1d_thread-1], 1, vector_lr, (coordinates[0]+1)*GRID_DIM_X+coordinates[1], 1, MPI_COMM_WORLD, &status);
	}
	if(coordinates[0] != 0)
	{
		MPI_Recv(&grid[grid_points_1d_thread], 1, vector_lr, (coordinates[0]-1)*GRID_DIM_X+coordinates[1], 1, MPI_COMM_WORLD, &status);
		MPI_Send(&grid[grid_points_1d_thread+1], 1, vector_lr, (coordinates[0]-1)*GRID_DIM_X+coordinates[1], 1, MPI_COMM_WORLD);
	}

	if(coordinates[1] != GRID_DIM_Y-1)
	{
		// int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
  		//            MPI_Comm comm)
		MPI_Send(&grid[grid_points_1d_thread*(grid_points_1d_thread-2)+1], grid_points_1d_thread-2, MPI_DOUBLE_PRECISION, 
				coordinates[0]*GRID_DIM_X+coordinates[1]+1, 1, MPI_COMM_WORLD);
		MPI_Recv(&grid[grid_points_1d_thread*(grid_points_1d_thread-1)+1], grid_points_1d_thread-2, MPI_DOUBLE_PRECISION, 
				coordinates[0]*GRID_DIM_X+coordinates[1]+1, 1, MPI_COMM_WORLD, &status);
	}
	if(coordinates[1] != 0)
	{
		MPI_Recv(&grid[1], grid_points_1d_thread-2, MPI_DOUBLE_PRECISION, 
				coordinates[0]*GRID_DIM_X+coordinates[1]-1, 1, MPI_COMM_WORLD, &status);
		MPI_Send(&grid[grid_points_1d_thread+1], grid_points_1d_thread-2, MPI_DOUBLE_PRECISION, 
				coordinates[0]*GRID_DIM_X+coordinates[1]-1, 1, MPI_COMM_WORLD);
	}

	HERE0();

	for (int i = 1; i < grid_points_1d_thread; i++)
	{
		for (int j = 1; j < grid_points_1d_thread; j++)
		{
			result[(i*grid_points_1d_thread)+j] =  (
							(4.0*grid[(i*grid_points_1d_thread)+j]) 
							- grid[((i+1)*grid_points_1d_thread)+j]
							- grid[((i-1)*grid_points_1d_thread)+j]
							- grid[(i*grid_points_1d_thread)+j+1]
							- grid[(i*grid_points_1d_thread)+j-1]
							) * (mesh_width*mesh_width);
		}	
	}

	MPI_Barrier(MPI_COMM_WORLD);
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
std::size_t solve(double* grid, double* b, std::size_t cg_max_iterations, double cg_eps)
{
	THD_0() std::cout << "Starting Conjugated Gradients" << std::endl;

	HERE0();

	double eps_squared = cg_eps*cg_eps;
	std::size_t needed_iters = 0;

	// define temporal vectors
	double* q = (double*)_mm_malloc(grid_points_1d_thread*grid_points_1d_thread*sizeof(double), 64);
	double* r = (double*)_mm_malloc(grid_points_1d_thread*grid_points_1d_thread*sizeof(double), 64);
	double* d = (double*)_mm_malloc(grid_points_1d_thread*grid_points_1d_thread*sizeof(double), 64);
	double* b_save = (double*)_mm_malloc(grid_points_1d_thread*grid_points_1d_thread*sizeof(double), 64);
			
	g_copy(q, grid);
	g_copy(r, grid);
	g_copy(d, grid);
	g_copy(b_save, b);
	
	double delta_0 = 0.0;
	double delta_old = 0.0;
	double delta_new = 0.0;
	double beta = 0.0;
	double a = 0.0;
	double residuum = 0.0;
	
	g_product_operator(grid, d);
	g_scale_add(b, d, -1.0);
	g_copy(r, b);
	g_copy(d, r);

	// calculate starting norm
	delta_new = g_dot_product(r, r);
	delta_0 = delta_new*eps_squared;
	residuum = (delta_0/eps_squared);
	
	THD_0() std::cout << "Starting norm of residuum: " << (delta_0/eps_squared) << std::endl;
	THD_0() std::cout << "Target norm:               " << (delta_0) << std::endl;

	while ((needed_iters < cg_max_iterations) && (delta_new > delta_0))
	{
		// q = A*d
		g_product_operator(d, q);

		// a = d_new / d.q
		a = delta_new/g_dot_product(d, q);
		
		// x = x + a*d
		g_scale_add(grid, d, a);
		
		if ((needed_iters % 50) == 0)
		{
			g_copy(b, b_save);
			g_product_operator(grid, q);
			g_scale_add(b, q, -1.0);
			g_copy(r, b);
		}
		else
		{
			// r = r - a*q
			g_scale_add(r, q, -a);
		}
		
		// calculate new deltas and determine beta
		delta_old = delta_new;
		delta_new = g_dot_product(r, r);
		beta = delta_new/delta_old;

		// adjust d
		g_scale(d, beta);
		g_scale_add(d, r, 1.0);
		
		residuum = delta_new;
		needed_iters++;
		THD_0() std::cout << "(iter: " << needed_iters << ")delta: " << delta_new << std::endl;
	}

	THD_0() std::cout << "Number of iterations: " << needed_iters << " (max. " << cg_max_iterations << ")" << std::endl;
	THD_0() std::cout << "Final norm of residuum: " << delta_new << std::endl;
	
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
	// For creating MPI Topology
	// MPI_Comm grid_comm;
	// int dim_sizes[2] = {GRID_DIM_X, GRID_DIM_Y};
	// int wrap_around[2] = {0, 0};
	// int reorder = 1;
	// int size;
	// int rank;
	// int coordinates[2];

	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	MPI_Cart_create (MPI_COMM_WORLD, 2, dim_sizes,
			wrap_around, reorder, &grid_comm);

	MPI_Comm_rank (grid_comm, &rank);
	MPI_Cart_coords (grid_comm, rank, 2, coordinates);

	// check if all parameters are specified
	if (argc != 4 && rank == 0)
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
	// Divide the inner points into grids and keep space for one edge on each border.
	grid_points_1d_thread = (grid_points_1d-2) / GRID_DIM_X + 2;
	
	// initialize the gird and rights hand side
	double* grid = (double*)_mm_malloc(grid_points_1d_thread*grid_points_1d_thread*sizeof(double), 64);
	double* b = (double*)_mm_malloc(grid_points_1d_thread*grid_points_1d_thread*sizeof(double), 64);
	init_grid(grid);
	store_grid(grid, "initial_condition1.gnuplot");
	init_b(b);
	store_grid(b, "b1.gnuplot");
	
	// solve Poisson equation using CG method
	timer_start();
	solve(grid, b, cg_max_iterations, cg_eps);
	double time = timer_stop();
	store_grid(grid, "solution1.gnuplot");
	
	THD_0() std::cout << std::endl << "Needed time: " << time << " s" << std::endl << std::endl;
	
	_mm_free(grid);
	_mm_free(b);

	return 0;
}

