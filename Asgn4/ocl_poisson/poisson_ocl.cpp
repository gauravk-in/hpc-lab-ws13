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

#define BLOCK_SIZE 32

#if defined __APPLE__ || defined(MACOSX)
    #include <OpenCL/opencl.h>
#else
    #include <CL/opencl.h>
#endif

//debugging variables
cl_int err;
cl_event event;

// list of platforms
cl_platform_id* platforms;
cl_uint numPlatforms;
int ch_plat;

// list of devices
cl_device_id* devices;
cl_uint numDevices;
int ch_dev;

cl_context context;
cl_command_queue command_queue;
cl_program program;
cl_kernel g_product_kernel;

size_t globalWorkSize[3];
size_t localWorkSize[3];

void selectDevice()
{
	int i;
	// Get number of platforms
    err = clGetPlatformIDs(0, NULL, &numPlatforms);
    //// printf("clGetPlatformIDs: %s\n", oclErrorString(err));
    printf("num of platforms = %u\n", (unsigned int)numPlatforms);

    // Get platform list
    platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * numPlatforms);
    err = clGetPlatformIDs(numPlatforms, platforms, NULL);
    //// printf("clGetPlatformIDs: %s\n", oclErrorString(err));

	char platform_name[20];
	for(i = 0; i < numPlatforms; i++)
	{
		err = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 20, (void*)platform_name, NULL);
		printf("platforms[%d] Name = %s\n", i, platform_name);
	}
	
	if(numPlatforms == 1)
	{
		printf("Using the only platform available\n");
		ch_plat = 0;
	}
	else
	{
		printf("\nWhich platform do you want to use? ");
		scanf("%d", &ch_plat);
	}

	// Get number of devices
	err = clGetDeviceIDs(platforms[ch_plat], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	//// printf("clGetDeviceIDs (get number of devices): %s\n", oclErrorString(err));
	printf("\nnum of devices = %u\n", (unsigned int)numDevices);

	// Get devices list
    devices = (cl_device_id*) malloc(sizeof(cl_device_id) * numDevices);
	err = clGetDeviceIDs(platforms[ch_plat], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
	//// printf("clGetDeviceIDs (create device list): %s\n", oclErrorString(err));

	// Print Device names
    {
    	for( i = 0; i < numDevices; i++)
    	{
	    	char device_name[20];
	    	int comp_avail;
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 20, (void*)device_name, NULL);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_COMPILER_AVAILABLE, 20, (void*)&comp_avail, NULL);
    		printf("devices[%d] Name = %s; Compiler = %s\n", i, device_name, comp_avail ? "Available" : "Not available");
    	}
	}

	if(numDevices == 0)
	{
		printf("No Device Available. Bye!\n");
		exit(1);
	}
	if(numDevices == 1)
	{
		printf("Using the only device available for this platform.\n");
		ch_dev = 0;
	}
	else
	{
		printf("\nWhich Device do you want to use? ");
		scanf("%d", &ch_dev);
	}

	return;
}

char* readProgramToBuffer(const char *path, size_t *program_length)
{
	FILE *fp = fopen(path, "r");
	char *buffer;

	if(!fp)
	{
		fprintf(stderr, "Unable to open %s for reading!\n", path);
		exit(1);
	}

	fseek(fp, 0, SEEK_END);
	*program_length = ftell(fp);
	fseek(fp, 0, SEEK_SET);

	buffer = (char *)malloc(*program_length * sizeof(char) + 1);
	*program_length = fread(buffer, 1, *program_length, fp);
	fclose(fp);
	buffer[*program_length] = '\0';

	return buffer;
}

void printBuildLog()
{
	cl_build_status build_status;
	err = clGetProgramBuildInfo(program, devices[ch_dev], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL);

	char *build_log;
	size_t ret_val_size;
	err = clGetProgramBuildInfo(program, devices[ch_dev], CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);

	build_log = (char *)malloc(ret_val_size+1);
	err = clGetProgramBuildInfo(program, devices[ch_dev], CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);
	build_log[ret_val_size] = '\0';
	printf("BUILD LOG: \n %s", build_log);
}

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
float timer_stop()
{
	gettimeofday(&end,(struct timezone *)0);
	float seconds, useconds;
	float ret, tmp;

	if (end.tv_usec >= begin.tv_usec)
	{
		seconds = (float)end.tv_sec - (float)begin.tv_sec;
		useconds = (float)end.tv_usec - (float)begin.tv_usec;
	}
	else
	{
		seconds = (float)end.tv_sec - (float)begin.tv_sec;
		seconds -= 1;					// Correction
		useconds = (float)end.tv_usec - (float)begin.tv_usec;
		useconds += 1000000;			// Correction
	}

	// get time in seconds
	tmp = (float)useconds;
	ret = (float)seconds;
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
void store_grid(float* grid, std::string filename)
{
	std::fstream filestr;
	filestr.open (filename.c_str(), std::fstream::out);
	
	// calculate mesh width 
	float mesh_width = 1.0/((float)(grid_points_1d-1));

	// store grid incl. boundary points
	for (int i = 0; i < grid_points_1d; i++)
	{
		for (int j = 0; j < grid_points_1d; j++)
		{
			filestr << mesh_width*i << " " << mesh_width*j << " " << grid[(i*grid_points_1d)+j] << std::endl;
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
float eval_init_func(float x, float y)
{
	return (x*x)*(y*y);
}

/**
 * initializes a given grid: inner points are set to zero
 * boundary points are initialized by calling eval_init_func
 *
 * @param grid the grid to be initialized
 */
void init_grid(float* grid)
{
	// set all points to zero
	for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
	{
		grid[i] = 0.0;
	}

	float mesh_width = 1.0/((float)(grid_points_1d-1));
	
	for (int i = 0; i < grid_points_1d; i++)
	{
		// x-boundaries
		grid[i] = eval_init_func(0.0, ((float)i)*mesh_width);
		grid[i + ((grid_points_1d)*(grid_points_1d-1))] = eval_init_func(1.0, ((float)i)*mesh_width);
		// y-boundaries
		grid[i*grid_points_1d] = eval_init_func(((float)i)*mesh_width, 0.0);
		grid[(i*grid_points_1d) + (grid_points_1d-1)] = eval_init_func(((float)i)*mesh_width, 1.0);
	}
}

/**
 * initializes the right hand side, we want to keep it simple and
 * solve the Laplace equation instead of Poisson (-> b=0)
 *
 * @param b the right hand side
 */
void init_b(float* b)
{
	// set all points to zero
	for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
	{
		b[i] = 0.0;
	}
}

/**
 * copies data from one grid to another
 *
 * @param dest destination grid
 * @param src source grid
 */
void g_copy(float* dest, float* src)
{
	for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
	{
		dest[i] = src[i];
	}
}

/**
 * calculates the dot product of the two grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid1 first grid
 * @param grid2 second grid
 */
float g_dot_product(float* grid1, float* grid2)
{
	float tmp = 0.0;

	for (int i = 1; i < grid_points_1d-1; i++)
	{
		for (int j = 1; j < grid_points_1d-1; j++)
		{
			tmp += (grid1[(i*grid_points_1d)+j] * grid2[(i*grid_points_1d)+j]);
		}
	}
	
	return tmp;
}

/**
 * scales a grid by a given scalar (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid grid to be scaled
 * @param scalar scalar which is used to scale to grid
 */
void g_scale(float* grid, float scalar)
{
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		for (int j = 1; j < grid_points_1d-1; j++)
		{
			grid[(i*grid_points_1d)+j] *= scalar;
		}
	}
}

/**
 * implements BLAS's Xaxpy operation for grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param dest destination grid
 * @param src source grid
 * @param scalar scalar to scale to source grid
 */
void g_scale_add(float* dest, float* src, float scalar)
{
	for (int i = 1; i < grid_points_1d-1; i++)
	{
		for (int j = 1; j < grid_points_1d-1; j++)
		{
			dest[(i*grid_points_1d)+j] += (scalar*src[(i*grid_points_1d)+j]);
		}
	}
}

/**
 * implements the the 5-point finite differences stencil (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 * 
 * @param grid grid for which the stencil should be evaluated
 * @param result grid where the stencil's evaluation should be stored
 */
void g_product_operator(float* grid, float* result)
{
	// 2. Allocate Memory on Host and Device
    cl_mem grid_buffer;
    cl_mem result_buffer;
    cl_mem grid_points_1d_buffer;
    grid_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, grid_points_1d * grid_points_1d * sizeof(float), grid, &err);
    result_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, grid_points_1d * grid_points_1d * sizeof(float), result, &err);

    // 3. Invoke Kernel
    err = clSetKernelArg(g_product_kernel, 0, sizeof(cl_mem), (void *) &grid_buffer);
    // err = clSetKernelArg(g_product_kernel, 1, sizeof(float) * BLOCK_SIZE * BLOCK_SIZE, 0);
    err = clSetKernelArg(g_product_kernel, 1, sizeof(cl_mem), (void *) &result_buffer);
    // err = clSetKernelArg(g_product_kernel, 3, sizeof(float) * BLOCK_SIZE * BLOCK_SIZE, 0);
    err = clSetKernelArg(g_product_kernel, 2, sizeof(cl_ulong), (void *) &grid_points_1d);
    clFinish(command_queue);
    globalWorkSize[0] = grid_points_1d-2;
    globalWorkSize[1] = grid_points_1d-2;
    localWorkSize[0] = BLOCK_SIZE;
    localWorkSize[1] = BLOCK_SIZE;
    err = clEnqueueNDRangeKernel(command_queue, g_product_kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
   	clFinish(command_queue);

   	// 4. Copy result from device
    err = clEnqueueReadBuffer(command_queue, result_buffer, CL_TRUE, 0, grid_points_1d * grid_points_1d * sizeof(float), result, 0, NULL, &event);
    clReleaseEvent(event);

    clFinish(command_queue);
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
void solve(float* grid, float* b, std::size_t cg_max_iterations, float cg_eps)
{
	std::cout << "Starting Conjugated Gradients" << std::endl;

	float eps_squared = cg_eps*cg_eps;
	std::size_t needed_iters = 0;

	// define temporal vectors
	float* q = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	float* r = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	float* d = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	float* b_save = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
			
	g_copy(q, grid);
	g_copy(r, grid);
	g_copy(d, grid);
	g_copy(b_save, b);
	
	float delta_0 = 0.0;
	float delta_old = 0.0;
	float delta_new = 0.0;
	float beta = 0.0;
	float a = 0.0;
	float residuum = 0.0;
	
	g_product_operator(grid, d);
	g_scale_add(b, d, -1.0);
	g_copy(r, b);
	g_copy(d, r);

	// calculate starting norm
	delta_new = g_dot_product(r, r);
	delta_0 = delta_new*eps_squared;
	residuum = (delta_0/eps_squared);
	
	std::cout << "Starting norm of residuum: " << (delta_0/eps_squared) << std::endl;
	std::cout << "Target norm:               " << (delta_0) << std::endl;

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
		std::cout << "(iter: " << needed_iters << ")delta: " << delta_new << std::endl;
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
	float mesh_width = atof(argv[1]);
	size_t cg_max_iterations = atoi(argv[2]);
	float cg_eps = atof(argv[3]);

	// calculate grid points per dimension
	grid_points_1d = (std::size_t)(1.0/mesh_width)+1;
	printf("Grid Points 1d = %d\n", grid_points_1d);
	
	
	// initialize the gird and rights hand side
	float* grid = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	float* b = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	init_grid(grid);
	store_grid(grid, "initial_condition.gnuplot");
	init_b(b);
	store_grid(b, "b.gnuplot");
	
	selectDevice();

	context = clCreateContext(0, 1, &devices[ch_dev], NULL, NULL, &err);
	command_queue = clCreateCommandQueue(context, devices[ch_dev], 0, &err);

	// Compiling Kernels
	size_t program_length;
	void *buffer;
	buffer = readProgramToBuffer("g_product_kernel.cl", &program_length);
	program = clCreateProgramWithSource(context, 1, (const char **) &buffer, &program_length, &err);
	free(buffer);
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printBuildLog();
		exit(1);
	}
	g_product_kernel = clCreateKernel(program, "g_product_operator_parallel", &err);

	// solve Poisson equation using CG method
	timer_start();
	solve(grid, b, cg_max_iterations, cg_eps);
	float time = timer_stop();
	store_grid(grid, "solution.gnuplot");
	
	std::cout << std::endl << "Needed time: " << time << " s" << std::endl << std::endl;
	
	_mm_free(grid);
	_mm_free(b);

	return 0;
}

