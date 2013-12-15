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
    printf("err = %d\n", err);
   	clFinish(command_queue);

   	// 4. Copy result from device
    err = clEnqueueReadBuffer(command_queue, result_buffer, CL_TRUE, 0, grid_points_1d * grid_points_1d * sizeof(float), result, 0, NULL, &event);
    printf("err = %d\n", err);
    clReleaseEvent(event);

    clFinish(command_queue);

    printf("result[grid_points_1d+2] = %f\n", result[5*grid_points_1d+5]);
}

/**
 * implements the the 5-point finite differences stencil (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 * 
 * @param grid grid for which the stencil should be evaluated
 * @param result grid where the stencil's evaluation should be stored
 */
void g_product_operator_serial(float* grid, float* result)
{
	float mesh_width = 1.0/((float)(grid_points_1d-1));

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


int main(int argc, char* argv[])
{
	if(argc == 1)
		grid_points_1d = 1024;
	if(argc > 1)
		grid_points_1d = atoi(argv[1]);

	float* grid = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	float* result = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	float* result_serial = (float*)_mm_malloc(grid_points_1d*grid_points_1d*sizeof(float), 64);
	init_grid(grid);
	store_grid(grid, "initial_condition.gnuplot");

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

	g_product_operator(grid, result);

	store_grid(result, "result_ocl.gnuplot");

	g_product_operator_serial(grid, result_serial);

	store_grid(result_serial, "result_serial.gnuplot");

	_mm_free(grid);
	_mm_free(result);
	_mm_free(result_serial);

	return 0;
}
