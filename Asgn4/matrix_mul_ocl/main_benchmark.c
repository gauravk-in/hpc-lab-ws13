#include <stdio.h>
//#include <sys/time.h>
#include "timer.h"

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
cl_kernel kernel;

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

	buffer = malloc(*program_length * sizeof(char) + 1);
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

	build_log = malloc(ret_val_size+1);
	err = clGetProgramBuildInfo(program, devices[ch_dev], CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);
	build_log[ret_val_size] = '\0';
	printf("BUILD LOG: \n %s", build_log);
}

void writeMatrix(char *fileName, float *mat, int width)
{
	FILE *fp = fopen(fileName, "w");
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			fprintf(fp, "%f ", mat[i * width + j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

int main(int argc, char **argv)
{
    size_t size;
    float *M_host;
    float *N_host;
    float *P_host;
    int width;
    int i, j;
    size_t globalWorkSize[3];
    size_t localWorkSize[3];
    
    time_marker_t time_full;
    time_marker_t time_kernel;
    float flops;
    
    if(argc == 1)
    {
    	width = 1024;
    }
    if(argc > 1)
    {
    	sscanf(argv[1], "%d", &width);
    }
    if(argc > 2)
    {
    	printf("Usage : %s <width>\n\tInitializes input matrices and performes multiplication", argv[0]);
    	exit(1);
    }

    flops = 2.0 * width * width * width / 1000000;

    size = width * width * sizeof(float);
	M_host = (float *) malloc(size);
	N_host = (float *) malloc(size);
	P_host = (float *) malloc(size);

    for (i = 0; i < width; i++){
            for (j = 0; j < width; j++){
                    *(M_host + i * width + j) = (float)i + (float)j;
                    *(N_host + i * width + j) = (float)(width - i) + (float)(width - j);
    	}      
    }
    memset(P_host, 0, size);

	// 1. Create Context and Command Queue
	selectDevice();

	// Start Time Here to consider overhead.
	time_full = get_time();

	context = clCreateContext(0, 1, &devices[ch_dev], NULL, NULL, &err);
	command_queue = clCreateCommandQueue(context, devices[ch_dev], 0, &err);

	// 1. Compile Kernel
	// printf("\n***\nCompiling Kernel for your device!\n***\n");
	size_t program_length;
	void *buffer;
	buffer = readProgramToBuffer("matrix_mul_kernel.cl", &program_length);
	program = clCreateProgramWithSource(context, 1, (const char **) &buffer, &program_length, &err);
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printBuildLog();
		exit(1);
	}
	kernel = clCreateKernel(program, "matrixMultiplyParallel", &err);

    // 2. Allocate Memory on Host and Device
    // printf("\n***\nAllocating Memory on the device!\n***\n");
    cl_mem M_buffer;
    cl_mem N_buffer;
    cl_mem P_buffer;
    cl_mem width_buffer;
    M_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, width * width * sizeof(float), M_host, &err);
    N_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, width * width * sizeof(float), N_host, &err);
    P_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, width * width * sizeof(float), NULL, &err);

    // 3. Invoke Kernel
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &M_buffer);
    err = clSetKernelArg(kernel, 1, sizeof(float) * BLOCK_SIZE * BLOCK_SIZE, 0);
    err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &N_buffer);
    err = clSetKernelArg(kernel, 3, sizeof(float) * BLOCK_SIZE * BLOCK_SIZE, 0);
    err = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &P_buffer);
    err = clSetKernelArg(kernel, 5, sizeof(cl_int), (void *) &width);
    clFinish(command_queue);
    globalWorkSize[0] = width;
    globalWorkSize[1] = width;
    localWorkSize[0] = BLOCK_SIZE;
    localWorkSize[1] = BLOCK_SIZE;
    
    // printf("Global Work Size = %d x %d\nLocal Work Size = %d x %d\n", globalWorkSize[0], globalWorkSize[1], localWorkSize[0], localWorkSize[1]);
    time_kernel = get_time();

    err = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
   	clFinish(command_queue);

   	printf("GFLOP/s for only kernel = ");
   	print_flops(flops, time_kernel);

    // 4. Copy result from device
    err = clEnqueueReadBuffer(command_queue, P_buffer, CL_TRUE, 0, width * width * sizeof(float), P_host, 0, NULL, &event);
    clReleaseEvent(event);

	printf("GFLOP/s with overhead = ");
   	print_flops(flops, time_full);

    return 0;
}
