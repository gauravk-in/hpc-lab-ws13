#include <stdio.h>

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

void readMatrix(char *fileName, float *mat, int width)
{
	FILE *fp = fopen(fileName, "r");
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			fscanf(fp, "%f ", &mat[i * width + j]);
		}
	}

	fclose(fp);
	return;
}

void writeInputMatrix(char *fileName, float value, int width)
{
	FILE *fp = fopen(fileName, "w");
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			fprintf(fp, "%f ", value);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
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

void printMatrix(float *mat, int width)
{
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			printf("%f ", mat[i * width + j]);
		}
		printf("\n");
	}
	return;
}

void fillMatrix(float *mat, int width, float value)
{
	int i, j;

	for(i = 0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			mat[i*width+j] = value;
		}
	}
}

int main(int argc, char **argv)
{
	size_t size;
    float *M_host;
    float *N_host;
    float *P_host;
    char *M_fileName;
    char *N_fileName;
    char *P_fileName;
    int width;
    size_t globalWorkSize[3];
    size_t localWorkSize[3];
    
    struct timeval tv_start;
    struct timeval tv_end;
    unsigned long us;
    struct timeval tv_start_kernel;
    struct timeval tv_end_kernel;
    unsigned long us_kernel;
    
    if(argc != 5)
    {
    	printf("Usage : %s <M_matrix.txt> <N_matrix.txt> <P_matrix.txt> <width>\n", argv[0]);
    	exit(1);
    }

    sscanf(argv[4], "%d", &width);

    size = width * width * sizeof(float);
	M_host = (float *) malloc(size);
	N_host = (float *) malloc(size);
	P_host = (float *) malloc(size);

	M_fileName = argv[1];
	N_fileName = argv[2];
	P_fileName = argv[3];

    readMatrix(M_fileName, M_host, width);
    readMatrix(N_fileName, N_host, width);

	printf("\n***\nInitialization of OpenCL object and Context Begins!\n***\n");

	// 1. Create Context and Command Queue
	printf(" Select platform and device to use\n");
	selectDevice();

	// Start Time Here
	gettimeofday(&tv_start, NULL);

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

    gettimeofday(&tv_start_kernel, NULL);
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
    err = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
   	clFinish(command_queue);

   	gettimeofday(&tv_end_kernel, NULL);

    // 4. Copy result from device
    // printf("\n***\nReading Result from Device Memory\n***\n");
    err = clEnqueueReadBuffer(command_queue, P_buffer, CL_TRUE, 0, width * width * sizeof(float), P_host, 0, NULL, &event);
    clReleaseEvent(event);

    // Stop Time Here
    gettimeofday(&tv_end, NULL);

    us = (tv_end.tv_sec - tv_start.tv_sec) * 1000000;
    us += ((tv_end.tv_usec - tv_start.tv_usec));

    us_kernel = (tv_end_kernel.tv_sec - tv_start_kernel.tv_sec) * 1000000;
    us_kernel += ((tv_end_kernel.tv_usec - tv_start_kernel.tv_usec));

    printf("\nTime taken with memory, and kernel compilation overhead = %lu usec\n", us);
    printf("\nTime taken for only kernel = %lu usec\n", us_kernel);

    writeMatrix(P_fileName, P_host, width);

    // printf("\n\n*** Successful, Result Stored in %s ***\n", P_fileName);

    // 4. Clean Up!

    return 0;
}