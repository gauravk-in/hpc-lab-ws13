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
cl_kernel reduce_kernel;

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

int main(int argc, char* argv[])
{
	float *input_buffer;
	float reduced_sum[2];
	int len = 10;

	input_buffer = malloc(sizeof(float) * len);
	for(i = 0; i<len; i++)
		input_buffer[i] = 1;

	selectDevice();

	context = clCreateContext(0, 1, &devices[ch_dev], NULL, NULL, &err);
	command_queue = clCreateCommandQueue(context, devices[ch_dev], 0, &err);

	// Compiling Kernels
	size_t program_length;
	void *buffer;
	buffer = readProgramToBuffer("trial_kernel.cl", &program_length);
	program = clCreateProgramWithSource(context, 1, (const char **) &buffer, &program_length, &err);
	free(buffer);
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printBuildLog();
		exit(1);
	}
	reduce_kernel = clCreateKernel(program, "reduce_parallel", &err);

	// 2. Allocate Memory on Host and Device
    cl_mem input_buffer;
    cl_mem reduce_tmp_buffer;
    input_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, len * sizeof(float), input_buffer, &err);
    reduce_tmp_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, len * sizeof(float), NULL, &err);

    // 3. Invoke Kernel
    err = clSetKernelArg(reduce_kernel, 0, sizeof(cl_mem), (void *) &input_buffer);
    err = clSetKernelArg(reduce_kernel, 1, sizeof(cl_mem), (void *) &reduce_tmp_buffer);
    clFinish(command_queue);
    globalWorkSize[0] = 1;
    localWorkSize[0] = BLOCK_SIZE;
    err = clEnqueueNDRangeKernel(command_queue, reduce_kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    printf("err = %d\n", err);
   	clFinish(command_queue);

   	// 4. Copy result from device
    err = clEnqueueReadBuffer(command_queue, reduce_tmp_buffer, CL_TRUE, 0, 2 * sizeof(float), result, 0, NULL, &event);
    printf("err = %d\n", err);
    clReleaseEvent(event);

    clFinish(command_queue);

    printf("Reduced Sum = %f\n", reduced_sum[1]);

	_mm_free(input_buffer);

	return 0;
}
