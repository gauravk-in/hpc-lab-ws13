#include <stdio.h>
#include <sys/time.h>


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
	    	cl_uint max_compute_units;
	    	size_t max_work_group_size;
	    	cl_uint max_work_item_dim;
	    	size_t* max_work_item_size;
	    	cl_ulong local_mem_size;
	    	cl_ulong global_mem_size;
	    	int j;

    		err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 20, (void*)device_name, NULL);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_COMPILER_AVAILABLE, 20, (void*)&comp_avail, NULL);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), (void*)&max_compute_units, NULL);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE , sizeof(size_t), (void*)&max_work_group_size, NULL);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), (void*)&max_work_item_dim, NULL);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), (void*)global_mem_size, NULL);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), (void*)local_mem_size, NULL);

    		max_work_item_size = (size_t *) malloc(sizeof(size_t) * (size_t)max_work_item_dim);
    		err = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * (size_t)max_work_item_dim, (void*)max_work_item_size, NULL);
    		printf("devices[%d] Name = %s; Compiler = %s\n", i, device_name, comp_avail ? "Available" : "Not available");
    		printf("Max Compute Units = %d\n", (unsigned int)max_compute_units);
    		printf("Max Work Item Dim = %d\n", (unsigned int)max_work_item_dim);
    		printf("Max Work Item Sizes = ");
    		for(j=0; j<max_work_item_dim; j++)
    			printf("%u ", max_work_item_size[j]);
    		printf("\n");
    		printf("Global Mem Size = %lu\n", (unsigned long)global_mem_size);
    		printf("Local Mem Size = %lu\n", (unsigned long)local_mem_size);
    	}
    	printf("\n");
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

int main(int argc, char **argv)
{
	selectDevice();
}