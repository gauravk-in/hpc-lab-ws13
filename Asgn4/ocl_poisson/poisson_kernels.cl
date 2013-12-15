#define BLOCKSIZE 32

__kernel void g_product_operator_parallel (	__global const float* grid,
						__global float* result,
						unsigned long grid_points_1d)
{
	float mesh_width = 1.0/((float)(grid_points_1d-1));
	
	int i = get_group_id(1) * BLOCKSIZE + get_local_id(1) + 1;
	int j = get_group_id(0) * BLOCKSIZE + get_local_id(0) + 1;

	result[(i*grid_points_1d)+j] =  (
					(4.0*grid[(i*grid_points_1d)+j]) 
					- grid[((i+1)*grid_points_1d)+j]
					- grid[((i-1)*grid_points_1d)+j]
					- grid[(i*grid_points_1d)+j+1]
					- grid[(i*grid_points_1d)+j-1]
					) * (mesh_width*mesh_width);
}

__kernel void g_dot_product_parallel (	__global const float* grid1,
						__global const float* grid2,
						__global float* reduce_tmp,
						unsigned long grid_points_1d)
{
	int x = 0;
	int y = get_group_id(0) * BLOCKSIZE + get_local_id(0) + 1;
	float tmp = 0.0;
	float i;
	float two = 2;
	float log_width = ceil(log2((float)grid_points_1d-two));

	for(x = 1; x < grid_points_1d-1; x++)
		reduce_tmp[y] += grid1[y * grid_points_1d + x] * grid2[y * grid_points_1d + x];

	for(i=log_width; i>0; i-=1)
	{
		if(y <= pown(two, i))
		{
			reduce_tmp[y % (unsigned long)pown(two, i-1)] += reduce_tmp[y];
		}
	}
}
