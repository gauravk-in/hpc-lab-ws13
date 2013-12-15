#define BLOCKSIZE 32

__kernel void g_product_operator_parallel (	__global const float* grid,
						__global float* result,
						unsigned long grid_points_1d)
{
	float mesh_width = 1.0/((float)(grid_points_1d-1));
	
	int bx = get_group_id(0);
	int by = get_group_id(1);

	int tx = get_local_id(0);
	int ty = get_local_id(1);

	int i = by * BLOCKSIZE + ty + 1;
	int j = bx * BLOCKSIZE + tx + 1;

//	result[(i*grid_points_1d)+j] = 1.0;

	result[(i*grid_points_1d)+j] =  (
					(4.0*grid[(i*grid_points_1d)+j]) 
					- grid[((i+1)*grid_points_1d)+j]
					- grid[((i-1)*grid_points_1d)+j]
					- grid[(i*grid_points_1d)+j+1]
					- grid[(i*grid_points_1d)+j-1]
					) * (mesh_width*mesh_width);
	barrier(CLK_LOCAL_MEM_FENCE);
}
