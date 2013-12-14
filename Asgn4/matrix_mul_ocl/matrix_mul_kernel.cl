#define BLOCK_SIZE 32
#define GLOBAL_WORK_SIZE (Width/BLOCK_SIZE)

__kernel void matrixMultiplyParallel ( 	__global const float* M,
										__local float* M_l,
										__global const float* N,
										__local float* N_l,
										__global float* P,
										int Width )
{
	int k;
	
	int bx = get_group_id(0);
	int by = get_group_id(1);

	int tx = get_local_id(0);
	int ty = get_local_id(1);

	int cur_block = 0;

	float product = 0;

	M_l[tx + ty * BLOCK_SIZE] = M[(by * GLOBAL_WORK_SIZE + ty) * Width + (bx * GLOBAL_WORK_SIZE + tx)];
	N_l[tx + ty * BLOCK_SIZE] = M[(by * GLOBAL_WORK_SIZE + ty) * Width + (bx * GLOBAL_WORK_SIZE + tx)];

	barrier(CLK_LOCAL_MEM_FENCE);

	#pragma unroll
	for(k = 0; k < BLOCK_SIZE; k++)
	{
		product += M_l[k + ty * BLOCK_SIZE] * N_l[tx + k * BLOCK_SIZE];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	P[(bx * BLOCK_SIZE + tx) + (by * BLOCK_SIZE + ty) * Width] = product;
}
