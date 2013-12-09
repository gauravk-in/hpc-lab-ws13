#define BLOCK_SIZE 32

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

	for(cur_block = 0; cur_block < Width/BLOCK_SIZE; cur_block++)
	{
		M_l[tx + ty * BLOCK_SIZE] = M[(cur_block * BLOCK_SIZE + tx) + (by * BLOCK_SIZE + ty) * Width];
		N_l[tx + ty * BLOCK_SIZE] = N[(bx * BLOCK_SIZE + tx) + (cur_block * BLOCK_SIZE + ty) * Width];

		barrier(CLK_LOCAL_MEM_FENCE);

		#pragma unroll
		for(k = 0; k < BLOCK_SIZE; k++)
		{
			product += M_l[k + ty * BLOCK_SIZE] * N_l[tx + k * BLOCK_SIZE];
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	}

	P[(bx * BLOCK_SIZE + tx) + (by * BLOCK_SIZE + ty) * Width] = product;
}
