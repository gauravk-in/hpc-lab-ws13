__kernel void reduce_parallel(__global const float* input,
								__global float* reduce_tmp)
{
	int y = get_group_id(0) * BLOCKSIZE + get_local_id(0) + 1;
	float i;
	float two = 2;
	float log_width = ceil(log2((float)grid_points_1d-two));

	reduce_tmp[y] = input[y];

    for(i=log_width; i>0; i-=1)
    {
        if(y <= pown(two, i) && y >= pown(two, i-1))
        {
            reduce_tmp[y % (unsigned long)pown(two, i-1) + 1] += reduce_tmp[y];
        }
    }

}