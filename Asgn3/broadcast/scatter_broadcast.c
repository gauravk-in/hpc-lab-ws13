#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include "timer.h"
#include <math.h>

#define N 200000

int main(int argc, char **argv)
{
	int rank, size;
	double *array, *recv;
	int i;
	long length, bytes_sent;
	MPI_Status status;
	time_marker_t time;
	double time_taken;
	int receiver;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	//printf("Rank = %d\n", rank);

	//Initialize array and insert padding if necessary
	if(N % size  > 0) length = N + (N % size);
	else length = N;
	
	array = (double*) malloc (sizeof(double) * length);
	
	recv = (double*) malloc (sizeof(double) * (length/size));
	
	
	if(rank == 0)
	{
		printf("Size = %d\n", size);
		printf("Length = %ld\n", length);
		printf("ChunkLength = %ld\n", length/size);
		for (i=0; i<N; i++)
		{
			array[i] = i;
		}
		for (i=N; i<length; i++)
		{
			array[i] = 0;
		}
		time = get_time();
	}
	
	//broadcast operation
	if(rank == 0)
	{
		MPI_Scatter(array, length/size, MPI_DOUBLE_PRECISION, recv, length/size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Allgather(recv, length/size, MPI_DOUBLE_PRECISION, array, length/size, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);
	
	
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0)
	{
		time_taken = get_ToD_diff_time (time);
		printf("Thread(0) : Time taken = %e\n", time_taken);
		bytes_sent = sizeof(double) * (length / size) * (size - 1)  + sizeof(double) * length *  (size - 1) ;
		printf("Thread(0) : Bytes sent = %ld\n", bytes_sent);
		printf("Thread(0) : Bandwidth = %e\n", bytes_sent/time_taken);
		
	}
	
	MPI_Finalize();
	return 0;
}
