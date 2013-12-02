#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include "timer.h"
#include <math.h>

#define N 200000

int main(int argc, char **argv)
{
	int rank, size;
	double *array;
	int i;
	MPI_Status status;
	time_marker_t time;
	double time_taken;
	long bytes_sent;
	int receiver;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	// printf("Rank = %d\n", rank);

	array = (double*) malloc (sizeof(double) * N);

	if(rank == 0)
	{
		printf("Size = %d\n", size);
		for (i=0; i<N; i++)
		{
			array[i] = i;
		}
		time = get_time();
	}
	

	for(i=0; i<ceil(log(size)/log(2)); i++)
	{
		if(rank < pow(2,i))
		{
			if(pow(2,i)+rank < size)
			{
				MPI_Send(array, N, MPI_DOUBLE_PRECISION, pow(2,i)+rank, 1, MPI_COMM_WORLD);
				// printf("Thread %d sent to %d\n", rank, (int)pow(2,i)+rank);
			}
		}
		if(rank >= pow(2,i) && rank <pow(2,i+1) && rank < size)
		{
			MPI_Recv(array, N, MPI_DOUBLE_PRECISION, rank-pow(2,i), 1, MPI_COMM_WORLD, &status);
			//printf("Recvd\n");
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0)
	{
		time_taken = get_ToD_diff_time (time);
		printf("Thread(0) : Time taken = %e\n", time_taken);
		bytes_sent = sizeof(double) * N * (size - 1);
		printf("Thread(0) : Bytes sent = %ld\n", bytes_sent);
		printf("Thread(0) : Bandwidth = %e\n", bytes_sent/time_taken);
	}
	
	MPI_Finalize();
	return 0;
}
