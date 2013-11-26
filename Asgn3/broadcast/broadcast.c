#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include "timer.h"

#define N 1000

int main(int argc, char **argv)
{
	int rank, size;
	double *array;
	int i;
	MPI_Status status;
	time_marker_t time;
	double time_taken;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	// printf("Rank = %d\n", rank);

	array = (double*) malloc (sizeof(double) * N);

	if(rank == 0)
	{
		for (i=0; i<N; i++)
		{
			array[i] = i;
		}

		time = get_time();
	}

	MPI_Bcast (array, N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
	
	if(rank == 0)
	{
		time_taken = get_ToD_diff_time (time);
		printf("Thread(0) : Time taken = %e\n", time_taken);
	}

	MPI_Finalize();
	return 0;
}