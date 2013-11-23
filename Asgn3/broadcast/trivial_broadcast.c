#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

#define N 10

int main(int argc, char **argv)
{
	int rank, size;
	double *array;
	int i;
	MPI_Status status;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	printf("Rank = %d\n", rank);

	array = (double*) malloc (sizeof(double) * N);

	if(rank == 0)
	{
		for (i=0; i<N; i++)
		{
			array[i] = i;
		}

		for (i=1; i<size; i++)
		{
			MPI_Send(array, N, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(array, N, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, &status);
		printf("Received by Thread = %d\n", rank);
	}
	
	MPI_Finalize();
	return 0;
}