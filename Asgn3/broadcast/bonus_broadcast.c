#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include "timer.h"
#include <math.h>

#define N 100000

int main(int argc, char **argv)
{
	int rank, size;
	double *array, *recv;
	int i, j;
	long bytes_sent, length;
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
	
	for(i = 1; i < size; i++){
		if(rank == 0)
		{
			MPI_Send((array+(length/size)*i), length/size, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD);
			//printf("Send operation %d completed\n", i);
		}
		else if(rank == i)
		{
			MPI_Recv(recv, length/size, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, &status);
			*(array + rank * (length/size)) = *(recv);
			//printf("Receive operation %d of process %d completed\n", i, rank);
		}
	}
	//printf("CheckPoint 1-%d reached\n", rank);
	for(i = 0; i < size; i++){
		for(j = 1; j < size; j++)
		{
			if(i == 0 && rank == 0)
			{
			 	MPI_Send(array, length/size, MPI_DOUBLE_PRECISION, j, 1, MPI_COMM_WORLD);
				//printf("Send operation %d of process %d to process %d completed\n", i, rank, j);
			}
			else if(rank == i && rank > 0 && i != j)
			{
				MPI_Send(recv, length/size, MPI_DOUBLE_PRECISION, j, 1, MPI_COMM_WORLD);
				//printf("Send operation %d of process %d to process %d completed\n", i, rank, j);
			}
			else if(rank == j && rank != i)
			{
				MPI_Recv((array + (length/size)*i), length/size, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, &status);
				//printf("Receive operation %d of process %d completed\n", i, rank);
				*(array + i * (length/size)) = *(recv);
			}
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
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