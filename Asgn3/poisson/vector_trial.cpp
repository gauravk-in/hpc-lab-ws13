#include <stdio.h>
#include <mpi.h>

#define N 4

int main(int argc, char **argv) {

    int myrank, nprocs;
    int i,j;

    int matrix[N][N]={{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}};
    int col[N];
    int matrix1[N][N] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    printf("Nprocs = %d\n", nprocs);

    // Define type "column"
    MPI_Datatype column;
    MPI_Type_vector(N,1,N,MPI_INT,&column);
    MPI_Type_commit(&column);

    if(myrank==0){
        j=1;
        MPI_Send(&matrix[0][j],1,column,1,99,MPI_COMM_WORLD);
    }

    if(myrank==1){
        // **** FIRST MODE: Don't use "column" type *****
        // MPI_Recv(col,N,MPI_INT,0,99,MPI_COMM_WORLD,&info);

        // **** SECOND MODE: Use "column" type *****
        MPI_Recv(matrix1,1,column,0,99,MPI_COMM_WORLD,&status);

        // printf("\nColumn: ");
        // for(j=0;j<N;j++)
        //     printf("\n %d",col[j]);

        printf("\nMatrix: ");
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                printf("| %2d |", matrix1[i][j]);
            printf("\n");
    }

    MPI_Type_free(&column);
    MPI_Finalize();

    return 0;
}