#include <stdio.h>
#include <stdlib.h>

void readMatrix(char *fileName, float *mat, int width)
{
	FILE *fp = fopen(fileName, "r");
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			fscanf(fp, "%f ", &mat[i * width + j]);
		}
	}

	fclose(fp);
	return;
}

void writeInputMatrix(char *fileName, float value, int width)
{
	FILE *fp = fopen(fileName, "w");
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			fprintf(fp, "%f ", value);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

void writeMatrix(char *fileName, float *mat, int width)
{
	FILE *fp = fopen(fileName, "w");
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			fprintf(fp, "%f ", mat[i * width + j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

void printMatrix(float *mat, int width)
{
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			printf("%f ", mat[i * width + j]);
		}
		printf("\n");
	}
	return;
}

void fillMatrix(float *mat, int width, float value)
{
	int i, j;

	for(i = 0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			mat[i*width+j] = value;
		}
	}
}

int main(int argc, char **argv)
{
	size_t size;
    float *M_host;
    float *N_host;
    float *P_host;
    char *M_fileName;
    char *N_fileName;
    char *P_fileName;
    int width;
    struct timeval tv_start_kernel;
    struct timeval tv_end_kernel;
    unsigned long us_kernel;
    int x, y, k;
    
    if(argc != 5)
    {
    	printf("Usage : %s <M_matrix.txt> <N_matrix.txt> <P_matrix.txt> <width>\n", argv[0]);
    	exit(1);
    }

    sscanf(argv[4], "%d", &width);

    size = width * width * sizeof(float);
	M_host = (float *) malloc(size);
	N_host = (float *) malloc(size);
	P_host = (float *) malloc(size);

	M_fileName = argv[1];
	N_fileName = argv[2];
	P_fileName = argv[3];

    readMatrix(M_fileName, M_host, width);
    readMatrix(N_fileName, N_host, width);

    gettimeofday(&tv_start_kernel, NULL);

    for(x = 0; x < width; x++)
	{
		for(y = 0; y < width; y++)
		{
			P_host[x + y * width] = 0;
			//#pragma omp parallel for
			for(k = 0; k < width; k++)
			{
				P_host[x + y * width] += M_host[k + y * width] * N_host[x + k * width];
			}
		}
	}

   	gettimeofday(&tv_end_kernel, NULL);

    us_kernel = (tv_end_kernel.tv_sec - tv_start_kernel.tv_sec) * 1000000;
    us_kernel += (tv_end_kernel.tv_usec - tv_start_kernel.tv_usec);

    printf("\nTime taken for only kernel = %lu usec\n", us_kernel);

    writeMatrix(P_fileName, P_host, width);

    return 0;
}
