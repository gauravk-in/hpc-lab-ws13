/** 
 * Quicksort implementation for practical course
 **/

#include "timer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cilk/cilk.h>
#include <math.h>


FILE *fp;

void print_list(double *data, int length){
	int i;
	for(i = 0; i < length; i++)
		printf("%e\t", data[i]);
	printf("\n");
}

void quicksort(double *data, int length){
	if (length <= 1) return;

	//print_list(data, length);

	double pivot = data[0];
	double temp;
	int left = 1;
	int right = length - 1;

	do {
		while(left < (length - 1) && data[left] <= pivot) left++;
		while(right > 0 && data[right] >= pivot) right--;
		
		/* swap elements */
		if(left < right){
			temp = data[left];
			data[left] = data[right];
			data[right] = temp;
		}

	} while(left < right);
	
	if(data[right] < pivot){
		data[0] = data[right];
		data[right] = pivot;
	}

	//print_list(data, length);

	/* recursion */
	cilk_spawn quicksort(data, right);
	quicksort(&(data[left]), length - left);

	cilk_sync;
}

int check(double *data, int length){
	int i;
	for(i = 1; i < length; i++)
		if(data[i] < data[i-1]) return 1;
	return 0;
}


int old_main(int length, int argc)
{
	//int length;
	double *data;

	int mem_size;

	int i, j, k;

	/*
	length = 10000000;
	if(argc > 1){
		length = atoi(argv[1]);
	}
	*/
	data = (double*)malloc(length * sizeof(double));
	if(0 == data){
		printf("memory allocation failed");
		return 0;
	}

	/* initialisation */
	srand(0);
	for (i = 0; i < length; i++){
		data[i] = (double)rand() / (double)RAND_MAX;
	}

	time_marker_t time = get_time();
	
	//print_list(data, length);

	quicksort(data, length);	

	/*print_list(data, length);*/
	if(check(data, length) != 0)
		printf("ERROR\n");

	printf("Size of dataset: %d, elapsed time[s] %e \n", length, get_ToD_diff_time(time));

	if(argc > 2){
		fprintf(fp,"%e;", get_ToD_diff_time(time));
	}

	return(0);
}

int main(int argc, char **argv)
{
        int length, l;
        int nworkers;
	int i, j, temp, temp2, temp3;
	char nthreads[4];

	l = 100000;

        if(argc > 1){
                l = atoi(argv[1]);
	}
	if(argc > 2){
		fp = fopen(argv[2], "a");
	}
	
	if(argc > 2)
	{
		fprintf(fp, ";");
		
		for(i = 1; i < 256; i *= 2)
		{
			fprintf(fp, "%d;", i);
		}

		fprintf(fp, "\n");
	}

	for(j = 1; j < 101; j++)
	{

		length = j * l;

		if(argc > 2)
                {
                	fprintf(fp, "%d;", length);
                }
		
		for(i = 1; i < 256; i *= 2)
		{
			__cilkrts_end_cilk();	
	
			sprintf(nthreads, "%d", i);	

			__cilkrts_set_param("nworkers", nthreads);

			nworkers = __cilkrts_get_nworkers();

			printf("Intel Cilk Plus # workers: %d\n", nworkers);
		

			old_main(length, argc);
					
		}	

		if(argc > 2)
                {
                	fprintf(fp, "\n");
                }
		
	}
	
	if(argc > 2) fclose(fp);

        return 0;
}

