/** 
 * Quicksort implementation for practical course
 **/

#include "timer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
	//#pragma omp parallel section
	
		#pragma omp task untied final(length < 1000)
		quicksort(data, right);
		#pragma omp task untied final(length < 1000)
		quicksort(&(data[left]), length - left);
	
}

int check(double *data, int length){
	int i;
	for(i = 1; i < length; i++)
		if(data[i] < data[i-1]) return 1;
	return 0;
}

int main(int argc, char **argv)
{
	int length;
	double *data;

	int mem_size;

	int i, j, k;

	length = 10000000;
	if(argc > 1){
		length = atoi(argv[1]);
	}

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

	/*print_list(data, length);
	if(check(data, length) != 0)
		printf("ERROR\n");*/

	printf("Size of dataset: %d, elapsed time[s] %e \n", length, get_ToD_diff_time(time));

	return(0);
}

