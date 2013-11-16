/** 
 * matrix matrix multiplication pattern for practical course
 **/

#include "timer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>


#define BLOCKSIZE 1024

int max(int a, int b){
        if(a >= b) return a;
        else return b;
}

inline int min(int a, int b){
        if(a <= b) return a;
        else return b;
}

//extern __m256d _mm256_loaddup_pd(double *);

int main(int argc, char **argv)
{
        int n;
        double * restrict a;
        double * restrict b;
        double * restrict c;

        long mem_size;

        int bk, bj, i, j, k, l;
        int b_base, c_base;
        double r;
	FILE *fp;

        /*char logfile_name[100];
        FILE *logfile_handle;*/

        n = 1024;
        if(argc > 1){
                n = atoi(argv[1]);
        }
	if(argc > 2){
		fp = fopen(argv[2], "w");
	}	

        //sprintf(logfile_name, "logfile_dgemm.txt");
        //logfile_handle = freopen(logfile_name, "w", stdout);

        mem_size = n * n * sizeof(double);
        a = (double*)malloc(mem_size);
        b = (double*)malloc(mem_size);
        c = (double*)malloc(mem_size);
        if(0 == a || 0 == b || 0 == c){
                printf("memory allocation failed");
                return 0;
        }

        /* initialisation */
        for (i = 0; i < n; i++){
                for (j = 0; j < n; j++){
                        *(a + i * n + j) = (double)i + (double)j;
                        *(b + i * n + j) = (double)(n - i) + (double)(n - j);
                }
        }
        memset(c, 0, mem_size);

        time_marker_t time = get_time();
        double flops;

        for(bk = 0; bk < n; bk+=BLOCKSIZE) {
                for(bj = 0; bj < n; bj+=BLOCKSIZE) {
                        for(i = 0; i < n; i++){
                                for(k = bk, l = bk; k < min(n, bk+BLOCKSIZE); k++, l+=4){
                                        //r = a[i * n + k];
                                        __m256d a_r = _mm256_load_pd(&a[i * n + l]);
                                        for(j = bj ; j < min(n, bj + BLOCKSIZE); j+=4){
                                                __m256d c_r = _mm256_load_pd(&c[i * n + j]);
                                                __m256d b_r = _mm256_load_pd(&b[k * n + j]);
                                                c_r = _mm256_add_pd(c_r, _mm256_mul_pd(a_r, b_r));
						_mm256_store_pd(&c[i * n + j], c_r);
                                                //c[i * n + j] += r * b[k * n + j];
                                        }
                                }
                        }
                }
        }

        flops = 2.0 * n * n * n / 1000000;

        print_flops(flops, time);

	if(argc > 2)
	{
		for(i = 0; i < n; i++)
			for(j = 0; j < n; j++)
				fprintf(fp, "%f ", c[i * n + j]);
	}

        return(0);
}
