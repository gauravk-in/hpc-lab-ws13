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

int main(int argc, char **argv)
{
        int n;
        double * restrict a;
        double * restrict b;
        double * restrict c;

        long mem_size;

        int bk, bj, i, j, k;
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

	__m256d a1_r;
	__m256d a2_r;
	__m256d a3_r;
	__m256d a4_r;
	__m256d b1_r;
	__m256d b2_r;
	__m256d b3_r;
	__m256d b4_r;
	__m256d c1_r;
	__m256d c2_r;
	__m256d c3_r;
	__m256d c4_r;


        for(bk = 0; bk < n; bk+=BLOCKSIZE) {
		#pragma omp parallel for schedule(static) private(a1_r, a2_r, a3_r, a4_r, b1_r, b2_r, b3_r, b4_r, c1_r, c2_r, c3_r, c4_r, i, k, j, bj)
                for(bj = 0; bj < n; bj+=BLOCKSIZE) {
                        for(i = 0; i < n; i+=2){
//				#pragma omp parallel for schedule(static) private(a1_r, a2_r, a3_r, a4_r, b1_r, b2_r, b3_r, b4_r, c1_r, c2_r, c3_r, c4_r)
                                for(k = bk; k < min(n, bk+BLOCKSIZE); k+=2){
                                                        a1_r = _mm256_set_pd(a[i * n + k], a[i * n + k], a[i * n + k], a[i * n + k]);
                                                        a2_r = _mm256_set_pd(a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1]);
                                                        a3_r = _mm256_set_pd(a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k]);
                                                        a4_r = _mm256_set_pd(a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1]);
//							#pragma omp parallel for private(a1_r, a2_r, a3_r, a4_r)
							for(j = bj ; j < min(n, bj + BLOCKSIZE); j+=8){
                                                                c1_r = _mm256_load_pd(&c[i * n + j]);
                                                                c2_r = _mm256_load_pd(&c[i * n + j + 4]);
                                                                c3_r = _mm256_load_pd(&c[(i+1) * n + j]);
                                                                c4_r = _mm256_load_pd(&c[(i+1) * n + j + 4]);
                                                                b1_r = _mm256_load_pd(&b[k * n + j]);
                                                                b2_r = _mm256_load_pd(&b[k * n + j + 4]);
                                                                b3_r = _mm256_load_pd(&b[(k+1) * n + j]);
                                                                b4_r = _mm256_load_pd(&b[(k+1) * n + j + 4]);
                                                                c1_r = _mm256_add_pd(c1_r, _mm256_add_pd(_mm256_mul_pd(a1_r, b1_r), _mm256_mul_pd(a2_r, b3_r)));
                                                                c2_r = _mm256_add_pd(c2_r, _mm256_add_pd(_mm256_mul_pd(a1_r, b2_r), _mm256_mul_pd(a2_r, b4_r)));
                                                                c3_r = _mm256_add_pd(c3_r, _mm256_add_pd(_mm256_mul_pd(a3_r, b1_r), _mm256_mul_pd(a4_r, b3_r)));
                                                                c4_r = _mm256_add_pd(c4_r, _mm256_add_pd(_mm256_mul_pd(a3_r, b2_r), _mm256_mul_pd(a4_r, b4_r)));
                                                                _mm256_store_pd(&c[i * n + j], c1_r);
                                                                _mm256_store_pd(&c[i * n + j + 4], c2_r);
                                                                _mm256_store_pd(&c[(i+1) * n + j], c3_r);
                                                                _mm256_store_pd(&c[(i+1) * n + j + 4], c4_r);
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
