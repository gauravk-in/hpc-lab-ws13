/** 
 * matrix matrix multiplication pattern for practical course
 **/

#include "timer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

//#define BLOCKSIZE 128 

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
	int BLOCKSIZE;
	int num_threads;
	float * a;
        float * b;
        float * c;

        long mem_size;

        int bk, bj, i, j, k;
        int b_base, c_base;
        float r;
	FILE *fp;

        /*char logfile_name[100];
        FILE *logfile_handle;*/

	num_threads = omp_get_max_threads();

        n = 1024;
        if(argc > 1){
                n = atoi(argv[1]);
        }
        if(argc > 2){
                fp = fopen(argv[2], "w");
        }

	if(num_threads > 1)
		BLOCKSIZE = min(256, n/num_threads);
	else
		BLOCKSIZE = min(1024, n);

	printf("Block Size = %d\n", BLOCKSIZE);        
	//sprintf(logfile_name, "logfile_dgemm.txt");
        //logfile_handle = freopen(logfile_name, "w", stdout);

        mem_size = n * n * sizeof(float);
        a = (float*)malloc(mem_size);
        b = (float*)malloc(mem_size);
        c = (float*)malloc(mem_size);
        if(0 == a || 0 == b || 0 == c){
                printf("memory allocation failed");
                return 0;
        }

        /* initialisation */
        for (i = 0; i < n; i++){
                for (j = 0; j < n; j++){
                        *(a + i * n + j) = (float)i + (float)j;
                        *(b + i * n + j) = (float)(n - i) + (float)(n - j);
                }
        }
        memset(c, 0, mem_size);

        time_marker_t time = get_time();
        float flops;

	__m256 a1_r;
	__m256 a2_r;
	__m256 a3_r;
	__m256 a4_r;
	__m256 b1_r;
	__m256 b2_r;
	__m256 b3_r;
	__m256 b4_r;
	__m256 c1_r;
	__m256 c2_r;
	__m256 c3_r;
	__m256 c4_r;


        for(bk = 0; bk < n; bk+=BLOCKSIZE) {
		#pragma omp parallel for schedule(static) private(a1_r, a2_r, a3_r, a4_r, b1_r, b2_r, b3_r, b4_r, c1_r, c2_r, c3_r, c4_r, i, k, j, bj)
                for(bj = 0; bj < n; bj+=BLOCKSIZE) {
                        for(i = 0; i < n; i+=2){
//				#pragma omp parallel for schedule(static) private(a1_r, a2_r, a3_r, a4_r, b1_r, b2_r, b3_r, b4_r, c1_r, c2_r, c3_r, c4_r)
                                for(k = bk; k < min(n, bk+BLOCKSIZE); k+=2){
                                                        a1_r = _mm256_set_ps(a[i * n + k], a[i * n + k], a[i * n + k], a[i * n + k], a[i * n + k], a[i * n + k], a[i * n + k], a[i * n + k]);
                                                        a2_r = _mm256_set_ps(a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1], a[i * n + k + 1]);
                                                        a3_r = _mm256_set_ps(a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k], a[(i+1) * n + k]);
                                                        a4_r = _mm256_set_ps(a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1], a[(i+1) * n + k+1]);
							for(j = bj ; j < min(n, bj + BLOCKSIZE); j+=16){
                                                                c1_r = _mm256_load_ps(&c[i * n + j]);
                                                                c2_r = _mm256_load_ps(&c[i * n + j + 8]);
                                                                c3_r = _mm256_load_ps(&c[(i+1) * n + j]);
                                                                c4_r = _mm256_load_ps(&c[(i+1) * n + j + 8]);
                                                                b1_r = _mm256_load_ps(&b[k * n + j]);
                                                                b2_r = _mm256_load_ps(&b[k * n + j + 8]);
                                                                b3_r = _mm256_load_ps(&b[(k+1) * n + j]);
                                                                b4_r = _mm256_load_ps(&b[(k+1) * n + j + 8]);
                                                                c1_r = _mm256_add_ps(c1_r, _mm256_add_ps(_mm256_mul_ps(a1_r, b1_r), _mm256_mul_ps(a2_r, b3_r)));
                                                                c2_r = _mm256_add_ps(c2_r, _mm256_add_ps(_mm256_mul_ps(a1_r, b2_r), _mm256_mul_ps(a2_r, b4_r)));
                                                                c3_r = _mm256_add_ps(c3_r, _mm256_add_ps(_mm256_mul_ps(a3_r, b1_r), _mm256_mul_ps(a4_r, b3_r)));
                                                                c4_r = _mm256_add_ps(c4_r, _mm256_add_ps(_mm256_mul_ps(a3_r, b2_r), _mm256_mul_ps(a4_r, b4_r)));
                                                                _mm256_store_ps(&c[i * n + j], c1_r);
                                                                _mm256_store_ps(&c[i * n + j + 4], c2_r);
                                                                _mm256_store_ps(&c[(i+1) * n + j], c3_r);
                                                                _mm256_store_ps(&c[(i+1) * n + j + 4], c4_r);
                                                                //c[i * n + j] += r * b[k * n + j];
                                                        }

                                }
                        }
                }
        }

        flops = 2.0 * n * n * n / 1000000000;

        print_flops(flops, time);

        if(argc > 2)
        {
                for(i = 0; i < n; i++)
                        for(j = 0; j < n; j++)
                                fprintf(fp, "%f ", c[i * n + j]);
        }

        return(0);
}
