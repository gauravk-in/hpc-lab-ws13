CC=gcc

all: matrix_mul writeInputMatrix matrix_mul_benchmark

matrix_mul: main.c
	$(CC) $^ -o $@ -I../opencl10 -L/usr/lib/nvidia-304 -lOpenCL

matrix_mul_benchmark: main_benchmark.c
	$(CC) $^ -o $@ -I../opencl10 -L/usr/lib/nvidia-304 -lOpenCL	

writeInputMatrix: writeInputMatrix.c
	$(CC) $^ -o $@

run: matrix_mul
	./matrix_mul M_matrix.txt N_matrix.txt P_matrix.txt 10

clean:
	rm -rf *.o matrix_mul writeInputMatrix matrix_mul_benchmark
