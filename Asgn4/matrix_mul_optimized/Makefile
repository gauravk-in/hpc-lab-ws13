CC=icc

all: matrix_mul writeInputMatrix matrix_mul_benchmark

matrix_mul: main.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)

matrix_mul_benchmark: main_benchmark.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)	

matrix_mul_seq_benchmark: main_seq_benchmark.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)

writeInputMatrix: writeInputMatrix.c
	$(CC) $^ -o $@

run: matrix_mul
	./matrix_mul M_matrix.txt N_matrix.txt P_matrix.txt 10

clean:
	rm -rf *.o matrix_mul writeInputMatrix matrix_mul_benchmark matrix_mul_seq_benchmark
