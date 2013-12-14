CC=icc

all: ocl_features matrix_mul writeInputMatrix matrix_mul_benchmark matrix_mul_seq_benchmark

ocl_features: ocl_features.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)

matrix_mul: main.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)

matrix_mul_benchmark: main_benchmark.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)	

matrix_mul_seq_benchmark: main_seq_benchmark.c
	$(CC) $^ -o $@ 

writeInputMatrix: writeInputMatrix.c
	$(CC) $^ -o $@

run: matrix_mul
	./matrix_mul M_matrix.txt N_matrix.txt P_matrix.txt 10

clean:
	rm -rf *.o matrix_mul writeInputMatrix matrix_mul_benchmark matrix_mul_seq_benchmark ocl_features
