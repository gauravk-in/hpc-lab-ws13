CC=icc

all: matrix_mul_benchmark matrix_mul_seq_benchmark 

timer.o: timer.c timer.h
	$(CC) $^ -c

ocl_features: ocl_features.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)

matrix_mul_benchmark: main_benchmark.c timer.c
	$(CC) $^ -o $@ $(OPENCL_SHLIB)	

matrix_mul_seq_benchmark: main_seq_benchmark.c timer.o
	$(CC) $^ -o $@ 

clean:
	rm -rf *.o matrix_mul_benchmark matrix_mul_seq_benchmark ocl_features
