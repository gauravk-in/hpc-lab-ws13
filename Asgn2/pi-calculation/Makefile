CC=gcc
CCFLAGS=-O3
LDFLAGS=-fopenmp

all: serial_pi omp_critical_pi omp_reduction_pi

timer.o: timer.c timer.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -c $^

serial_pi: serial_pi.c timer.c
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

omp_critical_pi: omp_critical_pi.c timer.c
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

omp_reduction_pi: omp_reduction_pi.c timer.c
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

clean: 
	rm -rf serial_pi omp_critical_pi omp_reduction_pi timer.h.gch *.o
