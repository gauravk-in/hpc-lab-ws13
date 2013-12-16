#!/bin/sh

CC=mpic++
CCC=icc
CCFLAGS=-O3

all: poisson lapiae lapiae2 lapsimd

poisson: poisson.cpp
	$(CCC) -o $@ $^

poisson_mpi: poisson_mpi.cpp
	$(CC) $(CCFLAGS) -o $@ $^

lapiae: iae_poisson.cpp
	$(CCC) -o $@ $^

lapiae2: iae2_poisson.cpp
	$(CCC) -o $@ $^

lapsimd: simd_poisson.cpp
	$(CCC) -O3 -restrict -vec_report3 -o $@ $^

clean:
	rm -rf *.o poisson poisson_mpi lapiae lapiae2 lapsimd *.gnuplot
