#!/bin/sh

CC=icpc
CCFLAGS=-O3

all: poisson poisson_ocl

poisson: poisson.cpp
	$(CC) $(CCFLAGS) -o $@ $^

poisson_ocl: poisson_ocl.cpp
	$(CC) $(CCFLAGS) -o $@ $^ $(OPENCL_SHLIB)

clean:
	rm -rf *.o poisson poisson_ocl *.gnuplot