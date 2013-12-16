#!/bin/sh

CC=icpc
CCFLAGS=-O3

all: poisson poisson_ocl trial

poisson: poisson.cpp
	$(CC) $(CCFLAGS) -o $@ $^

poisson_ocl: poisson_ocl.cpp
	$(CC) $(CCFLAGS) -o $@ $^ $(OPENCL_SHLIB)

trial: trial.cpp
	$(CC) $(CCFLAGS) -o $@ $^ $(OPENCL_SHLIB)

clean:
	rm -rf *.o poisson poisson_ocl *.gnuplot trial
