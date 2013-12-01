#!/bin/sh

CC=mpic++
CCFLAGS=-O3

all: poisson poisson_mpi

poisson: poisson.cpp
	$(CC) $(CCFLAGS) -o $@ $^

poisson_mpi: poisson_mpi.cpp
	$(CC) $(CCFLAGS) -o $@ $^

clean:
	rm -rf *.o poisson poisson_mpi *.gnuplot