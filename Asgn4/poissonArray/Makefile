#!/bin/sh

CC=mpic++
CCC=icc
CCFLAGS=-O3

all: poisson poisson_mpi lapiae

poisson: poisson.cpp
	$(CCC) -o $@ $^

poisson_mpi: poisson_mpi.cpp
	$(CC) $(CCFLAGS) -o $@ $^

lapiae: iae_poisson.cpp
	$(CCC) -o $@ $^

clean:
	rm -rf *.o poisson poisson_mpi lapiae *.gnuplot
