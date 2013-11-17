#!/bin/sh

CC=icc
CFLAGS=-std=c99 -O3
COMPFL=-openmp

all: quickser quickpar

quickser: quicksort.c timer.c
	$(CC) $(CFLAGS) -o $@ $^

quickpar: quicksort_parallel.c timer.c
	$(CC) $(CFLAGS) $(COMPFL) -o $@ $^

clean: 
	rm -rf quickser quickpar *.o
		  
