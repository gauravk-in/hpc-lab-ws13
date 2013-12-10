#!/bin/sh

CC=icc
CFLAGS=-std=c99 -O3
COMPFL=-openmp

all: quickser quickpar quickcilk quickcilk2

quickser: quicksort.c timer.c
	$(CC) $(CFLAGS) -o $@ $^

quickpar: par_quicksort.c timer.c
	$(CC) $(CFLAGS) $(COMPFL) -o $@ $^

quickcilk: cilk_quicksort.c timer.c
	$(CC) $(CFLAGS) -o $@ $^

quickcilk2: cilk2_quicksort.c timer.c
	$(CC) $(CFLAGS) -o $@ $^

clean: 
	rm -rf quickser quickpar quickcilk quickcilk2 *.o
		  
