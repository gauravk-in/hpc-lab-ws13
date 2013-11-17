#!/bin/sh

CC=icc
CFLAGS=-std=c99 -restrict -O3 -ip -ip_no_inlining -vec_report3
COMPFL=-openmp

quickser: quicksort.c timer.c
		  $(CC) $(CFLAGS) -o $@ $^

quickpar1: quicksort_parallel.c timer.c
		  $(CC) $(CFLAGS) $(COMPFL) -o $@ $^
		  
quickpar2: quicksort_parallel2.c timer.c
		  $(CC) $(CFLAGS) $(COMPFL) -o $@ $^  
		  
quickpar3: quicksort_parallel3.c timer.c
		  $(CC) $(CFLAGS) $(COMPFL) -o $@ $^

quickpar4: quicksort_parallel4.c timer.c
		  $(CC) $(CFLAGS) $(COMPFL) -o $@ $^
		  
