#!/bin/sh

CC=icc
CFLAGS=-std=c99 -restrict - O3 -ip -ip_no_inlining -vec_report3
COMPFL=-openmp

quickser: quicksort.c timer.c
		  $(CC) $(CFLAGS) -o $@ $^

quickpar1: quicksort_parallel.c timer.c
		  $(CC) $(COMPFL) -o $@ $^
		  
