#!/bin/sh

CC=icc
CFLAGS=-std=c99 -restrict -O3 -ip -ip_no_inlining -vec_report3 -opt_report_filedgemm.opt -opt_report_phaseall -fopenmp

all: dgemm

dgemm: dgemm.c timer.c
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf dgemm dgemm.opt cachegrind.out.*
