#!/bin/sh

CC=icc
CFLAGS=-std=c99 -restrict -O3 -ip -ip_no_inlining -vec_report3 -opt_report_filevector.opt -opt_report_phaseall

all: vector dgemm

vector: vector.c
	$(CC) $(CFLAGS) -o $@ $^

dgemm: dgemm.c timer.c
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf vector vector.opt dgemm dgemm.opt cachegrind.out.*
