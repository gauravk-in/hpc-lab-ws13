#!/bin/sh

CC=icc
CFLAGS=-std=c99 -restrict -O3 -ip -ip_no_inlining -vec_report3 -opt_report_filedgemm.opt -opt_report_phaseall -openmp -xAVX

all: dgemm

dgemm: dgemm.c timer.c
	$(CC) $(CFLAGS) -o $@ $^

run_test: dgemm
	./dgemm 1024 output_1024.txt > output_flops_1024.txt
	./dgemm 2048 output_2048.txt > output_flops_2048.txt
	./dgemm 4096 output_4096.txt > output_flops_4096.txt

clean:
	rm -rf dgemm dgemm.opt cachegrind.out.* output_*
