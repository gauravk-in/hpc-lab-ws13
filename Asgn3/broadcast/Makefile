#!/bin/sh

CC=mpicc
CCFLAGS=-O3

all: trivial_broadcast

trivial_broadcast: trivial_broadcast.c
	$(CC) $(CCFLAGS) -o $@ $^

clean:
	rm -rf trivial_broadcast *.o