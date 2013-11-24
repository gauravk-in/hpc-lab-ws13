#!/bin/sh

CC=mpicc
CCFLAGS=-O3

all: trivial_broadcast tree_broadcast

timer.o: timer.c
	$(CC) $(CCFLAGS) -c $^

trivial_broadcast: trivial_broadcast.c timer.o
	$(CC) $(CCFLAGS) -o $@ $^

tree_broadcast: tree_broadcast.c timer.o
	$(CC) $(CCFLAGS) -o $@ $^ -lm

clean:
	rm -rf trivial_broadcast tree_broadcast *.o