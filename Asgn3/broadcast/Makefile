#!/bin/sh

CC=mpicc
CCFLAGS=-O3

all: trivial_broadcast tree_broadcast broadcast

timer.o: timer.c
	$(CC) $(CCFLAGS) -c $^

broadcast: broadcast.c timer.o
	$(CC) $(CCFLAGS) -o $@ $^

trivial_broadcast: trivial_broadcast.c timer.o
	$(CC) $(CCFLAGS) -o $@ $^

tree_broadcast: tree_broadcast.c timer.o
	$(CC) $(CCFLAGS) -o $@ $^ -lm

clean:
	rm -rf trivial_broadcast tree_broadcast broadcast *.o