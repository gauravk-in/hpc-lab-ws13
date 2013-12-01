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

scatter_broadcast: scatter_broadcast.c timer.o
	$(CC) $(CCFLAGS) -o $@ $^

bonus_broadcast: bonus_broadcast.c timer.o
	$(CC) $(CCFLAGS) -o $@ $^

clean:
	rm -rf trivial_broadcast tree_broadcast broadcast scatter_broadcast bonus_broadcast *.o
