CC=g++
CFLAGS=-g -Wall -O3 -fast -fopenmp
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11
all: schottky

schottky: schottky.cc points.cc graphics.cc ifs.cc connected.cc interface.cc
	$(CC) $(CFLAGS) $(IFLAGS) -o schottky schottky.cc $(LFLAGS) -lm

clean:
	rm schottky
