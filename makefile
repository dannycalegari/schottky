CC=g++
CFLAGS=-g -Wall #-O3
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11
all: schottky

graphics.o: graphics.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graphics.cc

schottky.o: schottky.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c schottky.cc

ifs.o: ifs.cc ifs_draw.cc ifs_trap.cc ifs_interface.cc ifs_connected.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c ifs.cc

schottky: schottky.o graphics.o ifs.o
	$(CC) $(CFLAGS) -o schottky schottky.o graphics.o ifs.o $(LFLAGS) -lm

clean:
	rm *.o
	rm schottky
