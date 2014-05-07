CC=g++
CFLAGS=-g -Wall -Wextra -pedantic
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11
all: schottky

graphics.o: graphics.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graphics.cc

schottky.o: schottky.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c schottky.cc

trap_grid.o: trap_grid.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trap_grid.cc
	
movie.o: movie.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c movie.cc

ifs_gui.o: ifs_gui.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c ifs_gui.cc

ifs.o: ifs.cc ifs_draw.cc ifs_trap.cc ifs_interface.cc ifs_connected.cc ifs_trap_like.cc ifs_set_A.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c ifs.cc

schottky: schottky.o graphics.o ifs.o trap_grid.o movie.o ifs_gui.o
	$(CC) $(CFLAGS) -o schottky schottky.o graphics.o trap_grid.o movie.o ifs.o ifs_gui.o $(LFLAGS) -lm

clean:
	rm *.o
	rm schottky
