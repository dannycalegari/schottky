UNAME = $(shell uname -s)
ifeq ($(UNAME),Darwin)
  CC=clang++
else
  CC=g++
endif
CC_C=cc
CFLAGS=-O3 #-g -Wall -Wextra -pedantic
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11

# schottky      the interactive GUI (as before)
# certify_arc   headless trap / interval-certificate driver   (no X11 at all)
# funddom       limit traps on a fundamental domain of C^*/<sigma^b>  (plain C)
#
# certify_arc and funddom need neither the X11 headers nor libX11: the trap code
# is compiled a second time with -DIFS_NO_GRAPHICS, which drops the debug
# windows (see ifs.h).  So on a headless machine
#     make certify_arc funddom
# works with no X11 installed, while plain `make` still builds all three for
# anyone who has it.
all: schottky certify_arc funddom
headless: certify_arc funddom

graphics.o: graphics.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graphics.cc

# figure export (PNG/EPS/PDF).  Uses no X11 and no zlib, so the one object
# serves both the interactive program and the headless tools.
figure_export.o: figure_export.cc figure_export.h
	$(CC) $(CFLAGS) -c figure_export.cc

schottky.o: schottky.cc ifs.cc ifs_gui.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c schottky.cc

trap_grid.o: trap_grid.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trap_grid.cc

movie.o: movie.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c movie.cc

ifs_gui.o: ifs_gui.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c ifs_gui.cc

ifs.o: ifs.cc ifs_draw.cc ifs_trap.cc ifs_interface.cc ifs_connected.cc ifs_trap_like.cc ifs_set_A.cc ifs_set_B.cc ifs_nifs.cc ifs_gifs.cc ifs_2d.cc ifs_picture.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c ifs.cc

schottky: schottky.o graphics.o ifs.o trap_grid.o movie.o ifs_gui.o figure_export.o funddom_core.o
	$(CC) $(CFLAGS) -o schottky schottky.o graphics.o trap_grid.o movie.o ifs.o ifs_gui.o figure_export.o funddom_core.o $(LFLAGS) -lm

# the same three translation units again, without the debug-drawing code, so
# that the headless tools link without libX11 and compile without X11/Xlib.h
NOGFX=-DIFS_NO_GRAPHICS

ifs_nogfx.o: ifs.cc ifs_trap.cc ifs_connected.cc ifs_trap_like.cc ifs_set_A.cc ifs_set_B.cc ifs_nifs.cc ifs_gifs.cc ifs_2d.cc ifs_picture.cc
	$(CC) $(CFLAGS) $(NOGFX) -c ifs.cc -o ifs_nogfx.o

trap_grid_nogfx.o: trap_grid.cc
	$(CC) $(CFLAGS) $(NOGFX) -c trap_grid.cc -o trap_grid_nogfx.o

movie_nogfx.o: movie.cc
	$(CC) $(CFLAGS) $(NOGFX) -c movie.cc -o movie_nogfx.o

certify_arc: certify_arc.cc rigor.h ifs_nogfx.o trap_grid_nogfx.o movie_nogfx.o figure_export.o
	$(CC) $(CFLAGS) $(NOGFX) -o certify_arc certify_arc.cc ifs_nogfx.o trap_grid_nogfx.o movie_nogfx.o figure_export.o -lm

funddom: funddom.c funddom_core.h
	$(CC_C) $(CFLAGS) -DFUNDDOM_MAIN -o funddom funddom.c -lm

# The same file again without -DFUNDDOM_MAIN: no usage(), no main(), just the C API
# of funddom_core.h, so the GUI can call the limit-trap kernel instead of carrying a
# second copy of the mathematics.  Built with $(CC_C) because funddom.c is C99
# (it uses <complex.h>, which is not a C++ type -- hence the POD-only header).
funddom_core.o: funddom.c funddom_core.h
	$(CC_C) $(CFLAGS) -c funddom.c -o funddom_core.o

# Self test for figure_export.  The compressor is hand-rolled, so the check that
# matters is against a real, independent implementation: the C++ half writes the
# compressed corpus and Python's zlib round-trips it.  Needs no X11.
figure_export_test: figure_export_test.cc figure_export.o figure_export.h
	$(CC) $(CFLAGS) -o figure_export_test figure_export_test.cc figure_export.o -lm

test: figure_export_test
	@mkdir -p test_out
	./figure_export_test test_out
	@python3 -c 'import glob, sys, zlib;\
	bad = [f for f in sorted(glob.glob("test_out/defl_*.bin"))\
	       if zlib.decompress(open(f,"rb").read()) != open(f[:-4]+".orig","rb").read()];\
	print("  %-54s %s" % ("compressor round-trips through zlib", "FAIL "+str(bad) if bad else "ok"));\
	sys.exit(1 if bad else 0)'
	@echo "figures are in test_out/ ; view them, or check the vector ones with"
	@echo "  gs -dNOPAUSE -dBATCH -sDEVICE=nullpage -q test_out/fig.pdf test_out/fig.eps"

# Checks that the funddom C API (funddom_core.o) agrees with the funddom CLI
# pixel for pixel -- the invariant that keeps one copy of the mathematics honest.
# Generates the trap-like balls and the reference raster first, so it is slower
# than `make test`; hence its own target.
# rigor.h is the one file whose whole value is soundness, so it gets its own
# known-answer test: hand-computed answers and the certificate results in C++, plus an
# exact-rational containment check of every interval operation in Python.  See the header
# comment of rigor_test.cc for why the verifier is a separate program in another language.
rigortest: rigor_test.cc rigor.h ifs_nogfx.o trap_grid_nogfx.o movie_nogfx.o figure_export.o
	$(CC) $(CFLAGS) $(NOGFX) -o rigor_test rigor_test.cc ifs_nogfx.o trap_grid_nogfx.o movie_nogfx.o figure_export.o -lm
	@mkdir -p test_out
	./rigor_test test_out/rigor_corpus.txt
	python3 rigor_test.py test_out/rigor_corpus.txt

coretest: funddom_core_test.cc funddom_core.o funddom_core.h funddom certify_arc
	$(CC) $(CFLAGS) -o funddom_core_test funddom_core_test.cc funddom_core.o -lm
	@mkdir -p test_out
	./certify_arc dumptlbmany 45 20 1e-9 40 8 6 > test_out/tlb45.txt
	python3 prune_tlb.py test_out/tlb45.txt test_out/tlb45p.txt
	./funddom s0 ann 16 120 120 0.08 0 test_out/ref.bin 1 test_out/tlb45p.txt
	./funddom_core_test test_out/tlb45p.txt test_out/ref.bin 120

clean:
	rm -f *.o
	rm -f schottky certify_arc funddom figure_export_test funddom_core_test rigor_test
	rm -rf test_out
