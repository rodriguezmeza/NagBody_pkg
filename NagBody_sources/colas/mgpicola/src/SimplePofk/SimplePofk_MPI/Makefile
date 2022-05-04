SHELL := /bin/bash

# Set compiler
#CC = mpicc-openmpi-gcc6
CC = mpicc

# Deconvolve the windowfunction?
OPTIONS = -D_DECONVOLVEWINDOWFUNCTION

I = -I$(HOME)/NagBody_pkg/local/fftw3/include
L = -L$(HOME)/NagBody_pkg/local/fftw3/lib -lm -lfftw3_mpi -lfftw3
#I = -I/opt/local/include
#L = -L/opt/local/lib -lm -lfftw3_mpi -lfftw3

C = -O3 -Wall $(OPTIONS)

TARGETS := simplepofk_mpi
all: $(TARGETS)

OBJS = main.c
main.o		: main.c io_gadget.h Makefile

simplepofk_mpi: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

