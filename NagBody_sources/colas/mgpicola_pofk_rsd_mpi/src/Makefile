SHELL := /bin/bash

# Set compiler
#CC = mpicc-openmpi-gcc6 -std=c99
CC = mpicc -std=c99

# Deconvolve the windowfunction?
OPTIONS = -D_DECONVOLVEWINDOWFUNCTION -D_READVELOCITY

I = -I$(HOME)/NagBody_pkg/local/fftw3/include
L = -L$(HOME)/NagBody_pkg/local/fftw3/lib -lfftw3_mpi -lfftw3
#I = -I/opt/local/include
#L = -L/opt/local/lib -lfftw3_mpi -lfftw3


C = -O3 -Wall $(OPTIONS)

TARGETS := mgpicola_pofk_rsd_mpi
all: $(TARGETS)

# OBJECT FILES
OBJS = rsd_main.c

# DEPENDENCIES
main.o		: rsd_main.c global.h io_gadget.h Makefile

# HEADERS
HEADERS = io_gadget.h

mgpicola_pofk_rsd_mpi: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML) -lm

%.o: %.cpp $(HEADERS)
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o

