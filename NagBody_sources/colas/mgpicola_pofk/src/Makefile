SHELL := /bin/bash

# Set compiler
#CC = g++-mp-6
#CC = g++-mp-7
#CC = g++-mp-9
CC = g++

# Options
OPTIONS = -D_SUBTRACTSHOTNOISE -D_CIC

# Path to FFTW library
I = -I$(HOME)/NagBody_pkg/local/fftw3/include
L = -L$(HOME)/NagBody_pkg/local/fftw3/lib -lfftw3 -lfftw3_threads
#I = -I/opt/local/include
#L = -L/opt/local/lib -lfftw3 -lfftw3_threads

C = -O3 -Wall -fopenmp $(OPTIONS)

TARGETS := mgpicola_pofk
all: $(TARGETS)

OBJS = main.cpp
main.o: io_ascii.h  io_gadget.h io_ramses.h Makefile

mgpicola_pofk: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o
