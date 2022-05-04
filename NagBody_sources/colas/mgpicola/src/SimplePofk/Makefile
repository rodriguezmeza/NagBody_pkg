SHELL := /bin/bash

# Set compiler
#CC = g++-mp-6
#CC = g++-mp-7
CC = g++

# Options
OPTIONS = -D_SUBTRACTSHOTNOISE -D_CIC

# Path to FFTW library
#I = -I$(HOME)/local/include
#L = -L$(HOME)/local/lib -lfftw3 -lfftw3_threads
I = -I$(HOME)/NagBody_pkg/local/fftw3/include
L = -L$(HOME)/NagBody_pkg/local/fftw3/lib -lfftw3 -lfftw3_threads
#
#I = -I/opt/local/include
#I += -I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/
#I += -I/opt/local/lib/gcc9/gcc/x86_64-apple-darwin18/9.2.0/include
#L = -L/opt/local/lib -lfftw3 -lfftw3_threads
#L += -L/opt/local/lib/gcc9
#I += -I/usr/local/include
#L += -L/usr/local/lib/gcc/11

C = -O3 -Wall -fopenmp $(OPTIONS)

TARGETS := simplepofk
all: $(TARGETS)

OBJS = main.cpp
main.o: io_ascii.h  io_gadget.h io_ramses.h Makefile

simplepofk: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o
