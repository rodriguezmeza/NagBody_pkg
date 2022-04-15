# ----- MAKE FILE FOR nbody_n2 CODE -----
# Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.
#
# Define directory where particular machines are set
MACHINES_DIR = ../../../../Machines/
include $(MACHINES_DIR)machine.inc

# Set local options to this code
# 	OPT1 -- OPT10
# OPTIMIZE is to set some particular arquitecture.
# It is set in machine.inc o pmachine.inc
#
OPT1    =  -I ./
OPT2    =  -I ./general_libs/general
OPT3    =  -I ./general_libs/math
OPT4    =  -I ./general_libs/io
#OPT5    =  -DOPTIONNAME
#OPT6    =  -DOPTIONNAME
#OPT7    =  -I $(HOME)/plplot/include/plplot/
#OPT8    =  -I $(HOME)/NagBody_pkg/local/gsl/include
#OPT9    =  -DOPTIONNAME
#OPT10   =  -DOPTIONNAME

INCLUDES = $(OPT1) $(OPT2) $(OPT3) $(OPT4)

OPTIONS =  $(OPTIMIZE) \
	   $(OPT5) $(OPT6) $(OPT7) $(OPT8) \
	   $(OPT9) $(OPT10)

CPPFLAGS = $(INCLUDES)
CCFLAGS = $(OPTIONS)
LDFLAGS = $(LIBS)

#PRECISION = DOUBLEPREC

H_PATH = ./general_libs/

EXEC = nbody_n2

OBJS	= main.o nbody_n2_io.o direct_gravcalc.o timestep.o startrun.o \
	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o \
	$(H_PATH)math/mathfns.o $(H_PATH)math/numrec.o \
	$(H_PATH)io/inout.o

INCL	= globaldefs.h data_struc_defs.h cmdline_defs.h protodefs.h \
	$(H_PATH)general/getparam.h $(H_PATH)general/machines.h \
	$(H_PATH)math/mathfns.h $(H_PATH)general/stdinc.h \
	$(H_PATH)general/strings.h $(H_PATH)general/precision.h \
	$(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h \
	$(H_PATH)math/numrec.h \
	$(H_PATH)io/inout.h \
	$(H_PATH)general/switchs.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

