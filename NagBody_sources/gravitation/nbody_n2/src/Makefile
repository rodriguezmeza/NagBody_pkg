# ----- MAKE FILE FOR nbody_n2 CODE -----
# Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.
#
# Machine definitions
MACHINES_DIR = ./
include $(MACHINES_DIR)Makefile_machine

# OPTIMIZE is to set some particular architecture. See MACHINE section in Makefile_machine.

OPTIONS =  $(OPTIMIZE) $(OPT2)

CCFLAGS += $(OPTIONS) $(OMPFLAG)
CPPFLAGS = $(CCFLAGS) $(INCLUDES)
LDFLAGS += $(OPTIMIZE) $(OMPFLAG) $(LIBS) 

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
	$(CC) $(OBJS) $(LDFLAGS) -o $(EXEC) $(MLIBS)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

