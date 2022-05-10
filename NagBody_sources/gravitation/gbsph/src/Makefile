# ----- MAKE FILE FOR gbsph CODE -----
# Copyright (c) 2006-2011  M.A. Rodriguez-Meza, Mexico, D.F.
#
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

EXEC   = gbsph

OBJS   = main.o forcecalc.o gbsphio.o treeload.o \
	 startrun.o timestep.o \
	 $(H_PATH)general/clib.o $(H_PATH)general/getparam.o \
	 $(H_PATH)io/inout.o $(H_PATH)math/mathfns.o \
	 $(H_PATH)math/mathutil.o \
	 $(H_PATH)math/numrec.o $(H_PATH)general/lic.o 

INCL   = cmdline_defs.h \
	 protodefs.h switches.h globaldefs.h \
	 $(H_PATH)general/constant.h \
	 $(H_PATH)general/getparam.h \
	 $(H_PATH)io/inout.h $(H_PATH)general/machines.h \
	 $(H_PATH)math/mathfns.h $(H_PATH)math/mathutil.h \
	 $(H_PATH)math/numrec.h $(H_PATH)general/param.h \
	 $(H_PATH)general/precision.h $(H_PATH)general/stdinc.h \
	 $(H_PATH)general/strings.h \
	 $(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h \
	 $(H_PATH)general/times.h $(H_PATH)general/switchs.h $(H_PATH)general/lic.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC) $(MLIBS)

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

