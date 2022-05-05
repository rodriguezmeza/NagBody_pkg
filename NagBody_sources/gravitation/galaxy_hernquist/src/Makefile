# ----- MAKE FILE FOR galaxy_hernquist CODE -----
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

EXEC = galaxy_hernquist

####################################################


OBJS   = write_snap.o hernquist.o expdisk.o besselfunc.o spline.o qromb.o \
	trapzd.o polint.o nrutil.o \
	main.o startrun.o \
	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o


INCL   = spline.h nrutil.h xrandom.h Makefile \
	globaldefs.h data_struc_defs.h cmdline_defs.h protodefs.h \
	$(H_PATH)general/getparam.h \
	$(H_PATH)general/switchs.h $(H_PATH)general/machines.h \
	$(H_PATH)general/stdinc.h $(H_PATH)general/precision.h


CFLAGS = $(OPTIONS) $(GSL_CFLAGS)

LIBS   =   -g  $(GSL_LIBS) $(OPENMPFLAG) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL)


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)


