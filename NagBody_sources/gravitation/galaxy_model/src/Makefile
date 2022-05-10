# ----- MAKE FILE FOR galaxy_model CODE -----
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

EXEC = galaxy_model

OBJS	= main.o galmodio.o startrun.o galmod_general.o \
	galmod_bulge.o galmod_disk.o galmod_halo.o \
	galmod_bulgevel.o galmod_halovel.o galmod_satmod.o galmod_diskvel.o \
	$(H_PATH)NagBody/nagbody_gravforcecalc.o $(H_PATH)NagBody/nagbody_gravtreeload.o \
	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o \
	$(H_PATH)math/mathfns.o $(H_PATH)math/numrec.o $(H_PATH)io/inout.o

INCL	= globaldefs.h cmdlinedefs.h protodefs.h \
	$(H_PATH)NagBody/nagbody.h $(H_PATH)NagBody/nagbody_struct.h $(H_PATH)NagBody/nagbody_proto.h \
	$(H_PATH)general/stdinc.h $(H_PATH)math/vectdefs.h $(H_PATH)general/getparam.h \
	$(H_PATH)general/machines.h $(H_PATH)math/mathfns.h $(H_PATH)math/numrec.h \
	$(H_PATH)general/precision.h $(H_PATH)io/inout.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -o $(EXEC) $(MLIBS)

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

