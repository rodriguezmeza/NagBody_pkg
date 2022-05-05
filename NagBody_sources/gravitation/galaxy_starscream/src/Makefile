# ----- MAKE FILE FOR starscream CODE -----
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

EXEC = galaxy_starscream

OBJS	= starscream.o starscream_init.o starscream_io.o \
	main.o startrun.o \
	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o

INCL	= starscream.h Makefile \
	globaldefs.h cmdline_defs.h protodefs.h \
	$(H_PATH)general/getparam.h \
	$(H_PATH)general/switchs.h $(H_PATH)general/machines.h \
	$(H_PATH)general/stdinc.h $(H_PATH)general/precision.h
	

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -o $(EXEC)


$(OBJS): $(INCL)

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

