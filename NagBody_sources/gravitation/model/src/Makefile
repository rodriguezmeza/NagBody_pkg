# ----- MAKE FILE FOR model CODE -----
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

EXEC   = model

H_PATH = ./general_libs

OBJS   = main.o models_io.o models.o \
	startrun.o galaxy.o \
	$(H_PATH)/NagBody/nagbody.o \
	$(H_PATH)/NagBody/nagbody_tree.o \
	$(H_PATH)/NagBody/nagbody_io.o \
	$(H_PATH)/general/clib.o \
	$(H_PATH)/general/getparam.o \
	$(H_PATH)/io/inout.o \
	$(H_PATH)/math/mathfns.o \
	$(H_PATH)/math/mathutil.o \
	$(H_PATH)/math/numrec.o

INCL   = cmdline_defs.h \
	protodefs.h switches.h globaldefs.h models.h \
	$(H_PATH)/NagBody/nagbody.h \
	$(H_PATH)/NagBody/nagbody_struct.h \
	$(H_PATH)/NagBody/nagbody_proto.h \
	$(H_PATH)/general/constant.h \
	$(H_PATH)/general/getparam.h \
	$(H_PATH)/general/machines.h \
	$(H_PATH)/general/param.h \
	$(H_PATH)/general/precision.h \
	$(H_PATH)/general/stdinc.h \
	$(H_PATH)/general/strings.h \
	$(H_PATH)/general/switchs.h \
	$(H_PATH)/io/inout.h \
	$(H_PATH)/physics/physconstants.h \
	$(H_PATH)/math/mathfns.h \
	$(H_PATH)/math/mathutil.h \
	$(H_PATH)/math/numrec.h \
	$(H_PATH)/math/vectdefs.h \
	$(H_PATH)/math/vectmath.h

CFLAGS =   $(OPTIONS)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIONS) $(OBJS) $(LIBS)   -o  $(EXEC) $(MLIBS)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)
