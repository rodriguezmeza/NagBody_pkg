# ----- MAKE FILE FOR md_blj_n2_sird CODE -----
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

EXEC = md_blj_n2_sird

OBJS	= main.o md_blj_io.o timestep.o startrun.o models.o \
	 bljforcecalc.o \
	 $(H_PATH)NagBody/nagbody.o $(H_PATH)NagBody/nagbody_tree.o $(H_PATH)NagBody/nagbody_io.o \
	 $(H_PATH)general/clib.o $(H_PATH)general/getparam.o $(H_PATH)math/mathfns.o $(H_PATH)math/numrec.o \
	 $(H_PATH)io/inout.o $(H_PATH)math/mathutil.o

INCL	= globaldefs.h cmdline_defs.h protodefs.h \
	 $(H_PATH)NagBody/nagbody.h $(H_PATH)NagBody/nagbody_struct.h $(H_PATH)NagBody/nagbody_proto.h \
	 $(H_PATH)general/getparam.h $(H_PATH)general/machines.h $(H_PATH)io/inout.h $(H_PATH)general/precision.h \
	 $(H_PATH)math/mathfns.h $(H_PATH)math/mathutil.h $(H_PATH)general/stdinc.h $(H_PATH)general/strings.h \
	 $(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h $(H_PATH)math/numrec.h \
	 $(H_PATH)general/constant.h $(H_PATH)physics/physconstants.h $(H_PATH)general/switchs.h 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(OPTIONS) $(LDFLAGS) -o $(EXEC)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

