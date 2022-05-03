# ----- MAKE FILE FOR md_blj_n2 CODE -----
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
NEMOPATH = $(H_PATH)nemo/

EXEC = md_blj_n2

OBJS	= main.o md_blj_io.o timestep.o startrun.o models.o \
	 bljforcecalc.o \
	 $(H_PATH)nagbody.o $(H_PATH)nagbody_tree.o $(H_PATH)nagbody_io.o \
	 $(H_PATH)clib.o $(H_PATH)getparam.o $(H_PATH)mathfns.o $(H_PATH)numrec.o \
	 $(H_PATH)inout.o $(H_PATH)mathutil.o $(H_PATH)lic.o

INCL	= globaldefs.h cmdline_defs.h protodefs.h \
	 $(H_PATH)nagbody.h $(H_PATH)nagbody_struct.h $(H_PATH)nagbody_proto.h \
	 $(H_PATH)getparam.h $(H_PATH)machines.h $(H_PATH)inout.h $(H_PATH)precision.h \
	 $(H_PATH)mathfns.h $(H_PATH)mathutil.h $(H_PATH)stdinc.h $(H_PATH)strings.h \
	 $(H_PATH)vectdefs.h $(H_PATH)vectmath.h $(H_PATH)numrec.h $(H_PATH)lic.h \
	 $(H_PATH)constant.h $(H_PATH)physconstants.h $(H_PATH)switchs.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(OPTIONS) $(LDFLAGS) -o $(EXEC)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

