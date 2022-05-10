# ----- MAKE FILE FOR analysis_md CODE -----
# Copyright (c) 2006-2018  Mario A. Rodriguez-Meza, Ciudad de Mexico.
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

EXEC   = analysis_md

OBJS   = main.o analysis_md_io.o plplots.o\
	 startrun.o analysis.o analysis_anim.o \
	 $(H_PATH)nagbody.o $(H_PATH)nagbody_tree.o $(H_PATH)nagbody_io.o \
	 $(H_PATH)clib.o $(H_PATH)getparam.o \
	 $(H_PATH)inout.o $(H_PATH)mathfns.o \
	 $(H_PATH)mathutil.o \
	 $(H_PATH)numrec.o $(H_PATH)lic.o

INCL   = cmdline_defs.h \
	 protodefs.h switches.h globaldefs.h \
	 $(H_PATH)nagbody.h $(H_PATH)nagbody_struct.h $(H_PATH)nagbody_proto.h \
	 $(H_PATH)pldefs.h \
	 $(H_PATH)constant.h $(H_PATH)physconstants.h \
	 $(H_PATH)getparam.h \
	 $(H_PATH)inout.h $(H_PATH)machines.h \
	 $(H_PATH)mathfns.h $(H_PATH)mathutil.h \
	 $(H_PATH)numrec.h $(H_PATH)param.h \
	 $(H_PATH)precision.h $(H_PATH)stdinc.h \
	 $(H_PATH)strings.h \
	 $(H_PATH)vectdefs.h $(H_PATH)vectmath.h \
	 $(H_PATH)times.h $(H_PATH)switchs.h $(H_PATH)lic.h

$(EXEC): $(OBJS) 
	$(CC) $(OPTIONS) $(OBJS) $(LDFLAGS) -o $(EXEC) $(RPATHCMD) `$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)` $(MLIBS)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

