# ----- MAKE FILE FOR analysis_galaxy CODE -----
# Copyright (c) 2006-2022  Mario A. Rodriguez-Meza, Ciudad de Mexico.
#
# Machine definitions
MACHINES_DIR = ./
include $(MACHINES_DIR)Makefile_machine

# OPTIMIZE is to set some particular architecture. See MACHINE section in Makefile_machine.

OPTIONS =  $(OPTIMIZE) $(OPT2)

CCFLAGS += $(OPTIONS) $(OMPFLAG)
CPPFLAGS = $(CCFLAGS) $(INCLUDES)
LDFLAGS += $(OPTIMIZE) $(OMPFLAG) $(LIBS) 

EXEC   = analysis_galaxy

H_PATH = ./general_libs/
NEMOPATH = $(H_PATH)nemo/

OBJS   = main.o analysis_galaxy_io.o plplots.o\
	 startrun.o analysis.o analysis_anim.o \
	 $(H_PATH)NagBody/nagbody_gravforcecalc.o $(H_PATH)NagBody/nagbody_gravtreeload.o \
	 $(H_PATH)NagBody/nagbody.o $(H_PATH)NagBody/nagbody_tree.o $(H_PATH)NagBody/nagbody_io.o \
	 $(H_PATH)general/clib.o $(H_PATH)general/getparam.o \
	 $(H_PATH)io/inout.o $(H_PATH)math/mathfns.o \
	 $(H_PATH)math/mathutil.o \
	 $(H_PATH)math/numrec.o $(H_PATH)general/lic.o \
	 $(NEMOPATH)src/kernel/io/filesecret.o $(NEMOPATH)src/kernel/io/extstring.o \
	 $(NEMOPATH)src/kernel/io/history.o $(NEMOPATH)src/kernel/cores/bswap.o \
	 $(NEMOPATH)src/kernel/io/convert.o $(NEMOPATH)src/kernel/io/dprintf.o \
	 $(NEMOPATH)src/kernel/cores/error.o $(NEMOPATH)src/kernel/cores/strlib.o \
	 $(NEMOPATH)src/kernel/io/stropen.o $(NEMOPATH)src/kernel/io/getparam.o \
	 $(NEMOPATH)src/kernel/misc/nemoinp.o $(NEMOPATH)src/kernel/cores/date_id.o \
	 $(NEMOPATH)src/kernel/cores/file_size.o $(NEMOPATH)src/kernel/io/outdefv.o \
	 $(NEMOPATH)src/kernel/cores/filefn.o $(NEMOPATH)src/kernel/cores/burststring.o \
	 $(NEMOPATH)src/kernel/misc/herinp.o $(NEMOPATH)src/kernel/misc/setfblank.o

INCL   = cmdline_defs.h \
	 protodefs.h globaldefs.h \
	 $(H_PATH)NagBody/nagbody.h $(H_PATH)NagBody/nagbody_struct.h $(H_PATH)NagBody/nagbody_proto.h \
	 $(H_PATH)visual/pldefs.h \
	 $(H_PATH)general/constant.h \
	 $(H_PATH)general/getparam.h \
	 $(H_PATH)io/inout.h $(H_PATH)general/machines.h \
	 $(H_PATH)math/mathfns.h $(H_PATH)math/mathutil.h \
	 $(H_PATH)math/numrec.h $(H_PATH)general/param.h \
	 $(H_PATH)general/precision.h $(H_PATH)general/stdinc.h \
	 $(H_PATH)general/strings.h \
	 $(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h \
	 $(H_PATH)general/times.h $(H_PATH)general/switchs.h $(H_PATH)general/lic.h \
	 $(NEMOPATH)include/stdinc.h $(NEMOPATH)include/strlib.h \
	 $(NEMOPATH)include/filestruct.h $(NEMOPATH)include/history.h \
	 $(NEMOPATH)include/extstring.h $(NEMOPATH)include/getparam.h \
	 $(NEMOPATH)src/kernel/io/filesecret.h $(NEMOPATH)src/tools/unproto/stdarg.h \
	 $(NEMOPATH)src/kernel/misc/gipsyc.h \
	 $(NEMOPATH)include/snapshot/snapshot.h $(NEMOPATH)include/snapshot/body.h \
	 $(NEMOPATH)include/snapshot/put_snap.c $(NEMOPATH)include/snapshot/put_snap-old.c \
	 $(NEMOPATH)include/version.h $(NEMOPATH)include/options.h $(NEMOPATH)include/mathfns.h \
	 $(NEMOPATH)include/filefn.h $(NEMOPATH)lib/maxsizes.h $(NEMOPATH)lib/config.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -o $(EXEC) $(RPATHCMD) `$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)`

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

