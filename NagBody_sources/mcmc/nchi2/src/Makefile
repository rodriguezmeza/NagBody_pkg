# ----- MAKE FILE FOR nchi2 CODE -----
# Copyright (c) 2006-2020  Mario A. Rodriguez-Meza, Mexico, D.F.
#
# Se define el directorio en donde estan las
# definiciones de las arquitecturas particulares.
MACHINES_DIR = ../../../../Machines/
include $(MACHINES_DIR)machine.inc

MPINAGBODY = 1
# For MPI parallel version comment line below.
MPINAGBODY = 0


# Se definen opciones locales de este codigo
# 	OPT1 -- OPT10
# La serie de opciones OPT11 -- OPT20 y OPTIMIZE son para la arquitectura 
# particular y estan en el archivo machine.inc o pmachine.inc
#

OPTIMIZE =

#
# LINES TO not GNU version and Mac OS X version (Catalina 10.15.3 at least)
# Can be commented
#
#OPT1    =  -DNOGNU
# TO SOLVE LINKING ERROR:
# ld: library not found for -lSystem
# Fit OPT9 and OPT10 to your needs:
OPT9    =  -I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include
OPT10   =  -L/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib
#
#

#OPT2    =  -DOPTIONNAME
#OPT3    =  -DOPTIONNAME
#OPT4    =  -DOPTIONNAME
#OPT5    =  -DOPTIONNAME
#OPT6    =  -DOPTIONNAME
#OPT7    =  -I $(HOME)/plplot/include/plplot/
#OPT8    =  -DOPTIONNAME


OPTIONS =  $(OPTIMIZE) \
	   $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6) $(OPT7) $(OPT8) \
	   $(OPT9) $(OPT10) $(OPT11) $(OPT12) $(OPT13) $(OPT14) $(OPT15) \
           $(OPT16) $(OPT17) $(OPT18) $(OPT19) $(OPT20)


EXEC = nchi2

all: $(EXEC)

ifeq ($(MPINAGBODY),1)
CC = mpicc
MPICODE = -DMPICODE
OPTIONS += $(MPICODE)
else
# C compiler:
CC = gcc
# With icc (Intel):
#CC = icc
endif

# Default CFLAGS:
CFLAGS = -g -O3 $(OPTIONS)
#
# With OpenMP:
#CFLAGS = -g -O3 -fopenmp $(OPTIONS)
# For icc use instead:
#CFLAGS = -g -O3 -qopenmp $(OPTIONS)

LDFLAGS = $(LIBS)

#PRECISION = DOUBLEPREC

#
# Nothing to do below
#
#H_PATH = ../../../General_libs/
#
# This line is for GITHUB VERSION
H_PATH = general_libs/

OBJS	= main.o nchi2_io.o nchi2.o startrun.o \
	models.o nminchi2.o \
	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o \
	$(H_PATH)math/mathfns.o $(H_PATH)math/numrec.o \
	$(H_PATH)math/minpack.o $(H_PATH)math/linpack.o \
	$(H_PATH)io/inout.o $(H_PATH)math/mathutil.o $(H_PATH)math/diffeqs.o \
	$(H_PATH)math/quads.o

INCL	= globaldefs.h data_struc_defs.h cmdline_defs.h protodefs.h \
	models.h constants_defs.h model_user.h \
	models_io.h model_user_io.h \
	models/model_baocpl.h models/model_sncpl.h models/model_gseir_ci.h \
	models/model_rcbecho.h models/model_rcpiso.h \
	$(H_PATH)general/getparam.h $(H_PATH)general/machines.h \
	$(H_PATH)math/mathfns.h $(H_PATH)general/stdinc.h \
	$(H_PATH)math/mathfns.h $(H_PATH)general/stdinc.h \
	$(H_PATH)general/strings.h $(H_PATH)general/precision.h \
	$(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h \
	$(H_PATH)math/numrec.h $(H_PATH)math/minpack.h \
	$(H_PATH)io/inout.h $(H_PATH)math/mathutil.h \
	$(H_PATH)general/switchs.h $(H_PATH)math/diffeqs.h \
	$(H_PATH)math/quads.h $(H_PATH)math/linpack.h

$(EXEC): $(OBJS) 
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(CFLAGS) -o $@ -lm

# LINES FOR sMGPT GITHUB VERSION:
#$(EXEC): $(OBJS) 
#	($(CC) $(OPTIONS) $(OBJS) $(LIBS) $(CFLAGS) -o $@ -lm; cp $(EXEC) ../)

$(OBJS): $(INCL)

clean:
	(rm -f $(OBJS) $(EXEC); rm -fR nchi2.dSYM)


# LINES FOR sMGPT GITHUB VERSION:
#clean:
#	(rm -f $(OBJS) $(EXEC); rm -fR nchi2.dSYM; rm ../$(EXEC))


.PHONY : all clean check install uninstall

