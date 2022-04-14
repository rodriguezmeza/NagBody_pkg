# ----- MAKE FILE FOR nbody_n2 CODE -----
# Copyright (c) 2006-2014  M.A. Rodriguez-Meza, Mexico, D.F.
#
# Se define el directorio en donde estan las
# definiciones de las arquitecturas particulares.
MACHINES_DIR = ../../../../Machines/
include $(MACHINES_DIR)machine.inc

# Se definen opciones locales de este codigo
# 	OPT1 -- OPT10
# La serie de opciones OPT11 -- OPT20 y OPTIMIZE son para la arquitectura 
# particular y estan en el archivo machine.inc o pmachine.inc
#
#OPT1    =  -DOPTIONNAME
#OPT2    =  -DOPTIONNAME
#OPT3    =  -DOPTIONNAME
#OPT4    =  -DOPTIONNAME
#OPT5    =  -DOPTIONNAME
#OPT6    =  -DOPTIONNAME

#OPT7    =  -I $(HOME)/plplot/include/plplot/

OPT8    =  -I $(HOME)/NagBody_pkg/local/gsl/include
#OPT9    =  -DOPTIONNAME
#OPT10   =  -DOPTIONNAME

OPTIONS =  $(OPTIMIZE) \
	   $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6) $(OPT7) $(OPT8) \
	   $(OPT9) $(OPT10) $(OPT11) $(OPT12) $(OPT13) $(OPT14) $(OPT15) \
           $(OPT16) $(OPT17) $(OPT18) $(OPT19) $(OPT20)

CCFLAGS = $(OPTIONS)
LDFLAGS = $(LIBS)

#PRECISION = DOUBLEPREC

#H_PATH = ../../../General_libs/
H_PATH = ./general_libs/

EXEC = nbody_n2

OBJS	= main.o nbody_n2_io.o direct_gravcalc.o timestep.o startrun.o \
	$(H_PATH)clib.o $(H_PATH)getparam.o \
	$(H_PATH)mathfns.o $(H_PATH)numrec.o \
	$(H_PATH)inout.o

#	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o \
#	$(H_PATH)math/mathfns.o $(H_PATH)math/numrec.o \
#	$(H_PATH)io/inout.o

INCL	= globaldefs.h data_struc_defs.h cmdline_defs.h protodefs.h \
	$(H_PATH)getparam.h $(H_PATH)machines.h \
	$(H_PATH)mathfns.h $(H_PATH)stdinc.h \
	$(H_PATH)strings.h $(H_PATH)precision.h \
	$(H_PATH)vectdefs.h $(H_PATH)vectmath.h \
	$(H_PATH)numrec.h \
	$(H_PATH)inout.h \
	$(H_PATH)switchs.h

#	$(H_PATH)general/getparam.h $(H_PATH)general/machines.h \
#	$(H_PATH)math/mathfns.h $(H_PATH)general/stdinc.h \
#	$(H_PATH)general/strings.h $(H_PATH)general/precision.h \
#	$(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h \
#	$(H_PATH)math/numrec.h \
#	$(H_PATH)io/inout.h \
#	$(H_PATH)general/switchs.h


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

