# ----- MAKE FILE FOR gbsph CODE -----
# Copyright (c) 2006-2011  M.A. Rodriguez-Meza, Mexico, D.F.
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
#OPT7    =  -DBMAX            
OPT8    =  -I $(HOME)/NagBody_pkg/local/gsl/include
#OPT9    =  -DOPTIONNAME
#OPT10   =  -DOPTIONNAME

OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)\
           $(OPT7) $(OPT8) $(OPT9) $(OPT10) $(OPT11) $(OPT12) $(OPT13)\
	   $(OPT14) $(OPT15)

CFLAGS = $(OPTIONS)
#LDFLAGS = $(LIBS)

H_PATH = ./general_libs/
#H_PATH = ../../../General_libs/

EXEC   = gbsph

OBJS   = main.o forcecalc.o gbsphio.o treeload.o \
	 startrun.o timestep.o \
	 $(H_PATH)general/clib.o $(H_PATH)general/getparam.o \
	 $(H_PATH)io/inout.o $(H_PATH)math/mathfns.o \
	 $(H_PATH)math/mathutil.o \
	 $(H_PATH)math/numrec.o $(H_PATH)general/lic.o 

INCL   = cmdline_defs.h \
	 protodefs.h switches.h globaldefs.h \
	 $(H_PATH)general/constant.h \
	 $(H_PATH)general/getparam.h \
	 $(H_PATH)io/inout.h $(H_PATH)general/machines.h \
	 $(H_PATH)math/mathfns.h $(H_PATH)math/mathutil.h \
	 $(H_PATH)math/numrec.h $(H_PATH)general/param.h \
	 $(H_PATH)general/precision.h $(H_PATH)general/stdinc.h \
	 $(H_PATH)general/strings.h \
	 $(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h \
	 $(H_PATH)general/times.h $(H_PATH)general/switchs.h $(H_PATH)general/lic.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

