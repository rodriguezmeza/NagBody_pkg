# ----- MAKE FILE FOR templateBessel CODE -----
# Copyright (c) 2006-2011  M.A. Rodriguez-Meza, Mexico, D.F.
#
# Se define el directorio en donde estan las
# definiciones de las arquitecturas particulares.
MACHINES_DIR = ../../../../Machines/
include $(MACHINES_DIR)machine.inc

# Se definen opciones locales de este codigo
# 	OPTION1 -- OPTION10
# La serie de opciones OPTION11 -- OPTION20 y OPTIMIZE son para la arquitectura 
# particular y estan en el archivo machine.inc o pmachine.inc
#
#OPTION1    =  -DOPTIONNAME
#OPTION2    =  -DOPTIONNAME
#OPTION3    =  -DOPTIONNAME
#OPTION4    =  -DOPTIONNAME
#OPTION5    =  -DOPTIONNAME
#OPTION6    =  -DNDIM=3

#OPTION7    =  -I $(HOME)/plplot/include/plplot/

#OPTION8    =  -DDOUBLEPREC
#OPTION8    =  -DMIXEDPREC
#OPTION9    =  -DGRAVITY
#OPTION10   =  -DMACOSX

OPTIONS =  $(OPTIMIZE) \
	   $(OPTION1) $(OPTION2) $(OPTION3) $(OPTION4) $(OPTION5) \
	   $(OPTION6) $(OPTION7) $(OPTION8) \
	   $(OPTION9) $(OPTION10) $(OPTION11) $(OPTION12) $(OPTION13) \
	   $(OPTION14) $(OPTION15) \
           $(OPTION16) $(OPTION17) $(OPTION18) $(OPTION19) $(OPTION20)

CCFLAGS = $(OPTIONS)
LDFLAGS = $(LIBS)

H_PATH = ../../../General_libs/

EXEC = templateBessel

OBJS	= main.o besselfunc.o \
	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o 

INCL	= protodefs.h \
	$(H_PATH)general/getparam.h $(H_PATH)general/switchs.h $(H_PATH)general/machines.h \
	$(H_PATH)math/mathfns.h $(H_PATH)general/stdinc.h $(H_PATH)general/precision.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -o $(EXEC)


########################################################################
# Build individual object files.

main.o: main.c protodefs.h $(H_PATH)general/getparam.h \
	$(H_PATH)math/mathfns.h $(H_PATH)general/stdinc.h 
	$(CC) $(CCFLAGS) \
	   -o main.o -c main.c

besselfunc.o: besselfunc.c protodefs.h $(H_PATH)general/getparam.h \
	$(H_PATH)math/mathfns.h $(H_PATH)general/stdinc.h 
	$(CC) $(CCFLAGS) \
	   -o besselfunc.o -c besselfunc.c

$(H_PATH)clib.o: $(H_PATH)general/clib.c $(H_PATH)general/getparam.h \
	$(H_PATH)general/stdinc.h $(H_PATH)general/switchs.h $(H_PATH)general/machines.h 
	$(CC) $(CCFLAGS) \
	   -o $(H_PATH)general/clib.o -c $(H_PATH)general/clib.c

$(H_PATH)getparam.o: $(H_PATH)general/getparam.c $(H_PATH)general/getparam.h \
	$(H_PATH)general/stdinc.h
	$(CC) $(CCFLAGS) \
	   -o $(H_PATH)general/getparam.o -c $(H_PATH)general/getparam.c


########################################################################

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

