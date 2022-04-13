# ----- MAKE FILE FOR flmoon CODE -----
#
# Se define el directorio en donde estan las
# definiciones de las arquitecturas particulares.
#MACHINES_DIR = ../../../../Machines/
#include $(MACHINES_DIR)machine.inc

# Se definen opciones locales de este codigo
# 	OPT1 -- OPT10
# La serie de opciones OPT11 -- OPT20 y OPTIMIZE son para la arquitectura 
# particular y estan en el archivo machine.inc o pmachine.inc
#
#OPT1	=  -DOPTIONNAME     
#OPT2	=  -DOPTIONNAME     
#OPT3	=  -DOPTIONNAME     
#OPT4	=  -DOPTIONNAME     
#OPT5	=  -DOPTIONNAME     
#OPT6	=  -DOPTIONNAME    
#OPT7	=  -DOPTIONNAME
#OPT8	=  -DOPTIONNAME
#OPT9	=  -DOPTIONNAME
#OPT10	=  -DOPTIONNAME

OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)\
	   $(OPT7) $(OPT8) $(OPT11) $(OPT10) \
	   $(OPT11) $(OPT12) $(OPT13) $(OPT14) $(OPT15) $(OPT16) $(OPT17) \
	   $(OPT18) $(OPT19) $(OPT20)

CFLAGS = $(OPTIONS)
LDFLAGS = $(LIBS)


EXEC = xsplint

OBJS	= xsplint.o spline.o splint.o \
	nrutil_nransi.o

INCL	= nr_nransi.h nrutil_nransi.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

