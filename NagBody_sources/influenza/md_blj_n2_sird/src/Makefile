# ----- MAKE FILE FOR md_blj_n2_sird CODE -----
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
OPT10   =  -D_FORTIFY_SOURCE=0
# CON LA DEFINICION ANTERIOR EVITAMOS al correr el codigo
# EL Abort trap: 6 que aparece la rutina “output” cuando llama a strcp(mode,”w!”)

GSL_INCL= -I$(HOME)/NagBody_pkg/local/gsl/include
GSL_LIBS= -L$(HOME)/NagBody_pkg/local/gsl/lib

#
# LINES TO not GNU version and Mac OS X version (Catalina 10.15.3 at least)
#OPT1    =  -DNOGNU
# TO SOLVE LINKING ERROR:
# ld: library not found for -lSystem
OPT9    =  -I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include
OPT10   =  -L/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib

#OPTIMIZE = -O3 -g
# To avoid running error: Illegal instruction: 4
OPTIMIZE = -O2 -mmacosx-version-min=10.0

OPTIONS =  $(OPTIMIZE) \
	   $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6) $(OPT7) $(OPT8) \
	   $(OPT9) $(OPT10) $(OPT11) $(OPT12) $(OPT13) $(OPT14) $(OPT15) \
           $(OPT16) $(OPT17) $(OPT18) $(OPT19) $(OPT20) $(GSL_INCL) $(GSL_LIBS)

CFLAGS = $(OPTIONS)
LDFLAGS = $(LIBS) $(OPT10)

H_PATH = ../../../General_libs/
NEMOPATH = $(H_PATH)nemo/

EXEC = md_blj_n2_sird
#EXEC = md_blj_n2_2d_sird

OBJS	= main.o md_blj_io.o timestep.o startrun.o models.o \
	 bljforcecalc.o \
	 $(H_PATH)NagBody/nagbody.o $(H_PATH)NagBody/nagbody_tree.o NagBody/nagbody_io.o \
	 $(H_PATH)general/clib.o $(H_PATH)general/getparam.o math/mathfns.o math/numrec.o \
	 io/inout.o math/mathutil.o 

#$(H_PATH)general/lic.o 


# 

INCL	= globaldefs.h cmdline_defs.h protodefs.h \
	 NagBody/nagbody.h NagBody/nagbody_struct.h NagBody/nagbody_proto.h \
	 $(H_PATH)general/getparam.h $(H_PATH)general/machines.h io/inout.h $(H_PATH)general/precision.h \
	 math/mathfns.h math/mathutil.h $(H_PATH)general/stdinc.h $(H_PATH)general/strings.h \
	 math/vectdefs.h math/vectmath.h math/numrec.h \
	 $(H_PATH)general/constant.h $(H_PATH)physics/physconstants.h $(H_PATH)general/switchs.h 

#$(H_PATH)general/lic.h

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(OPTIONS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

