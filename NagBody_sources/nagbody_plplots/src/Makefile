# ----- MAKE FILE FOR nplot2d CODE --------------------------
# Copyright (c) 2006-2013  M.A. Rodriguez-Meza, Mexico, D.F.
#
# Se define el directorio en donde estan las
# definiciones de las arquitecturas particulares.
MACHINES_DIR = ../../../../Machines/
include $(MACHINES_DIR)machine.inc

# Se definen opciones locales de este codigo
# 	OPT1 -- OPT10
# La serie de opciones OPT11 -- OPT20 y OPTIMIZE son para la arquitectura 
# particular y estan en el archivo machine.inc o pmachine.inc

#OPT1    =  -DOPTIONNAME     
#OPT2    =  -DOPTIONNAME     
#OPT3    =  -DOPTIONNAME     
#OPT4    =  -DOPTIONNAME     
#OPT5    =  -DOPTIONNAME     
#OPT6    =  -DOPTIONNAME
#OPT7    =  -I $(HOME)/plplot/include/plplot/
# Last working
#OPT7    =  -I $(HOME)/NagBody_pkg/local/plplot/include/plplot/
#OPT7    =  -I $(HOME)/local/plplot/include/plplot/
#OPT7    =  -I /opt/local/lib/plplot510/include/plplot/
#
# Tests with plplot5.15.0 ::
# Working:
OPT7    =  -I /opt/local/include/plplot/
#OPT7    =  -I $(HOME)/NagBody_pkg/local/include/plplot/
#OPT7    =  -I $(HOME)/NagBody_pkg/local/plplot/include/plplot/
#
#
#OPT8    =  `plplot-config --cflagas --libs`

#OPT8    =  -I $(HOME)/NagBody_pkg/local/gsl/include

# GSL is installed in /opt/local
#GSL_INCL= -I$(HOME)/NagBody_pkg/local/gsl/include
#GSL_LIBS= -L$(HOME)/NagBody_pkg/local/gsl/lib

# PARA LA VERSION plplot-5.9.4 SE DEBE USAR:
#LIB_TAG = d
##PKG_CONFIG_ENV = PKG_CONFIG_PATH=$(HOME)/local/plplot-5.9.4/lib/pkgconfig
##PKG_CONFIG_ENV = PKG_CONFIG_PATH=$(HOME)/local/plplot/lib/pkgconfig
## Last working
##PKG_CONFIG_ENV = PKG_CONFIG_PATH=$(HOME)/NagBody_pkg/local/plplot/lib/pkgconfig
##PKG_CONFIG_ENV = PKG_CONFIG_PATH=/opt/local/lib/pkgconfig
##PKG_CONFIG_ENV = PKG_CONFIG_PATH=/usr/lib/pkgconfig
#PKG_CONFIG_ENV = PKG_CONFIG_PATH=/opt/local/lib/plplot510/lib/pkgconfig
##RPATHCMD = -Wl,-rpath -Wl,$(HOME)/local/plplot/lib
##RPATHCMD = -Wl,-rpath -Wl,$(HOME)/NagBody_pkg/local/plplot/lib
#RPATHCMD = 
#

# PARA LA VERSION plplot-5.14.x SE DEBE USAR:
# Tests with plplot5.15.0 ::
# Working:
PKG_CONFIG_ENV = PKG_CONFIG_PATH=/opt/local/lib/pkgconfig
#PKG_CONFIG_ENV = PKG_CONFIG_PATH=$(HOME)/NagBody_pkg/local/lib/pkgconfig
#PKG_CONFIG_ENV = PKG_CONFIG_PATH=$(HOME)/NagBody_pkg/local/plplot/lib/pkgconfig
LIB_TAG =
RPATHCMD = 

OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)\
	   $(OPT7) $(OPT8) $(OPT9) $(OPT10) \
	   $(OPT11) $(OPT12) $(OPT13) $(OPT14) $(OPT15) \
           $(OPT16) $(OPT17) $(OPT18) $(OPT19) $(OPT20) $(GSL_INCL) $(GSL_LIBS)

# para la version :: plplot-5.14.x usar
# Tests with plplot5.15.0 ::
# Working:
#CFLAGS = $(OPTIONS) -I$(HOME)/NagBody_pkg/local/plplot/include/plplot
CFLAGS = $(OPTIONS)
#CFLAGS = $(OPTIONS) -I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include
#LDFLAGS = $(LIBS) -L/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib

H_PATH = ./general_libs/
#H_PATH = ../../../General_libs/
#H_PLPATH = $(HOME)/plplot/include/plplot/
#H_PLPATH = $(HOME)/local/plplot/include/plplot/
#H_PLPATH = $(HOME)/NagBody_pkg/local/plplot/include/plplot/
#H_PLPATH = /opt/local/include/plplot/

EXEC = nplot2d

OBJS	= main.o nplot2dio.o nplot2d.o startrun.o \
	$(H_PATH)general/clib.o $(H_PATH)general/getparam.o $(H_PATH)general/lic.o \
	$(H_PATH)math/mathfns.o $(H_PATH)math/numrec.o $(H_PATH)io/inout.o

INCL	= globaldefs.h cmdlinedefs.h protodefs.h $(H_PATH)visual/pldefs.h \
	 $(H_PATH)general/stdinc.h $(H_PATH)math/vectdefs.h $(H_PATH)general/getparam.h \
	 $(H_PATH)general/lic.h $(H_PATH)general/machines.h $(H_PATH)general/switchs.h \
	 $(H_PATH)general/machines.h $(H_PATH)math/mathfns.h $(H_PATH)math/numrec.h $(H_PATH)math/vectmath.h \
	 $(H_PATH)general/precision.h $(H_PATH)NagBody/nagbody_struct.h \
	 $(H_PATH)io/inout.h $(H_PATH)NagBody/nagbody.h $(H_PATH)NagBody/nagbody_proto.h 


# PARA LA VERSION plplot-5.14.x SE DEBE USAR:
$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) $(LDFLAGS) -o $(EXEC) $(RPATHCMD) `$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)`

# PARA LA VERSION plplot-5.9.4 SE DEBE USAR:
#$(EXEC): $(OBJS) 
#	$(CC) $(OBJS) $(LIBS) $(LDFLAGS) -o $(EXEC) $(RPATHCMD) `$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)`

#$(EXEC): $(OBJS) 
#	$(CC) $(OBJS) -L$(HOME)/NagBody_pkg/local/plplot/lib $(LIBS) -o $(EXEC) $(RPATHCMD) `$(PKG_CONFIG_ENV) pkg-config --cflags --libs plplot$(LIB_TAG)`

# PARA LA VERSION plplot-5.6.1 SE DEBE USAR:
#$(EXEC): $(OBJS) 
#	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
#	$(CC) $(OBJS) $(LIBS) -o $(EXEC) `plplot-config --cflagas --libs`

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

