# ----- MAKE FILE FOR gdgt207galaxy CODE -----
# Se define el directorio en donde estan las 
# definiciones de las arquitecturas particulares,
MACHINES_DIR = ./
include $(MACHINES_DIR)makefile_machine

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

#OPT8    =  -DOPTIONNAME
#OPT9    =  -DOPTIONNAME
#OPT10   =  -DOPTIONNAME

#SIMULATIONTYPE="gassphere"
#SIMULATIONTYPE="galaxy"
#SIMULATIONTYPE="cluster"
#SIMULATIONTYPE="lcdm_gas"


# ESTAS OPCIONES PARA TRANSFERIR VARIABLES FUNCIONAN USANDO 
# make
#
# sin opciones É desde un nivel arriba —sea desde NagBody_sources/Gravitation/gadget207/
# usando el archivo de ah’ "NagBody" con "make" solito ...
#

#ifeq ($(simulation),gassphere)
#include makefile_gassphere
#else
#ifeq ($(simulation),galaxy)
include makefile_galaxy
#else
#ifeq ($(simulation),cluster)
#include makefile_cluster
#else
#ifeq ($(simulation),lcdm_gas)
#include makefile_lcdm_gas
#else
#fooprint:
#	@echo " "
#	@echo "No simulation type selected!"
#	@echo "use with the option: simulation=xxx, where xxx is one of"
#	@echo "gassphere, galaxy, cluster, lcdm_gas"
#	@echo " "
#endif
#endif
#endif
#endif


ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(OPT))) 
HDF5INCL =
HDF5LIB  =
endif


OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)\
           $(OPT7) $(OPT8) $(OPT9) $(OPT10) $(OPT11) $(OPT12) $(OPT13)\
	   $(OPT14) $(OPT15) $(OPT16) $(OPT17) $(OPT18) $(OPT19) $(OPT20) $(OPT)

# Otros parametros para el comportamiento del codigo (estan) estaran
# en el archivo switches.h. Se deberan hacer ahi los cambios para
# elegir entre galaxy, gassphere, lcdm_gas o cluster ...


H_PATH = ./general_libs/
#H_PATH = ../../../General_libs/

OBJS   = main.o  run.o  predict.o startrun.o global.o  \
	 timestep.o  init.o restart.o  io.o    \
	 accel.o   read_ic.o  ngb.o  \
	 system.o  allocate.o  density.o  \
	 gravtree.o hydra.o  driftfac.o  \
	 domain.o potential.o  \
         forcetree.o   peano.o gravtree_forcetest.o \
	 pm_periodic.o pm_nonperiodic.o longrange.o \
	 $(H_PATH)general/clib.o \
	 $(H_PATH)io/inout.o $(H_PATH)math/mathfns.o	\
	 $(H_PATH)math/mathutil.o $(H_PATH)general/getparam.o \
	 $(H_PATH)math/numrec.o $(H_PATH)general/lic.o \
	 $(H_PATH)mpi/endrun_mpi.o

#	 $(H_PATH)general/clib.o \
#	 $(H_PATH)mpi/endrun_mpi.o \
#	 $(H_PATH)io/inout.o $(H_PATH)math/mathfns.o \
#	 $(H_PATH)math/mathutil.o $(H_PATH)general/getparam.o \
#	 $(H_PATH)math/numrec.o $(H_PATH)general/lic.o


INCL   = globaldefs.h  protodefs.h  tags.h  Makefile \
	cmdline_defs.h switches.h \
	 $(H_PATH)general/constant.h \
	 $(H_PATH)io/inout.h $(H_PATH)general/machines.h \
	 $(H_PATH)math/mathfns.h $(H_PATH)math/mathutil.h $(H_PATH)general/getparam.h \
	 $(H_PATH)math/numrec.h $(H_PATH)general/param.h $(H_PATH)physics/physcons.h \
	 $(H_PATH)general/precision.h $(H_PATH)general/stdinc.h \
	 $(H_PATH)general/strings.h \
	 $(H_PATH)general/times.h $(H_PATH)physics/eosparam.h \
	 $(H_PATH)physics/units.h \
	 $(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h \
	 $(H_PATH)general/switchs.h $(H_PATH)general/lic.h \
	 $(H_PATH)mpi/mpi_proto.h

#	$(H_PATH)general/getparam.h $(H_PATH)general/stdinc.h $(H_PATH)general/lic.h \
#	$(H_PATH)math/mathfns.h $(H_PATH)math/numrec.h $(H_PATH)general/strings.h \
#	$(H_PATH)mpi/mpi_proto.h \
#	$(H_PATH)general/machines.h $(H_PATH)io/inout.h $(H_PATH)math/mathutil.h \
#	$(H_PATH)math/vectdefs.h $(H_PATH)math/vectmath.h $(H_PATH)general/switchs.h \
#	$(H_PATH)general/times.h $(H_PATH)general/param.h $(H_PATH)general/constant.h 


CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL)


ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif


LIBS   =   $(HDF5LIB) -g  $(MPICHLIB)  $(GSL_LIBS) -lgsl -lgslcblas -lm $(FFTW_LIB) 

#$(EXEC): $(OBJS) 
#	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

EXEC = gdgt207galaxy

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
#ifeq ($(simulation),gassphere)
#	rm -f $(OBJS) $(EXEC)
#else
#ifeq ($(simulation),galaxy)
	rm -f $(OBJS) $(EXEC)
#else
#ifeq ($(simulation),cluster)
#	rm -f $(OBJS) $(EXEC)
#else
#ifeq ($(simulation),lcdm_gas)
#	rm -f $(OBJS) $(EXEC)
#else
#	@echo " "
#	@echo "No simulation type selected!"
#	@echo "use with the option: simulation=xxx, where xxx is one of"
#	@echo "gassphere, galaxy, cluster, lcdm_gas"
#	@echo " "
#endif
#endif
#endif
#endif
