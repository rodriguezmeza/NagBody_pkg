#
# The compiler
# ------------
#
#FC = ifort
FC = gfortran
#
# The key words of the compiler 
# -----------------------------
#
# Activate the precompilation parameter -DDOOMP for calculations
# using openMP (this is mandatory for proper behavior of fftw2)
#
# Activate the precompilation parameter -DSIMPLE for avoiding
# conversion of the particles positions and weights to double precision
# to save memory. This option has not been tested sufficiently to be
# considered as "safe"
#
# Activate the precompilation parameter -DADD1US for
# fixing underscore problems while linking with fftw2.
#
# Activate big-endian/little-endian conversion compiling option when
# relevant to avoid problems while reading input simulations files.
#
#FCFLAGS = -O -fopenmp -DDOOMP -DSIMPLE
#FCFLAGS = -O -fopenmp -DDOOMP
#FCFLAGS = -O -DSIMPLE
FCFLAGS = -O 
#
# The path of FFTW2 (NOT 3) 
# -----------------
# 
# Be aware of subtleties with openMP: fftw2 should be installed with
# the right options
#
#LIBS = -L/softs/fftw/2.1.5-pureintel/lib -ldrfftw -ldfftw 
LIBS = -L$(HOME)/NagBody_pkg/local/fftw2/lib -ldrfftw -ldfftw 
# In principle, no need to change anything from here.
INC =

OBJS =  ../src/twopidef.o \
        ../src/fourier_tools.o \
        ../src/gadget_tools.o \
        ../src/fourier_taylor_tools.o \
        ../src/fourier_tools3D.o \
        ../src/powmes_common.o \
        ../src/powmes.o

all:   powmes

powmes:  $(OBJS)
	 $(FC) $(FCFLAGS)  $(OBJS) -o $@ $(INC) $(LIBS)

clean: 
	/bin/rm  *.o *.mod powmes

.SUFFIXES:
.SUFFIXES: .f90 .F90 .o
.f90.o:
	$(FC) $(FCFLAGS)  $(INC) -c $< -o $@
.F90.o:
	$(FC) $(FCFLAGS)  $(INC) -c $< -o $@
