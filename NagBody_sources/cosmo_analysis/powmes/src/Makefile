MACHINES_DIR = ../../../../Machines/
include $(MACHINES_DIR)machine.inc
#
# The compiler
# ------------
#
#FC = ifort
# VIVI: portland group
#FC = pgf90
#FC = gfortran-mp-4.6
#FC = gfortran-mp-7
#FC = gfortran-mp-9
FC = gfortran
#
# The key words of the compiler (optimization, etc)
# -----------------------------
#
# Activate the precompilation parameter OMP for calculations
# using openMP (-DOMP in most fortran90 compilers)
#
# Activate the precompilation parameter -DADD1US for
# fixing underscore problems while linking with fftw2.
#
# Activate big-endian/little-endian conversion compiling option when
# relevant to avoid problems while reading input simulations files.
#
#FCFLAGS=  -Vaxlib -CB -mcmodel=large -O -openmp -DOMP
FCFLAGS=  -O -fopenmp -DOMP
#
# The path of FFTW2 (NOT 3) 
# -----------------
# 
# Use -lfftw -lrfftw for standard calculation
#
# Use  -lrfftw_threads -lrfftw -lfftw_threads -lfftw -lpthread for
# calculations with openMP.
#
#LIBS = -lrfftw_threads -lrfftw -lfftw_threads -lfftw -lpthread
# VIVI: fftw2 compiled with pgf90
#LIBS = -L/home/mario/local/lib -lfftw -lrfftw -ldfftw -m
#LIBS = -L$(HOME)/NagBody_pkg/local/fftw2/lib -lfftw -lrfftw
#LIBS = -L/opt/local/lib/gcc6
LIBS = -L/opt/local/lib/gcc9
# TO SOLVE LINKING ERROR:
# ld: library not found for -lSystem
LIBS += -L/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib
#LIBS += -L$(HOME)/NagBody_pkg/local/fftw2/lib -lfftw -lrfftw
LIBS += -L/opt/local/lib -ldfftw -ldrfftw
#LIBS += -L/opt/local/lib -ldfftw -lrfftw

# In principle, no need to change anything from here.
# VIVI: fftw2 compiled with pgf90 … y puerto en /home/local ...
#INC =
#INC = -I/home/mario/local/include
#INC = -I$(HOME)/NagBody_pkg/local/fftw2/include
INC = -I/opt/local/include
INC += -I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include

OBJS =  twopidef.o \
	numrec_tools.o \
        random_number_tools.o \
        fourier_tools.o \
        gadget_ramses_tools.o \
        fourier_taylor_tools.o \
        fourier_tools3D.o \
        powmes_common.o \
        powmes.o

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
