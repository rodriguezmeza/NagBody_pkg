#EXEC   = 2LPTic
EXEC   = 2lptic

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile



#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing


#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type

OPT   +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles start from a glass (as opposed to grid)

#OPT += -DONLY_ZA # swith this on if you want ZA initial conditions (2LPT otherwise)



OPTIONS =  $(OPT)


#FFTW_INCL = -I/home/rs123/local/include
#FFTW_LIBS = -L/home/rs123/local/lib
FFTW_INCL= -I$(HOME)/NagBody_pkg/local/fftw2/include
FFTW_LIBS= -L$(HOME)/NagBody_pkg/local/fftw2/lib

#GSL_LIBS=   -L/usr/local/gsl/1.9/lib/
#GSL_INCL =  -I/usr/local/gsl/1.9/include/
GSL_INCL =  -I$(HOME)/NagBody_pkg/local/gsl/include
GSL_LIBS =  -L$(HOME)/NagBody_pkg/local/gsl/lib


#CC       =  /usr/mpi/gcc/openmpi-1.2.2-1/bin/mpicc
CC       =  mpicc   
#MPICHLIB = -L/usr/mpi/gcc/openmpi-1.2.2-1/lib64 -lmpi
MPICHLIB = -lmpi

OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)


FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



