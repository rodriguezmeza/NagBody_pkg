
# Mar
#FFTW_INCL= -I/opt/local/include
#FFTW_LIBS= -L/opt/local/lib -lfftw3f
#GSL_INCL =  -I/opt/local/include
#GSL_LIBS =  -L/opt/local/lib -lgsl -lgslcblas
# NagBody_pkg:
FFTW_INCL= -I$(HOME)/NagBody_pkg/local/fftw3/include
FFTW_LIBS= -L$(HOME)/NagBody_pkg/local/fftw3/lib -lfftw3f
GSL_INCL =  -I$(HOME)/NagBody_pkg/local/gsl/include
GSL_LIBS =  -L$(HOME)/NagBody_pkg/local/gsl/lib -lgsl -lgslcblas


CFLAGS= $(FFTW_INCL) $(GSL_INCL) -std=c99 -O3 -Wall 
LIBS= -lm $(FFTW_LIBS) $(GSL_LIBS) 

EXEC = cola

# Mar
COLAcode :  auxPM.o gadget_blob.o cola.o
		gcc ${CFLAGS} -o $(EXEC) main.c  auxPM.o gadget_blob.o cola.o ${LIBS}

auxPM.o :
		gcc ${CFLAGS} -c auxPM.c 

gadget_blob.o :
		gcc ${CFLAGS} -c gadget_blob.c 

# Mar
cola.o :
		gcc ${CFLAGS} -c cola.c

clean:
	rm -f *.o $(EXEC) 
