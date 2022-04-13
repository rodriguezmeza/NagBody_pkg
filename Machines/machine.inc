#
# Include file for building exe's for mac_osx
#
# C compilers:
CC = gcc
#CC = clang
#CC = gcc-4.2
#CC = gcc-mp-4.4
#CC = gcc-mp-4.6
#CC = gcc-mp-4.9
#CC = icc

# Fortran compilers:
#FC = gfortran-mp-4.6
#FC = gfortran-mp-4.9
FC = gfortran

#F90C = ifort
#F90C = gfortran-mp-4.6
#F90C = gfortran-mp-4.9
F90C = gfortran


#F90 = gfortran-mp-4.6
#F90 = gfortran-mp-4.9
F90 = gfortran

# If it is unsupported option
#CFLAGS = -O3 -fopenmp 
CFLAGS = -O3
#CFLAGS = -O3 -fopenmp -I/opt/local/lib/gcc8/gcc/x86_64-apple-darwin18/8.3.0/include/ssp/

OPTIMIZE = -O2
LIBS = -lm
