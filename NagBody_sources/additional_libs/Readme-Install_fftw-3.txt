All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


fftw-3 INSTALLATION

Dependencies: OpenMPI

1. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

2. The easy way for Mac OS X is:

sudo port install fftw-3 +openmpi
sudo port install fftw-3-single +openmpi

But need to modify Makefiles of the codes that need fftw-3 to set location of includes and libs.

2. We assume your using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following line

export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw3/lib:${DYLD_LIBRARY_PATH}

Then source env_config/nagbodyrc.sh.

3. If necessary, define the following environment variables:

For a standard Linux and Mac OS X:
export CC=gcc
export CXX=g++
export F77=gfortran
export FC=gfortran
export F90=gfortran

Portland group:
export CC=pgcc
export CXX=pgcc ?
export F77=pgf90
export FC=pgf90
export F90=pgf90

Intel compilers
export CC=icc
export CXX=icc
export F77=ifort
export FC=ifort
export F90=ifort


4. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/

5. Unpack file:

tar xvf fftw-3.3.7.tar.gz

6. Change to directory:

cd fftw-3.3.7

7. Configure, make and install:

First double without prefix:
./configure --prefix=$HOME/NagBody_pkg/local/fftw3 --enable-mpi --enable-threads 2>&1 | tee configure_gcc_gfortran_openmpi.log
make 2>&1 | tee make_gcc_gfortran_openmpi.log
make check
make install 2>&1 | tee install_gcc_gfortran_openmpi.log

Not WORKING: --enable-type-prefix ::::

Then in simple precision with prefix 's':
make distclean
./configure --prefix=$HOME/NagBody_pkg/local/fftw3 --enable-float --enable-type-prefix --enable-threads --enable-mpi 2>&1 | tee configure_gcc_gfortran_openmpi.log
make 2>&1 | tee make_gcc_gfortran_openmpi.log
make check
make install 2>&1 | tee install_gcc_gfortran_openmpi.log

We repeat for double precision
make distclean
./configure --prefix=$HOME/NagBody_pkg/local/fftw3 --enable-type-prefix --enable-threads --enable-mpi 2>&1 | tee configure_gcc_gfortran_openmpi.log
make 2>&1 | tee make_gcc_gfortran_openmpi.log
make check
make install 2>&1 | tee install_gcc_gfortran_openmpi.log


8. Clean directory:

make distclean

