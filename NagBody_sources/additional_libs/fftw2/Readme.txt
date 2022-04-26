All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


fftw-2 INSTALLATION

Dependencies: 

OpenMPI (or mpich)

See corresponding Readme file for instructions on how to install it:

$(NAGBODYDIR)/Readmes/additional_libs/.../Readme-Install_openmpi


1. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

2. The easy way is using NagBody_pkg patterns:

cd $HOME/NagBody_pkg
make -f NagBody install_fftw2

Cleaning and packing are:

make -f NagBody clean_fftw2
make -f NagBody packing_nagbody_fftw2

Or it may be installed in the other option for Mac OS X:

sudo port install fftw-2 +openmpi
sudo port install fftw-2-single +openmpi

But need to modify Makefiles of the codes that need fftw-2 to set location of includes and libs.

3. We assume you are using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following lines:

export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw2/lib:${DYLD_LIBRARY_PATH}

Then source env_config/nagbodyrc.sh.

That's it! 

No need to go further.

4. The other way to install the package is to follow the steps:

If necessary, define the following environment variables:

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

Intel compilers:
export CC=icc
export CXX=icc
export F77=ifort
export FC=ifort
export F90=ifort

MPI compilers:
export CC=mpicc
export CXX=mpicxx
export F77=mpif90
export FC=mpif90
export F90=mpif90

Also if MPI is not found (for possible location):
export PATH=/opt/local/libexec/openmpi-mp:${PATH}
export MANPATH=/opt/local/share/man:${MANPATH}
export DYLD_LIBRARY_PATH=/opt/local/lib/openmpi-mp:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=/opt/local/lib/openmpi-mp/openmpi:${DYLD_LIBRARY_PATH}


5. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/fftw2

6. Unpack file:

tar xvf fftw-2.1.5.tar.gz

7. Change to directory:

cd fftw-2.1.5

8. Configure, make and install:

First double without prefix:
./configure --prefix=$HOME/NagBody_pkg/local/fftw2 --enable-mpi 2>&1 | tee configure_gcc_gfortran_openmpi.log
make 2>&1 | tee make_gcc_gfortran_openmpi.log
make check
make install 2>&1 | tee install_gcc_gfortran_openmpi.log

Then in simple precision with prefix 's':
make distclean
./configure --prefix=$HOME/NagBody_pkg/local/fftw2 --enable-float --enable-type-prefix --enable-mpi 2>&1 | tee configure_gcc_gfortran_openmpi.log
make 2>&1 | tee make_gcc_gfortran_openmpi.log
make check
make install 2>&1 | tee install_gcc_gfortran_openmpi.log

We repeat for double precision
make distclean
./configure --prefix=$HOME/NagBody_pkg/local/fftw2 --enable-type-prefix --enable-mpi 2>&1 | tee configure_gcc_gfortran_openmpi.log
make 2>&1 | tee make_gcc_gfortran_openmpi.log
make check
make install 2>&1 | tee install_gcc_gfortran_openmpi.log

9. Clean directory:

make distclean


10. NOTA:
Instead of using export's, configuration can be done as:

./configure --prefix=$HOME/NagBody_pkg/local/fftw-2.1.5_intel CC=icc CXX=icc F77=ifort F90=ifort 2>&1 | tee configure.log
make
make install

then link as:

ln -s fftw-2.1.5_intel fftw2

and add the to the environment:

export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg//local/fftw2/lib:${DYLD_LIBRARY_PATH}

