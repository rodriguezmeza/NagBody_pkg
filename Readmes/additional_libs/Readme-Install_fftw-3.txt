All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


fftw-3 INSTALLATION

Dependencies: 

OpenMPI


Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.

The easy way for Mac OS X is:

sudo port install fftw-3 +openmpi
sudo port install fftw-3-single +openmpi


1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

> make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

If you installed a NagBody code before, this step is not necessary.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw3/lib:${DYLD_LIBRARY_PATH}

Then, refresh your terminal.

Note: 
In Mac OSX is .bash_profile.
In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv DYLD_LIBRARY_PATH ${HOME}/NagBody_pkg/local/fftw2/lib:${DYLD_LIBRARY_PATH}

Also note that in the file "env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is included in the profile files. So this step could not be necessary.

3. Define the following environment variables:

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

> make distclean

