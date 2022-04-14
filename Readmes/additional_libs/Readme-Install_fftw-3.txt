All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2013  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

fftw-3 INSTALLATION

Dependencies: 

OpenMPI

See corresponding read me file in directory "Readmes" for instructions on how to install it.

$(NAGBODYDIR)/Readmes/Additional_libs/Readme-Install_openmpi


Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.


1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

> make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

If you installed a NagBody code before, this step is not necessary.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw2/lib:${DYLD_LIBRARY_PATH}

Then, refresh your terminal.

Note: In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv DYLD_LIBRARY_PATH ${HOME}/NagBody_pkg/local/fftw2/lib:${DYLD_LIBRARY_PATH}

Also note that in the file "env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is included in the profile files. So this step could   not be necessary.

3. Define the following environment variables:

For a standard Linux and Mac OS X:
export CC=gcc
export CXX=g++
export F77=gfortran
export FC=gfortran
export F90=gfortran

For a specific Mac OS X:
export CC=gcc-mp-4.6
export CXX=g++-mp-4.6
export F77=gfortran-mp-4.6
export FC=gfortran-mp-4.6
export F90=gfortran-mp-4.6

export CC=gcc-mp-4.9
export CXX=g++-mp-4.9
export F77=gfortran-mp-4.9
export FC=gfortran-mp-4.9
export F90=gfortran-mp-4.9

VIVI:
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

> cd $HOME/NagBody_pkg/NagBody_sources/Additional_libs/Gadget

5. Unpack file:

> gunzip fftw-3.2.2.tar.gz
> tar xvf fftw-3.2.2.tar

6. Change to directory:

> cd fftw-3.2.2

7. Configure, make and install:

First double without prefix:
./configure --prefix=$HOME/NagBody_pkg/local/fftw3 --enable-mpi --enable-threads 2>&1 | tee configure_gcc_gfortran_openmpi.log
make 2>&1 | tee make_gcc_gfortran_openmpi.log
make check
make install 2>&1 | tee install_gcc_gfortran_openmpi.log

No FUNCIONA --enable-type-prefix ::::

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

8. Move the directory gel and make the link:

> mv fftw3 fftw-3.2.2_intel64
> ln -s fftw-3.2.2_intel64 fftw3


9. Clean directory:

> make distclean

