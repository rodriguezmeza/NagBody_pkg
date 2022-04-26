All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


gsl INSTALLATION

1. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

2. The easy way is using NagBody_pkg patterns:

cd $HOME/NagBody_pkg
make -f NagBody install_gsl

Cleaning and packing are:

make -f NagBody clean_gsl
make -f NagBody packing_nagbody_gsl

Or it may be installed in the other option for Mac OS X:

sudo port install fftw-3 +openmpi
sudo port install fftw-3-single +openmpi

But need to modify Makefiles of the codes that need fftw-3 to set location of includes and libs.

3. We assume you are using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following lines:

export PATH=${HOME}/NagBody_pkg/local/gsl/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/gsl/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/gsl/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/gsl/lib:${DYLD_LIBRARY_PATH}

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

Intel compilers:
export CC=icc
export CXX=icc
export F77=ifort
export FC=ifort
export F90=ifort


5. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/gsl

6. Unpack file:

tar xvf gsl-latest.tar.gz

7. Change to directory:

cd gsl-2.7.1

8. Configure, make and install:

./configure --prefix=$HOME/NagBody_pkg/local/gsl 2>&1 | tee configure_gsl.log
make 2>&1 | tee make_gsl.log
make check 2>&1 | tee check_gsl.log
make install 2>&1 | tee install_gsl.log

9. Clean directory:

make distclean

