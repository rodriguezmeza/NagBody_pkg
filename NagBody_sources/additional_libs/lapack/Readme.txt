All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.

http://www.netlib.org/lapack/

lapack INSTALLATION


1. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

2. The easy way is using NagBody_pkg patterns:

cd $HOME/NagBody_pkg
make -f NagBody install_lapack

Cleaning and packing are:

make -f NagBody clean_lapack
make -f NagBody packing_nagbody_lapack

3. We assume you are using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following line:

export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/lapack/lib:${LD_LIBRARY_PATH}

Then source env_config/nagbodyrc.sh:

cd $HOME/NagBody_pkg
source env_config/nagbodyrc.sh

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

5. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/lapack

6. Unpack file:

tar xvf lapack-3.3.0.tgz

7. Change to directory:

cd lapack-3.3.0

8. Configure, make and install:

cp INSTALL/make.inc.gfortran make.inc

Edit the file make.inc, to change platform, and if necessary the fortran compiler. And edit the file Makefile, to change lib definition.

make all 2>&1 | tee make_all.log
mkdir $HOME/NagBody_pkg/local/lapack
mkdir $HOME/NagBody_pkg/local/lapack/lib

mv blas_LINUX.a $HOME/NagBody_pkg/local/lapack/lib/libblas.a
mv lapack_LINUX.a $HOME/NagBody_pkg/local/lapack/lib/liblapack.a
mv tmglib_LINUX.a $HOME/NagBody_pkg/local/lapack/lib/libtmglib.a

Note: be careful with your platform definition. So instead of "_LINUX", if you have chosen "_MACOSX", then make the appropriate substitutions:

mv blas_OSX.a $HOME/NagBody_pkg/local/lapack/lib/libblas.a
mv lapack_OSX.a $HOME/NagBody_pkg/local/lapack/lib/liblapack.a
mv tmglib_OSX.a $HOME/NagBody_pkg/local/lapack/lib/libtmglib.a

9. Clean directory (it is recommended):

make clean


