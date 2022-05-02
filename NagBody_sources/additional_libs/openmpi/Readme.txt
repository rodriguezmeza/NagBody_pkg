All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


(The last version was downloaded from:
https://www.open-mpi.org/software/ompi/v4.1/)

openmpi INSTALLATION

Dependencies: fortran
(
For Mac OS X
sudo port install gcc46 +gfortran

Version can be different. Select the apropriate.
)

1. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

2. The easy way is using NagBody_pkg patterns:

cd $HOME/NagBody_pkg
make -f NagBody install_openmpi

Cleaning and packing are:

make -f NagBody clean_openmpi
make -f NagBody packing_nagbody_openmpi

3. We assume you are using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following lines:

export PATH=${HOME}/NagBody_pkg/local/openmpi/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/openmpi/share/man:${MANPATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/openmpi/lib:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/openmpi/lib/openmpi:${DYLD_LIBRARY_PATH}

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

For intel compilers:
export CC=icc
export CXX=icc
export F77=ifort
export FC=ifort
export F90=ifort


5. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/openmpi

6. Unpack file:

tar xvf openmpi-4.1.3.tar.gz

7. Change to directory:

cd openmpi-4.1.3

Note: read README

8. Configure, make and install:

./configure --prefix=$HOME/NagBody_pkg/local/openmpi 2>&1 | tee configure_openmpi.log
make 2>&1 | tee make_openmpi.log
make install 2>&1 | tee install_openmpi.log

9. Test:

make test


10. Clean directory:

make distclean

