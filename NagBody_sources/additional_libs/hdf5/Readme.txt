All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


(The last version was downloaded from:
https://www.hdfgroup.org/downloads/hdf5/source-code/)

hdf5 INSTALLATION

Dependencies:

1. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

2. The easy way is using NagBody_pkg patterns:

cd $HOME/NagBody_pkg
make -f NagBody install_hdf5

Cleaning and packing are:

make -f NagBody clean_hdf5
make -f NagBody packing_nagbody_hdf5

3. We assume you are using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following lines:

export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib:${DYLD_LIBRARY_PATH}

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

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/cfitsio

6. Unpack file:

tar xvf hdf5-1.12.2.tar.gz

7. Change to directory:

cd hdf5-1.12.2

Note: read README

8. Configure, make and install:

./configure --prefix=$HOME/NagBody_pkg/local/hdf5 2>&1 | tee configure_hdf5.log
make 2>&1 | tee make_hdf5.log
make install 2>&1 | tee install_hdf5.log

9. Test:

make test

10. Clean directory:

make distclean

