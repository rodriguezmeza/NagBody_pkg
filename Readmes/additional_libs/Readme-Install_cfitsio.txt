All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


(The last version was downloaded from:
http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)

CFITSIO is a library of C and Fortran subroutines for reading and writing data files in FITS (Flexible Image Transport System) data format. CFITSIO provides simple high-level routines for reading and writing FITS files that insulate the programmer from the internal complexities of the FITS format. CFITSIO also provides many advanced features for manipulating and filtering the information in FITS files.

cfitsio INSTALLATION

Dependencies: fortran
(
For Mac OS X
sudo port install gcc46 +gfortran
)

Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.

1. In the $NAGBODYDIR make bin, man, man/man1, tests, zip, and local directories:
with the instruction:

make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

If you installed a NagBody code before, this step is not necessary.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib:${DYLD_LIBRARY_PATH}

Then, refresh your terminal.

Note: 
In Mac OSX is .bash_profile.
In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv PKG_CONFIG_PATH ${HOME}/NagBody_pkg/local/cfitsio/lib/pkgconfig:${PKG_CONFIG_PATH}
setenv DYLD_LIBRARY_PATH ${HOME}/NagBody_pkg/local/cfitsio/lib:${DYLD_LIBRARY_PATH}

Also note that in the file "$NAGBODYDIR/env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is included in the profile files. So this step could   not be necessary.

3. Could be necessary to define the following environment variables:

For a standard Linux and Mac OS X:
export CC=gcc
export CXX=g++
export F77=gfortran
export FC=gfortran
export F90=gfortran


4. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/

5. Unpack file:

tar xvf cfitsio3390.tar.gz

6. Change to directory:

cd cfitsio

Note: read README

7. Configure, make and install:

./configure --prefix=$HOME/NagBody_pkg/local/cfitsio 2>&1 | tee configure_gcc.log
make 2>&1 | tee make_gcc.log
make install 2>&1 | tee install_gcc.log

8. Test:

make testprog
./testprog > testprog.lis
diff testprog.lis testprog.out
cmp testprog.fit testprog.std


9. Clean directory:

> make distclean

