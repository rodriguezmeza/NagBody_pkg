All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


lapack INSTALLATION


1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

If you installed a NagBody code before, this step is not necessary.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/lapack/lib:${LD_LIBRARY_PATH}

Then, refresh your terminal.

Note: 
In Mac OSX is .bash_profile.
In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv LD_LIBRARY_PATH ${HOME}/NagBody_pkg/local/lapack/lib:${LD_LIBRARY_PATH}

Also note that in the file "$NAGBODYDIR/env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is included in the profile files. So this step could   not be necessary.

3. Define the following environment variables:

For a standard Linux and Mac OS X:
export CC=gcc
export CXX=g++
export F77=gfortran
export FC=gfortran
export F90=gfortran


4. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/

5. Unpack file:

tar xvf lapack-3.3.0.tgz

6. Change to directory:

cd lapack-3.3.0

7. Configure, make and install:

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

8. Clean directory (it is recommended):

make clean


