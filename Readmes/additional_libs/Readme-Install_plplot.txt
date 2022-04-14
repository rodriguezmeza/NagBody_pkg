All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2012  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

plplot INSTALLATION

Dependencies: swig, LASi, cmake 

NOTE: If plplot libraries can be installed using apt-get or by other means then just make the following link:

> cd $HOME/NagBody_pkg/local
> ln -s /usr plplot

and the steps below are not necesary!! However make sure xwin and eps or ps drivers were installed!!


Note: cmake can be installed using the command (linux)

Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.

> sudo apt-get install cmake


1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

> make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

If you installed a NagBody code before, this step is not necessary.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export PATH=${HOME}/NagBody_pkg/local/plplot/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/plplot/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/plplot/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib/plplot5.9.10/driversd:${DYLD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib/plplot5.9.10/driversd:${LD_LIBRARY_PATH}

Then, refresh your terminal.

Note: In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv PATH ${HOME}/NagBody_pkg/local/plplot/bin:${PATH}
setenv MANPATH ${HOME}/NagBody_pkg/local/plplot/share/man:${MANPATH}
setenv PKG_CONFIG_PATH ${HOME}/NagBody_pkg/local/plplot/lib/pkgconfig:${PKG_CONFIG_PATH}
setenv DYLD_LIBRARY_PATH ${HOME}/NagBody_pkg/local/plplot/lib:${DYLD_LIBRARY_PATH}
setenv DYLD_LIBRARY_PATH ${HOME}/NagBody_pkg/local/plplot/lib/plplot5.9.5/driversd:${DYLD_LIBRARY_PATH}

Also note that in the file "$NAGBODYDIR/env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is included in the profile files. So this step could   not be necessary.

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
export CMAKE_INCLUDE_PATH=$HOME/NagBody_pkg/local/swig/share/swig/2.0.11

4. Go to directory: 

> cd $HOME/NagBody_pkg/NagBody_sources/Additional_libs/plplot

5. Unpack file:

> unzip plplot-5.9.7.tar.gz
> tar xvf plplot-5.9.7.tar

6. Change to directory:

> cd plplot-5.9.7

7. Configure, make and install:

mkdir build_dir 
cd build_dir

Con fortran 77 no funcion—:
cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_java:BOOL=OFF -DPLD_wxwidgets:BOOL=OFF -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=ON -DENABLE_f95:BOOL=ON -DBUILD_TEST=ON -DENABLE_DYNDRIVERS=ON -DPLD_aqt:BOOL=OFF -DENABLE_itcl:BOOL=OFF -DENABLE_python:BOOL=OFF -DPLD_jpeg:BOOL=ON -DPLD_gif:BOOL=ON -DPLD_png:BOOL=ON -DPLD_pdf:BOOL=ON -DPLD_pstex:BOOL=ON -DPLD_plmeta:BOOL=ON -DPLD_xterm:BOOL=OFF ../ >& cmake.out


Sin fortran 77:
cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_java:BOOL=ON -DPLD_wxwidgets:BOOL=ON -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=OFF -DENABLE_f95:BOOL=ON -DBUILD_TEST=ON -DENABLE_DYNDRIVERS=ON -DPLD_aqt:BOOL=OFF -DENABLE_itcl:BOOL=OFF -DENABLE_python:BOOL=ON -DPLD_jpeg:BOOL=ON -DPLD_gif:BOOL=ON -DPLD_png:BOOL=ON -DPLD_pdf:BOOL=ON -DPLD_pstex:BOOL=ON -DPLD_plmeta:BOOL=ON -DPLD_xterm:BOOL=ON ../ >& cmake.out


make >& make.out
make install >& make_install.out

mv plplot/ plplot-5.9.10_gcc-mp-4.9_gfortran-mp-4.9
ln -s plplot-5.9.10_gcc-mp-4.9_gfortran-mp-4.9/ plplot

cd ..
rm -fR build_dir


Para hacer pruebas ...

mkdir $HOME/Research/borrame/plplot
cd $HOME/Research/borrame/plplot
cp -fR /Users/mar/NagBody_pkg/local/plplot/share/plplot5.9.7/examples .
cd examples
make test_interactive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINUX UBUNTU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Primero se instala cmake con:

sudo apt-get install cmake

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options are found in:

http://www.miscdebris.net/plplot_wiki/index.php?title=CMake_options_for_PLplot