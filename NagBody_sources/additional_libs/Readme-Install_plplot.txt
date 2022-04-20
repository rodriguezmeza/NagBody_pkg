All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


plplot INSTALLATION

Dependencies: swig, LASi, cmake 

NOTE: If plplot libraries can be installed using apt-get or by other means then just make the following link:

cd $HOME/NagBody_pkg/local
ln -s /usr plplot

and the steps below are not necesary!! However make sure xwin and eps or ps drivers were installed!!


Note: cmake can be installed using the command (linux)

Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.

sudo apt-get install cmake


1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

If you installed a NagBody code before, this step is not necessary.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

Se incluyen las siguientes variables de ambiente en env_config/nagbodyrc:

export PATH=${HOME}/local/plplot/bin:${PATH}
export MANPATH=${HOME}/local/plplot/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/local/plplot/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/local/plplot/lib:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/local/plplot/lib/plplot5.15.0/driversd:${DYLD_LIBRARY_PATH}

Then, refresh your terminal.

Note: 
In Mac OSX is .bash_profile.
In some linux machines the file is .bash_profile or .bashrc.
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

If swig was installed under NagBody_pkg:
export CMAKE_INCLUDE_PATH=$HOME/NagBody_pkg/local/swig/share/swig/2.0.11

4. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/

5. Unpack file:

tar xvf plplot-5.15.0.tar.gz

6. Change to directory:

cd plplot-5.15.0

(directory name may change with the version)


7. Configure, make and install:

export CC="gcc -O2"
export CXX="g++ -O2"
export FC="gfortran -O2"

mkdir build_dir 
cd build_dir

Install withour fortran, qt, wxWidgets, Python, Java ...

cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_DYNDRIVERS=OFF -DPLD_psttfc=OFF -DPLD_plmeta=OFF .. 2>&1 | tee cmake.out

make 2>&1 | tee make.out
make install 2>&1 | tee make_install.out

cd ..
rm -fR build_dir

--------------------------
Other options:

mkdir build_dir
cd build_dir

We may install without Python, lua (luac), psttf.

In linux Ubuntu:

cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_lectures/NagBody_pkg/local/plplot -DENABLE_DYNDRIVERS=OFF -DPLD_psttfc=OFF -DPLD_plmeta=OFF .. >& cmake.out


In Mac OS X:

cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_lectures//NagBody_pkg/local/plplot -DENABLE_DYNDRIVERS=OFF -DPLD_psttfc=OFF -DPLD_plmeta=OFF -DENABLE_python:BOOL=OFF -DENABLE_lua:BOOL=OFF -DPLD_luac=OFF -DENABLE_cxx:BOOL=OFF -DPLD_psttf=OFF .. 2>&1 | tee cmake.out


Then execute:

make 2>&1 | tee make.out
make install 2>&1 | tee make_install.out

cd ..
rm -fR build_dir
--------------------------



To test:

cd $HOME/NagBody_pkg/local/plplot/share/plplot5.15.0/examples
make test_interactive


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINUX UBUNTU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Install first cmake:

sudo apt-get install cmake

and the following packages:

sudo apt-get install libx11-dev
sudo apt-get install gnuplot
sudo apt-get install xorg-dev
sudo apt-get install liblasi-dev
sudo apt-get install swig
sudo apt-get install tcl


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options are found in:

http://www.miscdebris.net/plplot_wiki/index.php?title=CMake_options_for_PLplot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Other options are:

With fortran 77 not worked:
cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_java:BOOL=OFF -DPLD_wxwidgets:BOOL=OFF -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=ON -DENABLE_f95:BOOL=ON -DBUILD_TEST=ON -DENABLE_DYNDRIVERS=ON -DPLD_aqt:BOOL=OFF -DENABLE_itcl:BOOL=OFF -DENABLE_python:BOOL=OFF -DPLD_jpeg:BOOL=ON -DPLD_gif:BOOL=ON -DPLD_png:BOOL=ON -DPLD_pdf:BOOL=ON -DPLD_pstex:BOOL=ON -DPLD_plmeta:BOOL=ON -DPLD_xterm:BOOL=OFF ../ >& cmake.out


Without fortran 77:
cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_java:BOOL=ON -DPLD_wxwidgets:BOOL=ON -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=OFF -DENABLE_f95:BOOL=ON -DBUILD_TEST=ON -DENABLE_DYNDRIVERS=ON -DPLD_aqt:BOOL=OFF -DENABLE_itcl:BOOL=OFF -DENABLE_python:BOOL=ON -DPLD_jpeg:BOOL=ON -DPLD_gif:BOOL=ON -DPLD_png:BOOL=ON -DPLD_pdf:BOOL=ON -DPLD_pstex:BOOL=ON -DPLD_plmeta:BOOL=ON -DPLD_xterm:BOOL=ON ../ >& cmake.out



Instalamos sin fortran, qt, wxWidgets, Python, Java ...

cmake -DCMAKE_INSTALL_PREFIX=$HOME/local/plplot -DENABLE_java:BOOL=OFF -DPLD_wxwidgets:BOOL=OFF -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=OFF -DENABLE_f95:BOOL=OFF -DBUILD_TEST=ON ../ >& cmake.out

Con los driveres dinamicos deshabilitados ...
cmake -DCMAKE_INSTALL_PREFIX=$HOME/local/plplot -DENABLE_java:BOOL=OFF -DPLD_wxwidgets:BOOL=OFF -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=OFF -DENABLE_f95:BOOL=OFF -DBUILD_TEST=ON -DENABLE_DYNDRIVERS=OFF ../ >& cmake.out

Finalmente la version que instalamos fue:

Sin Java, wxWidgets, qt (la ver instalada es de 32bits, necesitamos ver 64bits), tiene f77 y f90, drivers dinamicos ...
No tiene python ... tiene tcl y tk, tiene octave ...

cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_java:BOOL=OFF -DPLD_wxwidgets:BOOL=OFF -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=ON -DENABLE_f95:BOOL=ON -DBUILD_TEST=ON -DENABLE_DYNDRIVERS=ON -DPLD_aqt:BOOL=OFF -DENABLE_itcl:BOOL=OFF -DENABLE_python:BOOL=OFF -DPLD_jpeg:BOOL=ON -DPLD_gif:BOOL=ON -DPLD_png:BOOL=ON -DPLD_pdf:BOOL=ON -DPLD_pstex:BOOL=ON -DPLD_plmeta:BOOL=ON -DPLD_xterm:BOOL=OFF ../ >& cmake.out

cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_java:BOOL=ON -DPLD_wxwidgets:BOOL=ON -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DENABLE_f77:BOOL=ON -DENABLE_f95:BOOL=ON -DBUILD_TEST=ON -DENABLE_DYNDRIVERS=ON -DPLD_aqt:BOOL=OFF -DENABLE_itcl:BOOL=OFF -DENABLE_python:BOOL=ON -DPLD_jpeg:BOOL=ON -DPLD_gif:BOOL=ON -DPLD_png:BOOL=ON -DPLD_pdf:BOOL=ON -DPLD_pstex:BOOL=ON -DPLD_plmeta:BOOL=ON -DPLD_xterm:BOOL=OFF ../ >& cmake.out

