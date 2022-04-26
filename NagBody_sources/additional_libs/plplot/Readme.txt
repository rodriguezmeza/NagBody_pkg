All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


plplot INSTALLATION

===================================================
Dependencies: swig, LASi, pango, cmake, qt5, tcl...

In Mac OS X:

sudo port install cmake
sudo port install swig
sudo port install qt5
(
--->  Some of the ports you installed have notes:
 dbus has the following notes:
    Startup items (named 'dbus-system, dbus-session') have been generated that will aid in starting dbus with launchd. They are disabled by default. Execute the
    following command to start them, and to cause them to launch at startup:
    
        sudo port load dbus
  python310 has the following notes:
    To make this the default Python or Python 3 (i.e., the version run by the 'python' or 'python3' commands), run one or both of:
    
        sudo port select --set python python310
        sudo port select --set python3 python310
)
(
--->  Some of the ports you installed have notes:
  libomp has the following notes:
    To use this OpenMP library:
     * For clang-3.8+, or clang-3.7 with +openmp variant:
        add "-fopenmp" during compilation / linking.
     * For clang-3.7 without +openmp variant, use:
        "-I/opt/local/include/libomp -L/opt/local/lib/libomp -fopenmp"
)

(
--->  Some of the ports you installed have notes:
  tcl has the following notes:
    The Sqlite3 Tcl package is now being provided by the sqlite3-tcl port:
    sudo port install sqlite3-tcl
)


In linux Ubuntu:

sudo apt-get install cmake

and the following packages:

sudo apt-get install libx11-dev
sudo apt-get install gnuplot
sudo apt-get install xorg-dev
sudo apt-get install liblasi-dev
sudo apt-get install swig
sudo apt-get install tcl
===================================================


===================================================
Note: If plplot libraries can be installed using apt-get or by other means then just make the following link:

cd $HOME/NagBody_pkg/local
ln -s /usr plplot

In Mac OS X with Ports:

ln -s /opt/local plplot
(OLD one: plplot -> /opt/local/lib/plplot510/)


and the steps below are not necesary!! However make sure xwin, jpg (or png) and pdf drivers in plplot were installed!!
===================================================


===================================================
1. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

2. The easy way is using NagBody_pkg patterns:

cd $HOME/NagBody_pkg
make -f NagBody install_plplot

Cleaning and packing are:

make -f NagBody clean_plplot
make -f NagBody packing_nagbody_plplot

3. We assume you are using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following lines:

export PATH=${HOME}/local/plplot/bin:${PATH}
export MANPATH=${HOME}/local/plplot/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/local/plplot/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/local/plplot/lib:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/local/plplot/lib/plplot5.15.0/driversd:${DYLD_LIBRARY_PATH}

Then source env_config/nagbodyrc.sh:

cd $HOME/NagBody_pkg
source env_config/nagbodyrc.sh

That's it! 

No need to go further.

===================================================
4. The other way to install the package is to follow the steps:

If necessary, define the following environment variables:

For a standard Linux and Mac OS X:
export CC=gcc
export CXX=g++
export F77=gfortran
export FC=gfortran
export F90=gfortran

If swig was installed under NagBody_pkg:
export CMAKE_INCLUDE_PATH=$HOME/NagBody_pkg/local/swig/share/swig/2.0.11

5. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/plplot

6. Unpack file:

tar xvf plplot-5.15.0.tar.gz

7. Change to directory:

cd plplot-5.15.0

(directory name may change with the version)


8. Configure, make and install:

export CC="gcc -O2"
export CXX="g++ -O2"
export FC="gfortran -O2"

mkdir build_dir 
cd build_dir

cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DDEFAULT_NO_QT_DEVICES:BOOL=ON -DENABLE_qt:BOOL=OFF -DPLD_tiffqt=OFF -DPLD_bmpqt=OFF -DPLD_jpgqt=OFF -DPLD_pngqt=OFF -DPLD_epsqt=OFF -DPLD_pdfqt=OFF -DPLD_ppmqt=OFF -DPLD_svgqt=OFF -DPLD_qtwidget=OFF -DPLD_extqt=OFF -DBUILD_TEST=ON -DENABLE_python:BOOL=OFF .. 2>&1 | tee cmake.out

make 2>&1 | tee make.out
make install 2>&1 | tee make_install.out

To test:

cd $HOME/NagBody_pkg/local/plplot/share/plplot5.15.0/examples
make test_interactive

cd ..
rm -fR build_dir
===================================================




%%%%%%%%%%%%%%%%%%%% 
Other possibilities to configure:

Install withour fortran, qt, wxWidgets, Python, Java ...

mkdir build_dir 
cd build_dir

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

cmake -DCMAKE_INSTALL_PREFIX=$HOME/NagBody_pkg/local/plplot -DENABLE_DYNDRIVERS=OFF -DPLD_psttfc=OFF -DPLD_plmeta=OFF .. >& cmake.out


In Mac OS X:

cmake -DCMAKE_INSTALL_PREFIX=$HOME//NagBody_pkg/local/plplot -DENABLE_DYNDRIVERS=OFF -DPLD_psttfc=OFF -DPLD_plmeta=OFF -DENABLE_python:BOOL=OFF -DENABLE_lua:BOOL=OFF -DPLD_luac=OFF -DENABLE_cxx:BOOL=OFF -DPLD_psttf=OFF .. 2>&1 | tee cmake.out


Then execute:

make 2>&1 | tee make.out
make install 2>&1 | tee make_install.out

cd ..
rm -fR build_dir
--------------------------


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

