All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


(La œltima versi—n fue bajada de:
http://healpix.jpl.nasa.gov/)

HEALPix is an acronym for Hierarchical Equal Area isoLatitude Pixelization of a sphere. As suggested in the name, this pixelization produces a subdivision of a spherical surface in which each pixel covers the same surface area as every other pixel. The figure below shows the partitioning of a sphere at progressively higher resolutions, from left to right. The green sphere represents the lowest resolution possible with the HEALPix base partitioning of the sphere surface into 12 equal sized pixels. The yellow sphere has a HEALPix grid of 48 pixels, the red sphere has 192 pixels, and the blue sphere has a grid of 768 pixels (~7.3 degree resolution). 


Healpix INSTALLATION

Dependencies: 

-cfitsio

See corresponding Readme file in directory $HOME/NagBody_pkg/NagBody_sources/additional_libs/cfitsio for instructions on how to install it.

1. Follow the NagBody_pkg patterns:

cd $HOME/NagBody_pkg
make -f NagBody install_gsl

Install step only untar the file Healpix_3.30_2015Oct08.tar.gz and moves directory Healpix_3.30 to $HOME/NagBody_pkg/local. Then steps 4--6 below are not necessary.

Cleaning and packing are:

make -f NagBody clean_healpix
make -f NagBody packing_nagbody_healpix

2. See beginning of README.md in $HOME/NagBody_pkg for the initial steps configuring NagBody_pkg.

3. Define the following environment variables:

For a standard Linux and Mac OS X:
export CC=gcc
export CXX=g++
export F77=gfortran
export FC=gfortran
export F90=gfortran

4. Go to directory: 

cd $HOME/NagBody_pkg/NagBody_sources/additional_libs/healpix

5. Unpack file:

tar xvf Healpix_3.30_2015Oct08.tar.gz

6. Move directory Healpix_3.30 to $HOME/NagBody_pkg/local

mv Healpix_3.30 $HOME/NagBody_pkg/local/.

7. Change to that directory:

cd $HOME/NagBody_pkg/local/Healpix_3.30

8. Configure, make and install:

./configure

and choose options (2): Almost accept default, except when asked for cfitsio libraries and headers, there, type in the correct paths: Instead of "$HOME/NagBody_pkg/local/cfitsio/lib". Here use the full path (for example: /Users/mar/NagBody_pkg/local/cfitsio/lib). 

Option (3): use gfortran. No suffix selected. Create directories. Select default flags compilation and optimization. C compiler: gcc-mp-4.9. Enter default flags for optimization and library archiving. Enter default cfitsio options. No pgplot option selected. Use the parallel implementation. Use Position Independent Compilation. And no shared lib option also.

Then

make clean (if it was old or used)
make distclean
make
make test

Note: do not clean directory.

THE FOLLOWING IT IS ALREADY DONE IN THE PROCESS ABOVE DURING OPTION (2) if you answer yes:

Then, in linux we edit .profile and .bashrc to copy the last two lines, similar to:

# modifications by HEALPixAutoConf 2.15a
[ -r /Users/mar/.healpix/3_30_Darwin/config ] && . /Users/mar/.healpix/3_30_Darwin/config

from .profile to .bashrc. Refresh terminal to get the new environment.

