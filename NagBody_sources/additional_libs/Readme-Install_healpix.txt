All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.


(La œltima versi—n fue bajada de:
http://healpix.jpl.nasa.gov/)

HEALPix is an acronym for Hierarchical Equal Area isoLatitude Pixelization of a sphere. As suggested in the name, this pixelization produces a subdivision of a spherical surface in which each pixel covers the same surface area as every other pixel. The figure below shows the partitioning of a sphere at progressively higher resolutions, from left to right. The green sphere represents the lowest resolution possible with the HEALPix base partitioning of the sphere surface into 12 equal sized pixels. The yellow sphere has a HEALPix grid of 48 pixels, the red sphere has 192 pixels, and the blue sphere has a grid of 768 pixels (~7.3 degree resolution). 


Healpix INSTALLATION

Dependencies: 

-cfitsio

See corresponding read me file in directory "Readmes" for instructions on how to install it.

Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.


1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

If you installed a NagBody code before, this step is not necessary.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export PATH=${PATH}:${HOME}/NagBody_pkg/bin
export MANPATH=${MANPATH}:${HOME}/NagBody_pkg/bin

Then, refresh your terminal.

Note: 
In Mac OSX is .bash_profile.
In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv PATH ${PATH}:${HOME}/NagBody_pkg/bin
setenv MANPATH ${MANPATH}:${HOME}/NagBody_pkg/bin

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

tar xvf Healpix_3.30_2015Oct08.tar.gz

6. Move directory Healpix_3.30 to $HOME/NagBody_pkg/local

mv Healpix_3.30 $HOME/NagBody_pkg/local/.

6. Change to that directory:

cd $HOME/NagBody_pkg/local/Healpix_3.30

7. Configure, make and install:

./configure

and choose options (2): Almost accept default, except when asked for cfitsio libraries and headers, there, type in the correct paths: Instead of "$HOME/NagBody_pkg/local/cfitsio/lib". Here use the full path (for example: /Users/mar/NagBody_pkg/local/cfitsio/lib). 


Option (3): use gfortran-mp-4.9. No suffix selected. Create directories. Select default flags compilation and optimisation. C compiler: gcc-mp-4.9. Enter default flags for optimisation and library archiving. Enter default cfitsio options. No pgplot option selected. Use the parallel implementation. Use Position Independent Compilation. And no shared lib option also.


Then

make clean (if it was old or used)
make
make test

Note: do not clean directory.


THE FOLLOWING IT IS ALREADY DONE IN THE PROCESS ABOVE DURING OPTION (2) if you answer yes:

Then, in linux we edit .profile and .bashrc to copy the last two lines, similar to:

# modifications by HEALPixAutoConf 2.15a
[ -r /Users/mar/.healpix/2_15a_Darwin/config ] && . /Users/mar/.healpix/2_15a_Darwin/config

from .profile to .bashrc. Refresh terminal to get the new environment.

