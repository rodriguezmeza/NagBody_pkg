All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2014  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

gdgt207lcdm INSTALLATION

Dependencies: 

gsl-1.x
openmpi-1.4.x
fftw-2.x.x
szip-2.x
zlib-1.2.3
hdf5-1.8.x

gsl, szip, and zlib can be installed using "apt-get" (Linux) or "port" (Mac). HDF5 is optional and also szip and zlib. To install openmpi-1.4.x and fftw-2.x.x see the corresponding readme files.

1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

> make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export PATH=${PATH}:${HOME}/NagBody_pkg/bin
export MANPATH=${MANPATH}:${HOME}/NagBody_pkg/man

Then, refresh your terminal.

Note: In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv PATH ${PATH}:${HOME}/NagBody_pkg/bin
setenv MANPATH ${MANPATH}:${HOME}/NagBody_pkg/man

Also note that in the file "$NAGBODYDIR/env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is used by the profile files. So this step could not be necessary.

3. To install the binary and man page gdgt207lcdm files execute the command

> make -f NagBody clean_gdgt207lcdm
> make -f NagBody install_gdgt207lcdm

With this command you will have a running binary and its man page. Now, you may run 

> mpirun -np 4 gdgt207lcdm parameterfile > output &

The parameter file must be consistent with the type of simulation you generate.

You may test the code É go to directory "NagBody_Tests/gdgt2_derived_codes/gdgt207lcdm/". There read the readme file for instructions.



