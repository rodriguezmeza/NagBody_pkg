All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2012  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

cosmmomc INSTALLATION

Dependencies: 

-cfitsio (needed by camb)
-Healpix (needed by camb)
-likelihood_v4 (THIS MUST BE INSTALLED BEFORE CosmoMC)
-Planck clp-2.0 and data clp_2.0

See corresponding read me files in directory "Readmes" for instructions on how to install them.

1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:
with the instruction:

> make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is $HOME/NagBody_pkg.

2. Modify your profile file. We assume the use of bash shell. Then edit 
the .profile file and make it to contain the following lines

export PATH=${PATH}:${HOME}/NagBody_pkg/bin
export MANPATH=${MANPATH}:${HOME}/NagBody_pkg/bin

Then, refresh your terminal.

Note: In some linux machines the file is .bash_profile or .bashrc.
Or if you are using tcsh, the file is .tcshrc, the above lines have to be:

setenv PATH ${PATH}:${HOME}/NagBody_pkg/bin
setenv MANPATH ${MANPATH}:${HOME}/NagBody_pkg/bin

Also note that in the file "$NAGBODYDIR/env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is used by the profile files. So this step could not be necessary.

3. To install the binary and man page cosmmomc files execute the command

> make -f NagBody install_cosmmomc

With this command you will have a running binary and its man page. Now, you may invoke the man page using:

> man cosmmomc

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/CosmoMC_CAMB/cosmmomc

and there run the cosmmomc code with the default options:

> mpirun -np 4 cosmmomc params.ini > output &

When finish running you will have several files. Then you may run:

getdist distparams.ini



Note: take care of the Makefile in the "cosmic/src/source" directory. It must have well defined the lapack, cfitsio, likelihood (WMAP), gls libraries.
