All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2019  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

planck2018 INSTALLATION

Dependencies: 

-cfitsio
-lapack

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

3. To install the binary and man page camb files execute the command:

Edit Makefile to define fortran compiler to gfortran, define flags, paths to blas, lapack, and cfitsio libraries. Be careful with the order of lapack libraries, have to be: "-llapack -lblas".


Then, from $NAGBODYDIR execute:

> make -f NagBody install_planck

and test it:

> make -f NagBody install_planck_test


With this command you will have a running binary and its man page, and will be tested. Now, you may invoke the man page using:

> man planck

to get the corresponding man page. 

