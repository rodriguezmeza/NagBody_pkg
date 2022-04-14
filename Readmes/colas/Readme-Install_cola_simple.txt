All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2012  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

cola_only INSTALLATION

Dependencies: 

-fftw-2
-gls

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

3. To install the binary and man page cola_only files execute the command

> make -f NagBody cola_only

With this command you will have a running binary and its man page. 

Now, you may go to directory NagBody_Tests/colas/cola_only run 

> mpirun -np 2 cola_only in.param > output &




