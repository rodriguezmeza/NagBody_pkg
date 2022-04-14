All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2012  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

nbody_n2 INSTALLATION

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

Also note that in the file "$NAGBODYDIR/env_config/nagbodyrc", this environment variables were already defined and the first time you install a NagBody code this file is used by the profile files. So this step could not be necessary. In the case that it necessary do the following:

For a Mac OS X:

> cp .profile profile_nagbody_bak
> cat profile_nagbody_bak NagBody_pkg/env_config/profilerc .profile

For a Linux:

> cp .bashrc bashrc_nagbody_bak
> cat bashrc_nagbody_bak NagBody_pkg/env_config/profilerc > .bashrc

Open a new terminal to refresh your environment. 

3. Change to NagBody directory:

> cd NagBody_pkg

4. Tune the file "Machines/machine.inc" or "Machines/pmachines.inc" in order to use the right compilers and libraries. See in this directory the files in there for several examples. 

For now, in Mac OS X:

> cp Machines/mac_osx.in Machines/machine.inc
> cp Machines/pmac_osx_mpi.in Machines/pmachine.inc

and for Linux:

> cp Machines/linux.in Machines/machine.inc
> cp Machines/plinux.in Machines/pmachine.inc


Note: I you unpack the file of a NagBody code when you have already unpacked another one, normally this process overwrite some configuration files, like the ones in the directory "Machines". If you already tune file "Machine/machine.inc" or "Machine/pmachine.inc" to work with your compilers and libraries, then you have to tune the files again or repeat steps as described in 4, above.

5. To install the binary and man page nbody_n2 files execute the command

> make -f NagBody install_nbody_n2

With this command you will have a running binary and its man page. Now,
you may run 

> nbody_n2 -help

to get the nbody_n2 help menu or

> man nbody_n2

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/Gravitation/nbody_n2

and there run the nbody_n2 code with the default options:

> nbody_n2

When finish running you will have several files. You may see some
plots using the nplot2d and analysis_grav codes. Read the corresponding
man pages for info on how to use these codes.


