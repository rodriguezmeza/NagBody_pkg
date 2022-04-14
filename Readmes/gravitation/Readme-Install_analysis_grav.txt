All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2011  M.A. Rodriguez-Meza, Mexico, D.F.

analysis_grav INSTALLATION

Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.

1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories
with the instruction:

> make -f NagBody install_dirs

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is NagBody_pkg.

In your home directory ($HOME) make:

For a Mac OS X:

> cp .profile profile_nagbody_bak
> cat profile_nagbody_bak NagBody_pkg/env_config/profilerc .profile

For a Linux:

> cp .bashrc bashrc_nagbody_bak
> cat bashrc_nagbody_bak NagBody_pkg/env_config/profilerc > .bashrc


2. To install the binary and man page analysis_grav files execute the command

> make -f NagBody install_analysis_grav

With this command you will have a running binary and its man page. Now,
you may run 

> analysis_grav -help

to get the analysis_grav help menu or

> man analysis_grav

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/Gravitation/analysis_grav

and there run the analysis_grav code with the default options:

> analysis_grav

