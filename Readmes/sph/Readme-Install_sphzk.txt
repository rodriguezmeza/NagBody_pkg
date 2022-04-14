All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2009-2014  M.A. Rodriguez-Meza, Mexico, D.F.

SPHZK INSTALLATION

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


2. To install the binary and man page sphzk files execute the command

> make -f NagBody install_sphzk

With this command you will have a running binary and its man page. Now,
you may run 

> sphzk -help

to get the sphzk help menu or

man sphzk

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/sph/sphzk

and there run the sphzk code with the default options:

> sphzk



