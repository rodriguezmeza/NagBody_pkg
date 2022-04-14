All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2009-2012  M.A. Rodriguez-Meza, Mexico, D.F.

Note: when you see below execution lines starting with character ">", this
symbol means the unix prompt and should not be type it in.

mgpt INSTALLATION

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

2. To install the binary and man page mgpt files execute the command
(See note below to install binary and man page files precompiled)

> make -f NagBody install_mgpt

With this command you will have a running binary and its man page. Now,
you may run 

> mgpt -help

to get the mgpt help menu or

> man mgpt

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/perturbations/mgpt

and there run the mgpt code with the default options:

> mgpt

When funish running you will have several files.

