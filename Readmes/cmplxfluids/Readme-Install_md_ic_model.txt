All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2009-2012  M.A. Rodriguez-Meza, Mexico, D.F.

MD_IC_MODEL INSTALLATION

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


2. To install the binary and man page md_ic_model files execute the command

> make -f NagBody install_md_ic_model

With this command you will have a running binary and its man page. Now,
you may run 

> md_ic_model -help

to get the md_ic_model help menu or

man md_ic_model

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/Cmplxfluids/md_ic_model

and there run the md_ic_model code with the default options:

> md_ic_model

When funish running you will have several files.


