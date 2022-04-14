All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2009-2011  M.A. Rodriguez-Meza, Mexico, D.F.

MD_BLJ INSTALLATION

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


2. To install the binary and man page md_blj files execute the command

> make -f NagBody install_md_blj

With this command you will have a running binary and its man page. Now,
you may run 

> md_blj -help

to get the md_blj help menu or

man md_blj

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/Cmplxfluids/md_blj

and there run the md_blj code with the default options:

> md_blj

When funish running you will have several files. You may also run as

> md_blj computeRdf=true

Run the nplot2d code as:

> nplot2d in=gdr.exp,rdf.dat ws=true,false symboltype=4 pj=false,true uc=1:2,1:2 linecolor=1

It is shown the simulation radial distribution function as a line and
the corresponding experimental radial distribution function as circles
(Experiment reference: J.L. Yarnell, et al., Phys. Rev. A 7, 2130 (1973)).


