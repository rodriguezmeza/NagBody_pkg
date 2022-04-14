All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2011  M.A. Rodriguez-Meza, Mexico, D.F.

nplot2d INSTALLATION

DEPENDENCIES: This code needs plplot libraries. Install them first. There are in the Additional_libs/plplot directory that is packed with nplot2d code. See the corresponding Readme file in the Readmes/Additional_libs directory.

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

2. To install the binary and man page nplot2d files execute the command

> make -f NagBody install_nplot2d

With this command you will have a running binary and its man page. Now,
you may run 

> nplot2d -help

to get the nplot2d help menu or

> man nplot2d

to get the corresponding man page. Or you may go to the directory

> cd NagBody_Tests/NagBody_PLPlots/nplot2d

and there run the nplot2d code with the default options:

> nplot2d inputfile=data.dat


ADDITIONAL NOTES:

1. Installing precompiled binaries and man page files. You just
change to directory above "NagBody_pkg" and axecute the command:

gunzip nplot2d_bin.tar.gz
tar xvf nplot2d_bin.tar

and you will get install the files.

2. If you are installing a NagBody code which needs plplot library
then you have to install it. Go to the directory where this library
is, normally it should be in: 

$NAGBODYDIR/NagBody_sources/NagBody_PLPlots/nplot2d

There you have to unzipped and unpacked

gunzip plplot-5.6.1.tar.gz
tar xvf plplot-5.6.1.tar 

Change to directory plplot-5.6.1 and run

./configure --prefix=$HOME/local/plplot-5.6.1
make
make install

Edit your .profile or .bashrc file as explained in sec. 2 above to
add the plplot path:

export PATH=${PATH}:${HOME}/bin:${HOME}/local/plplot/bin

and make the symbolic link in the directory $HOME/local:

ln -s plplot-5.6.1 plplot

Start a new terminal and that's it.


Can be the case that plplot library was installed using MacPorts (or apt-get in linux). Then simply use:

ln -s /opt/local plplot

Also edit file pldefs.h in $HOME/NagBody_pkg/NagBody_sources/General_libs/visual

// En la version plplot-5.9.10 _c_plwid aparece como simbolo no definido... debe cambiarse a c_plwidth
// /Users/mar/NagBody_pkg/local/plplot/include/plplot/plplot.h:755:0: note: this is the location of the previous definition
// #define    plwid                    c_plwid
//
#define plwid plwidth           // Version plplot-5.9.10 (NEW 2013_11_27)
//
// Ver el manual plplot-5.9.10.pdf para ver como se usa plwidth (Y no est‡ plwid)
//
