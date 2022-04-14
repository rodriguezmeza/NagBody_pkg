All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2010  M.A. Rodriguez-Meza, Mexico, D.F.

analysis_galaxy INSTALLATION

1. In the $NAGBODYDIR make bin and man directories:

mkdir bin
mkdir man
mkdir man/man1

Here $NAGBODYDIR is the directory which contains the NagBody files. 
For example, if you unpacked the NagBody zipped file in your $HOME dir
then the $NAGBODYDIR is NagBody_pkg.

In your home directory ($HOME) make the simbolic links:

ln -s $NAGBODYDIR/bin bin
ln -s $NAGBODYDIR/man man

For example, if you unpacked the NagBody zipped file in your $HOME dir
then the above commands are:

ln -s NagBody_pkg/bin bin
ln -s NagBody_pkg/man man

2. Modify your profile file. We assume the use of bash SHELL. Then edit 
the .profile file and make it to contain the following lines

export PATH=${PATH}:${HOME}/bin

Note: In some linux machines the file is .bash_profile.

If you are going to use gnuplot under X Windows edit .bashrc file to 
contain the same line above.

3. To install the binary and man page analysis_galaxy files execute the command

make -f NagBody install_analysis_galaxy

With this command you will have a running binary and its man page. Now,
you may run 

analysis_galaxy -help

to get the analysis_galaxy help menu or

man analysis_galaxy

to get the corresponding man page. Or you may go to the directory

cd NagBody_Tests/Gravitation/analysis_galaxy

and there run the analysis_galaxy code with the default options:

analysis_galaxy

