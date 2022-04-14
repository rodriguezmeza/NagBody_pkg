All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2010  M.A. Rodriguez-Meza, Mexico, D.F.

model INSTALLATION

DEPENDENCIES: This code needs gsl libraries. Install them first. They are in the “gadget207_Addlibs.tar.gz” package. See the corresponding Readme file.


Note: steps 1 and 2 are not necessary if you have installed another NagBody code already.


1. In the $NAGBODYDIR make bin, man, man/man1, zip, and local directories:

mkdir bin
mkdir man
mkdir man/man1
mkdir zip
mkdir local

with the instruction

make -f NagBody install_dirs

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

export PATH=${PATH}:${HOME}/NagBody_pkg/bin

Note: In some linux machines the file is .bash_profile.

If you are going to use gnuplot under X Windows edit .bashrc file to 
contain the same line above.

3. To install the binary and man page model files execute the command

make -f NagBody install_model

With this command you will have a running binary and its man page. Now,
you may run 

model -help

to get the model help menu or

man model

to get the corresponding man page. Or you may go to the directory

cd NagBody_Tests/Gravitation/model

and there run the models code with the default options:

model

