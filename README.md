# NagBody_pkg
Numerical algoritms for general body dynamics

In your home:

git clone https://github.com/rodriguezmeza/NagBody_pkg.git

Installation:

$ cd $HOME/NagBody_pkg

It is assumed that the trajectory to NagBody_pkg is under the your home directory. Then

$ make -f NagBody install_dirs

to make directories: bin, local, man, tests and zip. Let us do

$ make -f NagBody install_nbody_n2

Then move to tests:

$ cd tests

and run:

$ ../bin/nbody_n2

To get help:

$ ../bin/nbody_n2 -help

or

$ man ../man/man1/nbody_n2.1

and return to NagBody_pkg dir:

$ cd ..

Ok. Now let source:

$ source env_config/nagbodyrc.sh

Go again to tests directory and do the above tests:

$ nbody_n2
$ nbody_n2 -help
$ man nbody_n2

You may see that now you do not need to add ../bin, etc. That is because in the env_config/nagbodyrc.sh there are some additions to the environment variables and you activated by sourcing it.



