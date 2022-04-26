All INSTALL PARTICLE DYNAMICS PROJECT (NagBody)
Copyright (c) 2006-2022  M.A. Rodriguez-Meza, Mexico, D.F.

OpenBLAS dependencies: gmake. Install it:

sudo port install gmake

In a similar way in Linux with apt-get.

configure the package by executing in the OpenBLAS folder.

gmake CC=gcc FC=gfortran

Install the package via

make install PREFIX=$HOME/NagBody_pkg/local/OpenBLAS