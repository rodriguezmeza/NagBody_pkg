
Must be a directory named "ICs".

Edit file ics.param

Nmesh 	64
Nsample 64

Then

$ mpirun -np 2 ngenic ics.param

Number of processors must be equal to NumFilesWrittenInParallel. In this case "2".

cd ICs

$ gadgetviewer ics.0
$ gadgetviewer ics.1

$ nplot2d in=inputspec_ics.txt,inputspec_ics.txt plottype=3 uc=1:2,3:4 ws=0,1 symbolsize=0.2

Shows spectrum of each on of the snaps.