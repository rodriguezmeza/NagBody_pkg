
Must be a directory named "ICs".

Edit file 2lpt_Carmen.param

Nmesh 	64
Nsample 64

GlassTileFac    128

This value gives 2097152 total number of particles

Then

$ mpirun -np 1 ngenic ics.param

Number of processors must be equal to NumFilesWrittenInParallel. In this case "2".

cd ICs

$ gadgetviewer ICs/ics_Carmen_2001

Edit inputspec_ics_Carmen_2001.txt and comment with a "#" first line. And also remove the zeros at the end. Then

$ nplot2d in=ICs/inputspec_ics_Carmen_2001.txt,ICs/inputspec_ics_Carmen_2001.txt plottype=3 uc=1:2,3:4 ws=0,1 symbolsize=0.2

Shows spectrum of each on of the snaps.

Compare with

$ nplot2d in=input_spectrum_LasDamas.dat

which is the input spectrum, already in log-log.
