
Must be a directory named "files_out"

Edit files/run_parameters.dat

OutputDir                   files_out/
Nsample          256

which gives the total number of particles. This value gives 256^3 = 16,777,216

mpirun -np 4 lpicola files/run_parameters.dat

gadgetviewer files_out/example_filename_z0p000.2
