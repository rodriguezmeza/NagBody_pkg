################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_colas.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_cfitsio
make -f NagBody packing_fftw2
make -f NagBody packing_fftw3
make -f NagBody packing_gsl
make -f NagBody packing_Healpix
make -f NagBody packing_lapack
make -f NagBody packing_OpenBLAS
make -f NagBody packing_plplot

