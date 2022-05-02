################################################################################
#NagBody: Installation: all NagBody codes
# make_all/install_sh/install_all_additional_libs.sh 2>&1 |tee make_all/install_all.log
################################################################################

make -f NagBody install_cfitsio
make -f NagBody install_fftw2
make -f NagBody install_fftw3
make -f NagBody install_gsl
make -f NagBody install_Healpix
make -f NagBody install_lapack
make -f NagBody install_OpenBLAS
make -f NagBody install_plplot

