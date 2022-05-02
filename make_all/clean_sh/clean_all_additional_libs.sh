################################################################################
#NagBody: Installation: all NagBody codes
# make_all/clean_sh/clean_all_additional_libs.sh 2>&1 |tee make_all/clean_all.log
################################################################################

make -f NagBody clean_cfitsio
make -f NagBody clean_fftw2
make -f NagBody clean_fftw3
make -f NagBody clean_gsl
make -f NagBody clean_Healpix
make -f NagBody clean_lapack
make -f NagBody clean_OpenBLAS
make -f NagBody clean_plplot

