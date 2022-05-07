################################################################################
#NagBody: Installation: all NagBody codes
# make_all/clean_sh/clean_all_gravitation.sh 2>&1 |tee make_all/clean_all.log
################################################################################

make -f NagBody clean_analysis_galaxy
make -f NagBody clean_analysis_grav
make -f NagBody clean_galaxy_hernquist
make -f NagBody clean_galaxy_model
make -f NagBody clean_galaxy_starscream
make -f NagBody clean_gbsph
make -f NagBody clean_model
make -f NagBody clean_nbody_n2
