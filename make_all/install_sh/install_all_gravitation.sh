################################################################################
#NagBody: Installation: all NagBody codes
# make_all/install_sh/install_all_gravitation.sh 2>&1 |tee make_all/install_all.log
################################################################################

make -f NagBody install_analysis_galaxy
make -f NagBody install_analysis_grav
make -f NagBody install_galaxy_hernquist
make -f NagBody install_galaxy_model
make -f NagBody install_galaxy_starscream
make -f NagBody install_gbsph
make -f NagBody install_model
make -f NagBody install_nbody_n2
