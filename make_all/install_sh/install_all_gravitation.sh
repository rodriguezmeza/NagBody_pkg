################################################################################
#NagBody: Installation: all NagBody codes
# make_all/install_sh/install_all_gravitation.sh 2>&1 |tee make_all/install_all.log
################################################################################

make -f NagBody install_nbody_n2
make -f NagBody install_gbsph
make -f NagBody install_model
