################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_gravitation.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_analysis_galaxy
make -f NagBody packing_analysis_grav
make -f NagBody packing_galaxy_hernquist
make -f NagBody packing_galaxy_model
make -f NagBody packing_galaxy_starscream
make -f NagBody packing_gbsph
make -f NagBody packing_model
make -f NagBody packing_nbody_n2
