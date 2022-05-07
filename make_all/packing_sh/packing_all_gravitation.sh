################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_gravitation.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_nbody_n2
make -f NagBody packing_gbsph
make -f NagBody packing_model
