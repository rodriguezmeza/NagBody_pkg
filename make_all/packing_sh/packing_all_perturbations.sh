################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_perturbations.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_mgpt
make -f NagBody packing_TreeCorr
