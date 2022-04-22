################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_class.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_class
make -f NagBody packing_hiclass
make -f NagBody packing_montepython
make -f NagBody packing_planck2018

