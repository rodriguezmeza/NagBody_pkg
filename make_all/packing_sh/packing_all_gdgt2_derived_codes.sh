################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_gdgt2_derived_codes.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_galic
make -f NagBody packing_gdgt207cluster
make -f NagBody packing_gdgt207galaxy
make -f NagBody packing_gdgt207gassphere
make -f NagBody packing_gdgt207lcdm
