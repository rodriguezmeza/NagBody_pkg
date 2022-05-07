################################################################################
#NagBody: Installation: all NagBody codes
# make_all/install_sh/install_all_gdgt2_derived_codes.sh 2>&1 |tee make_all/install_all.log
################################################################################

make -f NagBody install_galic
make -f NagBody install_gdgt207cluster
make -f NagBody install_gdgt207galaxy
make -f NagBody install_gdgt207gassphere
make -f NagBody install_gdgt207lcdm
