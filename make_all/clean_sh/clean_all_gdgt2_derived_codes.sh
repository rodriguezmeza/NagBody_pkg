################################################################################
#NagBody: Installation: all NagBody codes
# make_all/clean_sh/clean_all_gravitation.sh 2>&1 |tee make_all/clean_all.log
################################################################################

make -f NagBody clean_galic
make -f NagBody clean_gdgt207cluster
make -f NagBody clean_gdgt207galaxy
make -f NagBody clean_gdgt207gassphere
make -f NagBody clean_gdgt207lcdm
