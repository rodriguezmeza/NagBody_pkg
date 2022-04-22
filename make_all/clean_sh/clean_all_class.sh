################################################################################
#NagBody: Installation: all NagBody codes
# make_all/clean_sh/clean_all_class.sh 2>&1 |tee make_all/clean_all.log
################################################################################

make -f NagBody clean_class
make -f NagBody clean_hiclass
make -f NagBody clean_montepython
make -f NagBody clean_planck2018

