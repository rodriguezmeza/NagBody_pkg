################################################################################
#NagBody: Installation: all NagBody codes
# make_all/install_sh/install_all_class.sh 2>&1 |tee make_all/install_all.log
################################################################################

make -f NagBody install_class
make -f NagBody install_hiclass
make -f NagBody install_montepython
make -f NagBody install_planck2018

