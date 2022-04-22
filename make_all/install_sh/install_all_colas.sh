################################################################################
#NagBody: Installation: all NagBody codes
# make_all/install_sh/install_all_colas.sh 2>&1 |tee make_all/install_all.log
################################################################################

make -f NagBody install_cola
make -f NagBody install_cola_simple
make -f NagBody install_lpicola
make -f NagBody install_mgpicola
make -f NagBody install_mgpicola_fofr
make -f NagBody install_mgpicola_pofk
make -f NagBody install_mgpicola_pofk_mpi
make -f NagBody install_mgpicola_pofk_rsd_mpi

