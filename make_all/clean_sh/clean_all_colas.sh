################################################################################
#NagBody: Installation: all NagBody codes
# make_all/clean_sh/clean_all_colas.sh 2>&1 |tee make_all/clean_all.log
################################################################################

make -f NagBody clean_cola
make -f NagBody clean_cola_simple
make -f NagBody clean_lpicola
make -f NagBody clean_mgpicola
make -f NagBody clean_mgpicola_fofr
make -f NagBody clean_mgpicola_pofk
make -f NagBody clean_mgpicola_pofk_mpi
make -f NagBody clean_mgpicola_pofk_rsd_mpi

