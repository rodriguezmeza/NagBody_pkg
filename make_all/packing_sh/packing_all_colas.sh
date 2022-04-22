################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_colas.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_cola
make -f NagBody packing_cola_simple
make -f NagBody packing_lpicola
make -f NagBody packing_mgpicola
make -f NagBody packing_mgpicola_fofr
make -f NagBody packing_mgpicola_pofk
make -f NagBody packing_mgpicola_pofk_mpi
make -f NagBody packing_mgpicola_pofk_rsd_mpi

