################################################################################
#NagBody: Packing: all NagBody codes
#
# IT IS RECOMENDED TO CLEAN NagBody codes BEFORE PACKING!
#
# make_all/packing_sh/packing_all_cmplxfluids.sh 2>&1 |tee make_all/packing_all.log
################################################################################

make -f NagBody packing_analysis_md
make -f NagBody packing_md_blj
make -f NagBody packing_md_blj_n2
make -f NagBody packing_md_ic_model
make -f NagBody packing_md_lj_tree
