################################################################################
#NagBody: Installation: all NagBody codes
# make_all/install_sh/install_all_cmplxfluids.sh 2>&1 |tee make_all/install_all.log
################################################################################

make -f NagBody install_analysis_md
make -f NagBody install_md_blj
make -f NagBody install_md_blj_n2
make -f NagBody install_md_ic_model
make -f NagBody install_md_lj_tree
