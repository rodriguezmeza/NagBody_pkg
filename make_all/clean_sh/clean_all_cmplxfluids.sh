################################################################################
#NagBody: Installation: all NagBody codes
# make_all/clean_sh/clean_all_cmplxfluids.sh 2>&1 |tee make_all/clean_all.log
################################################################################

make -f NagBody clean_analysis_md
make -f NagBody clean_md_blj
make -f NagBody clean_md_blj_n2
make -f NagBody clean_md_ic_model
make -f NagBody clean_md_lj_tree
