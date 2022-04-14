This is Planck data release 2018

In directory $HOME/NagBody_pkg run:

make -f NagBody install_planck2018

To use with class/montepython run:

source $HOME/NagBody_pkg/bin/clik_profile.sh

Include above line in .bash_profile to have clik environment always available. Check the environment:

env |grep CLIK

Also make the link in montepython/data

cd $HOME/NagBody_pkg/NagBody_sources/class/montepython/montepython_public/data
ln -s $HOME/NagBody_pkg/NagBody_sources/class/planck2018/data/baseline/plc_3.0/ clik_14.0



%%%%%
Installation using waf:

./waf --cfitsio_include=$HOME/NagBody_pkg/local/cfitsio/include --cfitsio_lib=$HOME/NagBody_pkg/local/cfitsio/lib configure

./waf install

You may uninstall or clean:

./waf uninstall
./waf clean

(Could be useful to clean thoroughly:)
make clean


(Other options of waf:)
./waf distclean


source /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/code/plc_3.0/plc-3.01/bin/clik_profile.sh

or

source $HOME/NagBody_pkg/bin/clik_profile.sh


%%%%%
Tests:

./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22_TT.clik/
./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik/
./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/hi_l/plik_lite/plik_lite_v22_TT.clik/
./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik/

./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_CMBmarged.clik_lensing/
./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing/

./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_BB_Aplanck_B.clik/
./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik/
./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EEBB_Aplanck_B.clik/

./clik_example_py /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik/



