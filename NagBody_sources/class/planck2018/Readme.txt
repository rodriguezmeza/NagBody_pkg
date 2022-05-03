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


#####
Problems using Homebrew in Mac:

FileNotFoundError: [Errno 2] No such file or directory: '/usr/local/bin/cython'

Or using make:

make clean
make install
make install_python

Traceback (most recent call last):
  File "/Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/class/planck2018/src/code/plc_3.0/plc-3.01/setup.py", line 3, in <module>
    from Cython.Distutils import build_ext
ModuleNotFoundError: No module named 'Cython'
make: *** [install_python] Error 1

Was solved

$ brew install cython

==> Caveats
cython is keg-only, which means it was not symlinked into /usr/local,
because this formula is mainly used internally by other formulae.
Users are advised to use `pip` to install cython.

If you need to have cython first in your PATH, run:
  echo 'export PATH="/usr/local/opt/cython/bin:$PATH"' >> /Users/mar/.bash_profile

==> Summary
/usr/local/Cellar/cython/0.29.28: 325 files, 7.0MB
==> Running `brew cleanup cython`...
Disable this behaviour by setting HOMEBREW_NO_INSTALL_CLEANUP.
Hide these hints with HOMEBREW_NO_ENV_HINTS (see `man brew`).

AND then executing

(base) lapmar:plc-3.01 mar$ echo 'export PATH="/usr/local/opt/cython/bin:$PATH"' >> /Users/mar/.bash_profile

LAST History:

  557  conda uninstall Cython
  558  make clean
  559  ./waf configure
  560  history |grep ./waf
  561  ./waf --cfitsio_include=$HOME/NagBody_pkg/local/cfitsio/include --cfitsio_lib=$HOME/NagBody_pkg/local/cfitsio/lib configure
  562  brew uninstall cython
  563  conda install Cython
  564  ./waf --cfitsio_include=$HOME/NagBody_pkg/local/cfitsio/include --cfitsio_lib=$HOME/NagBody_pkg/local/cfitsio/lib configure
  565  ./waf install
  566  make clean
  567  make
  568  make install
  569  make install_python
  570  brew install cython
  571  echo 'export PATH="/usr/local/opt/cython/bin:$PATH"' >> /Users/mar/.bash_profile
  572  history



