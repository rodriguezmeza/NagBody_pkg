# this code cannot be run directly
# do 'source env_config/nagbodyrc.sh' from your sh shell or put it in your profile
#
#

#
#To solve some configuration problems related to zsh visit:
#https://scriptingosx.com/2019/06/moving-to-zsh-part-2-configuration-files/
#
#
# NagBody Defintion of variables ...

# To activate nnbodykit environment (added to load notebooks in jupyter, does not work)
#source /anaconda2/bin/activate /anaconda2/envs/nbodykit-env/

# NagBody:
export PATH=${HOME}/NagBody_pkg/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/man:${MANPATH}

# CFITSIO (NEEDED BY CAMB):
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib:${DYLD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib:${LD_LIBRARY_PATH}

# FFTW-2:
# INSTALADO CON PORT ... fftw, fftw-single, falto que se incluya con openmpi ...
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw2/lib:${DYLD_LIBRARY_PATH}

# FFTW3:
# Instalar con port ... fftw-3 ... tambien la precision simple fftw-3-single
export PATH=${HOME}/NagBody_pkg/local/fftw3/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/fftw3/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/fftw3/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw3/lib:${DYLD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw3/lib:${LD_LIBRARY_PATH}

# GSL:
# INSTALADO CON PORT ...
export PATH=${HOME}/NagBody_pkg/local/gsl/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/gsl/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/gsl/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/gsl/lib:${LD_LIBRARY_PATH}
#

# LAPACK (Instalada con ports en Mac):
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/lapack/lib:${LD_LIBRARY_PATH}

###################################################################
# PLPlot (ver. 5.15.0):
export PATH=${HOME}/NagBody_pkg/local/plplot/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/plplot/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/plplot/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib/plplot5.15.0/driversd:${DYLD_LIBRARY_PATH}

###################################################################
# HDF5:
# INSTALADO CON PORT ... hdf5-18 ...
#export PATH=${HOME}/NagBody_pkg/local/hdf5/bin:${PATH}
#export MANPATH=${HOME}/NagBody_pkg/local/hdf5/share/man:${MANPATH}
#export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/hdf5/lib:${LD_LIBRARY_PATH}
#export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/hdf5/lib:${DYLD_LIBRARY_PATH}

# PHDF5 (Parallel):
#export PATH=${HOME}/NagBody_pkg/local/phdf5/bin:${PATH}
#export MANPATH=${HOME}/NagBody_pkg/local/phdf5/share/man:${MANPATH}
#export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/phdf5/lib:${LD_LIBRARY_PATH}
#export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/phdf5/lib:${DYLD_LIBRARY_PATH}

# GADGETVIEWER:
export PATH=${HOME}/NagBody_pkg/local/gadgetviewer/bin:${PATH}

###################################################################
# Planck 2018:
# To activate Planck 2018 for CLASS:
#source $HOME/NagBody_pkg/NagBody_sources/class/planck2018/src/code/plc_3.0/plc-3.01/bin/clik_profile.sh
source $HOME/NagBody_pkg/bin/clik_profile.sh



###################################################################
# Python environment (ver 3.9):
#export PYTHONPATH=${HOME}/.local/lib/python3.9/site-packages:${PYTHONPATH}

