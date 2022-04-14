

This is MontePython version 3.5.0

(Download a new version and substitute montepython_public with the new one)

(This is version downloaded 2022-04-05)

1. Add to $HOME/.bash_profile the line

Not working! ::
export PATH=/path/to/MontePython/montepython/:$PATH
export PATH=$HOME/NagBody_pkg/NagBody_sources/class/montepython/montepython_public/montepython:$PATH

to be able to call the program from anywhere.


2. Planck2018:

Edit $HOME/.bash_profile and add the line:

source $HOME/NagBody_pkg/NagBody_sources/class/planck/Planck2018/code/plc_3.0/plc-3.01/bin/clik_profile.sh

In directory $HOME/NagBody_pkg/NagBody_sources/class/montepython/montepython_public/data make the link

ln -s $HOME/NagBody_pkg/NagBody_sources/class/planck/Planck2018/baseline/plc_3.0 clik_14.0


3. In directory $HOME/NagBody_pkg/NagBody_sources/class/montepython/runs make the links:

ln -s ../montepython_public/data data
ln -s ../montepython_public/covmat/ covmat
ln -s ../montepython_public/bestfit/ bestfit

Can also be useful to make a link to montepython_public/montepython folder

ln -s ../montepython_public/montepython montepython


4. Finally edit default.conf in order to have class and Planck2018 paths defined properly. This file is used by default. If you have other config files they can be called bye montepython using the flat "--conf".

5. The runs/input folder is just a copy of the one in montepython_public folder. It has several input parameters files.

Examples:

python montepython/MontePython.py run --conf default.conf  -o results/jla -p input/jla.param -N 20000
python montepython/MontePython.py info results/jla --bins=10

6. Install getdist to have better confidence plots.

From directory runs

python ../python/GetDist.py getdist/distparams_Union2.ini


%%%%%%

See Readme.txt in tests or in runs if you make this copy.

%%%%%%
If you modify class could be better to copy src directory as

cp -pR src/ src_MODELNAME

where MODELNAME is a string characterizing the new class model.

Remember to refresh classy by:

cd $HOME/NagBody_pkg
make -f NagBody clean_class
cd $HOME/NagBody_pkg/NagBody_sources/class/class/src_MODELNAME

make clean
make all
cd python
python setup.py build
python setup.py install --user
cd ..


