
NORMAL NAGBODY and GNU are WORKING

This is CLASS version dowloaded 2022-02-06

Produce with Anaconda 2. With Anaconda 3 is now working.

Go to src/class/class/python dir and execute:

$ python setup.py build
$ python setup.py install --user

in order to produce the python wrapper for class (classy).


%%%%%%
To test class copy tests directory as run

cp -pR tests/ runs

cd runs
class default.ini

Run Jupyter:

jupyter notebook

and open notebooks/warmup.ipynb 

Execute all cells.

%%%%%%
To test class in montepython 

cd $HOME/NagBody_pkg/NagBody_sources/class/montepython
cp -pR tests/ runs
cd runs

Download JLA data:

The JLA data files were not found. Please download the following link
http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v6.tgz, extract
it, and copy all files present in `jla_likelihood_v6/data` to
`your_montepython/data/JLA`

cp jla_likelihood_v6/data/* ./data/JLA/.


Edit default.conf properly and run

python montepython/MontePython.py run --conf default.conf  -o results/jla -p input/jla.param -N 20000
python montepython/MontePython.py info results/jla --bins=10



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

%%%%%%
Check which classy are in site-packages:

ls -lht /opt/anaconda3/lib/python3.9/site-packages/classy*
ls -lht /Users/mar/.local/lib/python3.9/site-packages/

rm /opt/anaconda3/lib/python3.9/site-packages/classy*
rm /Users/mar/.local/lib/python3.9/site-packages/classy*

%%%%%%
Remove classy from site-packages:

rm /Users/mar/.local/lib/python3.9/site-packages/classy*
rm /opt/anaconda3/lib/python3.9/site-packages/classy*
