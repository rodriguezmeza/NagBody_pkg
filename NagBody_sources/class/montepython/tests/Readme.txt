
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

(Run twice to have better samplings::)
python montepython/MontePython.py run --conf default.conf  -o results/jla -p input/jla.param -N 20000
python montepython/MontePython.py run --conf default.conf  -o results/jla -p input/jla.param -N 20000
python montepython/MontePython.py info results/jla --bins=10

cd results

cp -pR chains chains_jla

Copy both chains into chains_jla directory as chain_JLA_1.txt and chain_JLA_2.txt
Copy also the ???.paramnames to chains_jla directory.


cp jla/2022-04-08_20000_.paramnames chains_jla/chain_JLA.paramnames
cp jla/2022-04-08_20000__1.txt chains_jla/chain_JLA_1.txt
cp jla/2022-04-08_20000__2.txt chains_jla/chain_JLA_2.txt

cd chains_jla

Edit chain_JLA.paramnames if it is necessary to have math symbols in LaTeX, then
edit distparams.in 

file_root = ./chain_JLA
no_plots = F

and execute:

getdist distparams.ini
python chain_jla.py
python chain_jla_tri.py

Finally execute

jupyter notebook

and load getdist-plots_stat.ipynb

In directories confidenceRegions and posteriors you will have two pdf files.



%%%%%%
See also:

../xother_files/From_cosmomc_camb_purgar-y-mantener to find other python scripts and readme.

Getdist guidelines are in:

https://getdist.readthedocs.io/en/latest/index.html


%%%%%%
To see how many classy packs are in python do:

ls -lht /opt/anaconda3/lib/python3.9/site-packages/classy*

Could be useful activate an environment.

%%%%%%
For Planck 2015 (updated to python3, see # Mar) make this link:

ln -s /Users/mar/NagBody_pkg/NagBody_sources/class/planck2015/src/plc_2.0 clik

For Planck 2018 make this link:

ln -s /Users/mar/NagBody_pkg/NagBody_sources/class/planck2018/src/baseline/plc_3.0/ clik_14.0

%%%%%%
Other tests using Input files

(
python montepython/MontePython.py run --help
python montepython/MontePython.py run -h [short version]
)

(Planck2018 input files gives Segmentation fault)

python montepython/MontePython.py run --conf default.conf -o results/test_gaussian -p input/test_gaussian.param -N 20000
python montepython/MontePython.py info results/test_gaussian --bins=10

python montepython/MontePython.py run --conf default.conf -o results/test -p input/test.param -N 20000
python montepython/MontePython.py info results/test --bins=10

python montepython/MontePython.py run --conf default.conf -o results/test -p input/test_is.param -N 20000
python montepython/MontePython.py info results/test_is --bins=10

(Gives Segmentation fault)
python montepython/MontePython.py run --conf default.conf -o results/sdss_lrgDR7 -p input/sdss_lrgDR7.param -N 20000
python montepython/MontePython.py info results/sdss_lrgDR7 --bins=10

(Gives Segmentation fault: 11)
python montepython/MontePython.py run --conf default.conf -o results/base2018TTTEEE -p input/base2018TTTEEE.param -N 20000 --covmat covmat/base2018TTTEEE.covmat

python montepython/MontePython.py run --conf default_planck2015.conf -o results/base2018TTTEEE -p input/base2018TTTEEE.param -N 20000

python montepython/MontePython.py run --conf default.conf -o results/base2018TTTEEE_lite -p input/base2018TTTEEE_lite.param -N 20000 --bestfit bestfit/base2018TTTEEE_lite.bestfit --covmat covmat/base2018TTTEEE_lite.covmat --j fast --method NS

python montepython/MontePython.py run --conf default.conf -o results/base2018TTTEEE_lite -p input/base2018TTTEEE_lite.param -N 20000 --bestfit bestfit/base2018TTTEEE_lite.bestfit --covmat covmat/base2018TTTEEE_lite.covmat --j fast --method Fisher

(
-j global, sequential, fast
--jumping global
)
