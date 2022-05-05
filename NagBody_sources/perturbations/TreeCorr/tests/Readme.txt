
It is best to leave this directory clean. Copy the whole tests_from-src directory

cp -pR tests_from-src ../runs

Then move there

cd ../runs

and execute:

time ./run_all_tests 2>&1 | tee run_all_tests.log


Also consider:

If you want to run the unit tests, you can do the following::

pip install -r test_requirements.txt
cd tests
nosetests
