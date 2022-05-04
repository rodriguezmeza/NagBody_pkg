
Dependencies:

fftw-3: see Readme-Install_fftw-3.txt to install it under $HOME/NagBody_pkg
gsl: see Readme-Install_gsl.txt to install it under $HOME/NagBody_pkg

To install:

cd $HOME/NagBody_pkg
make -f NagBody install_cola

Or install it in its own directory:

cd src
make
cd ..

Use it:

$ man ./doc/cola.1

to see the manual. Also visit:

https://bitbucket.org/tassev/colacode

Test it:

cd tests
../src/cola 


