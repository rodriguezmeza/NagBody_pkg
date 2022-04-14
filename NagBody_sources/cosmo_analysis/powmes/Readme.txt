
https://www.vlasix.org/index.php?n=Main.Powmes

########################################################################
PARA COMPILAR E INSTALAR powmes.0.2

SE NECESITA instalar OpenMPI y FFTW2. Vienen en Additional libs:

$(NAGBODY_ROOT_DIR)/NagBody_sources/Additional_libs/Gadget

Las anteriores librerias se empaqueten en $(NAGBODY_ROOT_DIR) usando:

make -f NagBody packing_gadget207_Addlibs

DespuŽs entre en el directorio: $(NAGBODY_ROOT_DIR)/Readmes/Additional_libs
Y lea los readmes correspondientes y siga las instrucciones.


POR DEFAULT EL COMPILADOR ES INTEL.
SE NECESITA FFTW2 no FFTW3, compilado con intel y sin sufijos.

La generacion de FFTW2 tiene un detalle particular ... ????
Genermos OpenMPI con gcc e ifort y despues generamos 
FFTW2 sin prefijos (double) y con gcc e ifort ... (despues de todo powmes es f90).
Lo anterior lo hicimos sin hacer 'mkdir build_dir ... cd build_dir' ya que 
esto genero un problema para encontrar un header...

La ligas que usamos en local fueron:

openmpi -> compilers/openmpi-1.4_gcc_ifort
fftw2 -> nbody/fftw-2.1.5_gcc_ifort_openmpi


cd powmes.0.2/bin

make

Nota: nos da el warning option '-mcmodel' not supported

Hicimos la prueba:

cd test
../bin/powmes powmes.config.template
diff powspec.dat powspec.dat.template

Las diferencias seran en las ultimas cifras ... lo cual esta bien ...

Finalmente para tenerlo accesible lo pusimos en local/bin ...

cp powmes $HOME/local/bin/.

Para limpiar se hace:

make clean

########################################################################


