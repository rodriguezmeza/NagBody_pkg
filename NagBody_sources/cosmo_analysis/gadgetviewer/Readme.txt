See the latest reales:
https://github.com/jchelly/gadgetviewer/releases

Install GTK+

sudo port install gtk3 

(
  dbus has the following notes:
    Startup items (named 'dbus-system, dbus-session') have been generated that will aid in starting dbus with launchd. They are disabled by default. Execute the following
    command to start them, and to cause them to launch at startup:
    
        sudo port load dbus
)

No funciona gtk3

SE NECESITA instalar GTK2+:

$ sudo port install gtk2

THEN GO TO gadgetviewer_2022-03-30

autoreconf -i
autoscan
./configure --prefix=$HOME/NagBody_pkg/local/gadgetviewer
make
make install

We got error at "make"

make distclean
sudo port install hdf5



%%%%%%%%%%%%%%%%%%%%%%%%%

òltima versi—n instalada desde:

xother_files/gadgetviewer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FINAL VERSION Working

SE NECESITA instalar GTK2+:

$ sudo port install gtk2


$ autoreconf -i
$ autoscan
$ export CC=clang
$ ./configure --prefix=$HOME/NagBody_pkg/local/gadgetviewer FCFLAGS=-I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include LDFLAGS="-L/usr/lib -F/Library/Frameworks -F/System/Library/Frameworks"
$ make
$ make install

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



########################################################################
PARA COMPILAR E INSTALAR 

SE NECESITA instalar GTK2+:

$ sudo port install gtk2

$ autoreconf -i
$ autoscan
$ ./configure --prefix=$HOME/NagBody_pkg/local/gadgetviewer
$ make


############################
(cc1: error: argument to '-O' should be a non-negative integer, 'g', 's' or 'fast'

al hacer :
libtool: compile: tick_interval.c

)

Funciona (evitando PLPLOT_FLAGS):
./configure --prefix=$HOME/NagBody_pkg/local/gadgetviewer --without-plplot 2>&1 | tee configure_gcc.log
make 2>&1 | tee make_gcc.log
make install
make distclean


########################################################################

Possibilities to solve the problem:
ld: framework not found CoreFoundation


 ./configure FCFLAGS=-I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include LDFLAGS=-L/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib

Not functioning


-F${frameworks.CoreFoundation}/Library/Frameworks -framework CoreFoundation $NIX_LDFLAGS

This does the work!

"-L/usr/lib -F/Library/Frameworks -F/System/Library/Frameworks" ::

./configure --prefix=$HOME/NagBody_pkg/local/gadgetviewer --with-plplo=/opt/local/lib --with-hdf5=/opt/local/lib FCFLAGS=-I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include LDFLAGS="-L/usr/lib -F/Library/Frameworks -F/System/Library/Frameworks"

hdf5 and plplot not working
