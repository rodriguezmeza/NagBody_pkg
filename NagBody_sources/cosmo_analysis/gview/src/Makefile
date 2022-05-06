#CPP=g++
#CPP=g++-mp-4.6
CPP=g++

VERSION       = 0.8

#INCLUDE       = # -I/usr/local/include
#INCLUDE       = -I/opt/local/include
INCLUDE       = -I/usr/X11/include
WARNINGS      = -Wall -Winit-self -pedantic
OPTIMIZATIONS = -O3

OPTIONS = -DVERSION=$(VERSION) $(WARNINGS) $(OPTIMIZATIONS) $(INCLUDE) -DNDEBUG # -g

DEP=$(shell ls *.h) Makefile
SRC=$(shell ls *.cpp)
OBJ=$(SRC:.cpp=.oo)

BIN=gview

$(BIN): $(OBJ) $(DEP)
	$(CPP) $(OPTIONS) -L/usr/X11/lib -lGL -lGLU -lglut $(OBJ) -o $@


#	$(CPP) $(OPTIONS) -L/opt/local/lib -lglut $(OBJ) -o $@

#	$(CPP) $(OPTIONS) -L/lib -lglut $(OBJ) -o $@

%.oo: %.cpp $(DEP)
	@ echo $(CPP) $<
	@ $(CPP) $(OPTIONS) -c $< -o $@

.PHONY:sync
sync:
	@ cp ../Speciale/array.h ../Speciale/configfile.h ../Speciale/exception.h ../Speciale/triple.h ../Speciale/ver.h .

.PHONY:install
install: $(BIN)
	@ echo copying $+ to install dir ~/Bin/cosmo_local
	@ cp $+ ~/Bin/cosmo_local # --attribute=-

.PHONY:pub
pub: clean sync install $(BIN)
	@ echo making publish version
	@ rm -rf $(OBJ) $(BIN) *~
	@ cp ../../LaTeX/Extra/gview.1 .
	@ touch *
	@ pushd .. ; tar czf gadget_addon_gview-$(VERSION).tgz gadget_addon_gview ; mv gadget_addon_gview-$(VERSION).tgz pub ; popd

.PHONY:clean
clean:
	@ rm -rf $(OBJ) $(BIN)
