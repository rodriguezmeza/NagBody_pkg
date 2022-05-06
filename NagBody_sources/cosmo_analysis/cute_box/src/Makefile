####################################################
###          User definable stuff
#Monopole
NB_R = 128
R_MAX = 100.
#Logarithmic binning
N_LOGINT = 5
#Behavior options
DEFINEOPTIONS = -D_VERBOSE -D_LONGIDS
# Mar
#DEFINEOPTIONS += -D_HAVE_OMP #Comment this out if you don't have the OpenMP headers
#DEFINEOPTIONS += -D_LOGBIN
### End of user-definable stuff
####################################################

#Transform to inverse scales
I_R_MAX = $(shell echo "scale=5;1./$(R_MAX)" | bc)
LOG_R_MAX = $(shell echo "scale=9;l($(R_MAX))/l(10)" | bc -l)

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)
DEFINEFLAGSCPU += -DNB_R=$(NB_R) -DI_R_MAX=$(I_R_MAX) -DLOG_R_MAX=$(LOG_R_MAX)
DEFINEFLAGSCPU += -DN_LOGINT=$(N_LOGINT)

# COMPILER AND OPTIONS
# Mar
COMPCPU = gcc
OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)
#COMPCPU = clang
#OPTCPU = -Wall -O3 $(DEFINEFLAGSCPU)
#OPTCPU += -I/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include

#INCLUDES AND LIBRARIES
# Mar
#INCLUDECOM = -I./src
#SRC = ./src
SRC = ./
INCLUDECOM = -I$(SRC)
LIBCPU = -lm
#LIBCPU += -L/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib


#.o FILES
#CUTE
DEF = $(SRC)/define.o
COM = $(SRC)/common.o
PM = $(SRC)/pm.o
TREE = $(SRC)/tree.o
IO = $(SRC)/io.o
NEIGH = $(SRC)/neighbors.o
CORR = $(SRC)/correlator.o
MAIN = $(SRC)/main.c
OFILES = $(DEF) $(COM) $(PM) $(TREE) $(IO) $(NEIGH) $(CORR) $(MAIN)

#FINAL GOAL
EXE = cute_box

#RULES
default : $(EXE)
#RULE TO MAKE .o's FROM .c's
$(DEF) : $(SRC)/define.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(COM) : $(SRC)/common.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(CORR) : $(SRC)/correlator.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(TREE) : $(SRC)/tree.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(IO) : $(SRC)/io.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(NEIGH) : $(SRC)/neighbors.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(PM) : $(SRC)/pm.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES)
	$(COMPCPU) $(OPTCPU) $(OFILES) -o $(EXE) $(INCLUDECOM) $(LIBCPU)

#CLEANING RULES
clean :
	rm -f $(SRC)*.o $(EXE)

cleaner :
	rm -f $(SRC)*.o .$(SRC)*~ *~ CUTE_box
