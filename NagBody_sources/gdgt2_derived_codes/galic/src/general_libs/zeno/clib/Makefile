# ~/lib/clib/Makefile: build and test Zeno C library package.
# Copyright (c) 1999  Joshua E. Barnes  Honolulu, Hawai'i.

OBJFILES = allocate.o besselfunc.o burststring.o cputime.o datatypes.o\
	   error.o extstring.o filestruct.o getparam.o history.o\
	   mathfns_d.o mathfns_f.o phatstruct.o pickpnt_d.o pickpnt_f.o\
	   pickpnt_m.o scanopt.o setrange_d.o setrange_f.o setrange_m.o\
	   spline_d.o spline_f.o spline_m.o stropen.o strset.o\
	   vectmath_d.o vectmath_f.o vectmath_m.o within.o xrandom.o

TESTCODE = bessel_test burststring_test datatypes_test\
	   error_test extstring_test filestruct_test getparam_test\
	   mathfns_test_d mathfns_test_f mathfns_test_m phatstruct_test\
	   scanopt_test setrange_test spline_test_d spline_test_f\
	   spline_test_m stropen_test strset_test vectmath_test_d\
	   vectmath_test_f vectmath_test_m within_test

ZENOPATH = $(HOME)/NagBody/NagBody_sources/General_libs/zeno
CCFLAGS = -DMACOSX -I$(ZENOPATH)/inc

# Assemble object library.

libClib.a: $(OBJFILES)
	ar ruv libClib.a $(OBJFILES)

tidy:
	rm -f $(OBJFILES) $(TESTCODE)

# Create individual object files.

allocate.o: allocate.c stdinc.h
	$(CC) $(CCFLAGS) -c allocate.c

besselfunc.o: besselfunc.c stdinc.h
	$(CC) $(CCFLAGS) -c besselfunc.c

burststring.o: burststring.c stdinc.h
	$(CC) $(CCFLAGS) -c burststring.c

cputime.o: cputime.c stdinc.h
	$(CC) $(CCFLAGS) -c cputime.c

datatypes.o: datatypes.c datatypes.h stdinc.h
	$(CC) $(CCFLAGS) -c datatypes.c

error.o: error.c stdinc.h
	$(CC) $(CCFLAGS) -c error.c

extstring.o: extstring.c stdinc.h
	$(CC) $(CCFLAGS) -c extstring.c

filestruct.o: filestruct.c filestruct.h stdinc.h
	$(CC) $(CCFLAGS) -c filestruct.c

getparam.o: getparam.c getparam.h stdinc.h
	$(CC) $(CCFLAGS) -c getparam.c

history.o: history.c stdinc.h
	$(CC) $(CCFLAGS) -c history.c

mathfns_d.o: mathfns.h mathfns.c stdinc.h
	$(CC) $(CCFLAGS) -c mathfns.c -DDOUBLEPREC -o mathfns_d.o

mathfns_f.o: mathfns.h mathfns.c stdinc.h
	$(CC) $(CCFLAGS) -c mathfns.c -DSINGLEPREC -o mathfns_f.o

phatstruct.o: phatstruct.c phatstruct.h stdinc.h
	$(CC) $(CCFLAGS) -c phatstruct.c

pickpnt_d.o: pickpnt.c stdinc.h
	$(CC) $(CCFLAGS) -c pickpnt.c -DDOUBLEPREC -o pickpnt_d.o

pickpnt_f.o: pickpnt.c stdinc.h
	$(CC) $(CCFLAGS) -c pickpnt.c -DSINGLEPREC -o pickpnt_f.o

pickpnt_m.o: pickpnt.c stdinc.h
	$(CC) $(CCFLAGS) -c pickpnt.c -DMIXEDPREC  -o pickpnt_m.o

scanopt.o: scanopt.c stdinc.h
	$(CC) $(CCFLAGS) -c scanopt.c

setrange_d.o: setrange.c stdinc.h
	$(CC) $(CCFLAGS) -c setrange.c -DDOUBLEPREC -o setrange_d.o

setrange_f.o: setrange.c stdinc.h
	$(CC) $(CCFLAGS) -c setrange.c -DSINGLEPREC -o setrange_f.o

setrange_m.o: setrange.c stdinc.h
	$(CC) $(CCFLAGS) -c setrange.c -DMIXEDPREC  -o setrange_m.o

spline_d.o: spline.c stdinc.h
	$(CC) $(CCFLAGS) -c spline.c -DDOUBLEPREC -o spline_d.o

spline_f.o: spline.c stdinc.h
	$(CC) $(CCFLAGS) -c spline.c -DSINGLEPREC -o spline_f.o

spline_m.o: spline.c stdinc.h
	$(CC) $(CCFLAGS) -c spline.c -DMIXEDPREC  -o spline_m.o

stropen.o: stropen.c stdinc.h
	$(CC) $(CCFLAGS) -c stropen.c

strset.o: strset.c strset.h stdinc.h
	$(CC) $(CCFLAGS) -c strset.c

vectmath_d.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(CC) $(CCFLAGS) -c vectmath.c -DDOUBLEPREC -o vectmath_d.o

vectmath_f.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(CC) $(CCFLAGS) -c vectmath.c -DSINGLEPREC -o vectmath_f.o

vectmath_m.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(CC) $(CCFLAGS) -c vectmath.c -DMIXEDPREC  -o vectmath_m.o

within.o: within.c stdinc.h
	$(CC) $(CCFLAGS) -c within.c

xrandom.o: xrandom.c stdinc.h
	$(CC) $(CCFLAGS) -c xrandom.c

# Test library routines.

tests: libClib.a
	$(MAKE) $(TESTCODE)

bessel_test: besselfunc.c stdinc.h
	$(CC) $(CCFLAGS) -o bessel_test -DTESTBED besselfunc.c libClib.a -lm

burststring_test: burststring.c stdinc.h
	$(CC) $(CCFLAGS) -o burststring_test -DTESTBED burststring.c libClib.a

datatypes_test: datatypes.h datatypes.c stdinc.h
	$(CC) $(CCFLAGS) -o datatypes_test -DTESTBED datatypes.c libClib.a

error_test: error.c stdinc.h
	$(CC) $(CCFLAGS) -o error_test -DTESTBED error.c

extstring_test: extstring.c stdinc.h
	$(CC) $(CCFLAGS) -o extstring_test -DTESTBED extstring.c libClib.a

filestruct_test: filestruct_test.c stdinc.h
	$(CC) $(CCFLAGS) -o filestruct_test filestruct_test.c libClib.a

getparam_test: getparam.c stdinc.h
	$(CC) $(CCFLAGS) -o getparam_test -DTESTBED getparam.c libClib.a

mathfns_test_d: mathfns.c mathfns.h stdinc.h
	$(CC) $(CCFLAGS) -o mathfns_test_d -DDOUBLEPREC -DTESTBED mathfns.c \
	   libClib.a -lm

mathfns_test_f: mathfns.c mathfns.h stdinc.h
	$(CC) $(CCFLAGS) -o mathfns_test_f -DSINGLEPREC -DTESTBED mathfns.c \
	   libClib.a -lm

mathfns_test_m: mathfns.c mathfns.h stdinc.h
	$(CC) $(CCFLAGS) -o mathfns_test_m -DMIXEDPREC  -DTESTBED mathfns.c \
	   libClib.a -lm

phatstruct_test: phatstruct.c phatstruct.h stdinc.h
	$(CC) $(CCFLAGS) -o phatstruct_test -DTESTBED phatstruct.c libClib.a

scanopt_test: scanopt.c stdinc.h
	$(CC) $(CCFLAGS) -o scanopt_test -DTESTBED scanopt.c libClib.a

setrange_test: setrange.c stdinc.h
	$(CC) $(CCFLAGS) -o setrange_test -DTESTBED setrange.c libClib.a -lm

spline_test_d: spline.c stdinc.h
	$(CC) $(CCFLAGS) -o spline_test_d -DDOUBLEPREC -DTESTBED spline.c \
	   libClib.a -lm

spline_test_f: spline.c stdinc.h
	$(CC) $(CCFLAGS) -o spline_test_f -DSINGLEPREC -DTESTBED spline.c \
	   libClib.a -lm

spline_test_m: spline.c stdinc.h
	$(CC) $(CCFLAGS) -o spline_test_m -DMIXEDPREC  -DTESTBED spline.c \
	   libClib.a -lm

stropen_test: stropen.c stdinc.h
	$(CC) $(CCFLAGS) -o stropen_test -DTESTBED stropen.c libClib.a

strset_test: strset.c strset.h stdinc.h
	$(CC) $(CCFLAGS) -o strset_test -DTESTBED strset.c libClib.a

vectmath_test_d: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(CC) $(CCFLAGS) -o vectmath_test_d -DDOUBLEPREC vectmath_test.c \
	   libClib.a -lm

vectmath_test_f: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(CC) $(CCFLAGS) -o vectmath_test_f -DSINGLEPREC vectmath_test.c \
	   libClib.a -lm

vectmath_test_m: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(CC) $(CCFLAGS) -o vectmath_test_m -DMIXEDPREC  vectmath_test.c \
	   libClib.a -lm

within_test: within.c stdinc.h
	$(CC) $(CCFLAGS) -o within_test -DTESTBED within.c libClib.a
	

.PHONY : clean
clean:
	rm -f *.o libClib.a ../lib/libClib.a


