# ~/nbody/library/Makefile: build and test NBody library package.
# Copyright (c) 1997  Joshua E. Barnes  Honolulu, Hawai'i.

OBJFILES = bodyio_f3.o bodyio_d3.o buildmap.o phatbody.o\
	   snapcenter_f3.o snapcenter_d3.o snapkinem_f3.o snapkinem_d3.o

BINFILES = buildmap

ZENOPATH = $(HOME)/NagBody/NagBody_sources/General_libs/zeno
CCFLAGS = -DMACOSX -I$(ZENOPATH)/inc

# Assemble object library.

libNBody.a: $(OBJFILES)
	ar ruv libNBody.a $(OBJFILES)

tidy:
	rm -f $(OBJFILES) snapmap_s3.o snapmap_m3.o snapmap_d3.o

# Create individual object files.

bodyio_f3.o: bodyio.c phatbody.h
	$(CC) $(CCFLAGS) -DSINGLEPREC -DNDIM=3 -o bodyio_f3.o -c bodyio.c

bodyio_d3.o: bodyio.c phatbody.h
	$(CC) $(CCFLAGS) -DDOUBLEPREC -DNDIM=3 -o bodyio_d3.o -c bodyio.c

buildmap.o: buildmap.c mapdefs.h
	$(CC) $(CCFLAGS) -c buildmap.c

phatbody.o: phatbody.c phatbody.h
	$(CC) $(CCFLAGS) -c phatbody.c

snapcenter_f3.o: snapcenter.c phatbody.h
	$(CC) $(CCFLAGS) -DSINGLEPREC -DNDIM=3 -o snapcenter_f3.o \
	      -c snapcenter.c

snapcenter_d3.o: snapcenter.c phatbody.h
	$(CC) $(CCFLAGS) -DDOUBLEPREC -DNDIM=3 -o snapcenter_d3.o \
	      -c snapcenter.c

snapkinem_f3.o: snapkinem.c phatbody.h
	$(CC) $(CCFLAGS) -DSINGLEPREC -DNDIM=3 -o snapkinem_f3.o \
	      -c snapkinem.c

snapkinem_d3.o: snapkinem.c phatbody.h
	$(CC) $(CCFLAGS) -DDOUBLEPREC -DNDIM=3 -o snapkinem_d3.o \
	      -c snapkinem.c

# Main program object files for snapshot mapping.

snapmap_s3.o: snapmap.c phatbody.h
	$(CC) $(CCFLAGS) -DSINGLEPREC -DNDIM=3 -o snapmap_s3.o -c snapmap.c

snapmap_m3.o: snapmap.c phatbody.h
	$(CC) $(CCFLAGS) -DMIXEDPREC  -DNDIM=3 -o snapmap_m3.o -c snapmap.c

snapmap_d3.o: snapmap.c phatbody.h
	$(CC) $(CCFLAGS) -DDOUBLEPREC -DNDIM=3 -o snapmap_d3.o -c snapmap.c

# Utility program for snapshot mapping.

buildmap: buildmap.c mapdefs.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED -o buildmap buildmap.c -lClib


.PHONY : clean
clean:
	rm -f *.o libNBody.a ../lib/libNBody.a


