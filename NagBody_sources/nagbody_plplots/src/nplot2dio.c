/* =============================================================================
	MODULE: nplot2d_io.c			[nplot2d]
	Written by: Mario A. Rodriguez-Meza
	Starting date:	May 2006
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#include "globaldefs.h"


void StartOutput(void)
{
	gd.comment = "nplot2d comment about a particular running";
    printf("\n  \t -- %s --\n", gd.comment);
}

#define IN 1
#define OUT 0
#define SI 1
#define NO 0

void InputData(string filename, int col1, int col2, int *npts)
{
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
	short int *lineQ;

    instr = stropen(filename, "r");

	fprintf(stdout,
		"\nReading columns %d and %d from file %s... ",col1,col2,filename);

	state = OUT;
	nl = nw = nc = 0;
	while ((c = getc(instr)) != EOF) {
		++nc;
		if (c=='\n')
			++nl;
		if (c==' ' || c=='\n' || c=='\t')
			state = OUT;
		else if (state == OUT) {
			state = IN;
			++nw;
		}
	}
	printf("\n\nGeneral statistics : ");
	printf("number of lines, words, and characters : %d %d %d\n", nl, nw, nc);

	rewind(instr);

	lineQ = (short int *) allocate(nl * sizeof(short int));
	for (i=0; i<nl; i++) lineQ[i]=FALSE;

	nw = nrow = ncol = nwxc = 0;
	state = OUT;
	salto = NO;

	i=0;

	while ((c = getc(instr)) != EOF) {

		if(c=='%' || c=='#') {
			while ((c = getc(instr)) != EOF)
				if (c=='\n') break;
			++i;
			continue;
		}

		if (c=='\n' && nw > 0)
			if (salto==NO) {
				++nrow; 
				salto=SI;
				if (ncol != nwxc && nrow>1) {
					printf("\nvalores diferentes : ");
					error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
						nrow, ncol, nwxc);
				}
				ncol = nwxc;
				lineQ[i]=TRUE;
				++i;
				nwxc=0;
			} else {
				++i;
			}

		if (c==' ' || c=='\n' || c=='\t')
			state = OUT;
		else 
			if (state == OUT) {
				state = IN;
				++nw; ++nwxc;
				salto=NO;
			}
	}
	printf("\nValid numbers statistics : ");
	printf("nrow, ncol, nvalues : %d %d %d\n", nrow, ncol, nw);

	rewind(instr);

	npoint=nrow;
	row = (realptr) allocate(ncol*sizeof(real));

	*npts = npoint;
	gd.xval = (real *) allocate(npoint * sizeof(real));
	gd.yval = (real *) allocate(npoint * sizeof(real));

	ip = 0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			gd.xval[ip] = row[col1-1];
			gd.yval[ip] = row[col2-1];
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

    fclose(instr);

	fprintf(stdout,"\n... done.\n");
}

void InputData_3c(string filename, int col1, int col2, int col3, 
	int *npts)
{
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
	short int *lineQ;

    instr = stropen(filename, "r");

	fprintf(stdout,
		"\nReading columns %d, %d, and %d from file %s... ",col1,col2,col3,filename);

	state = OUT;
	nl = nw = nc = 0;
	while ((c = getc(instr)) != EOF) {
		++nc;
		if (c=='\n')
			++nl;
		if (c==' ' || c=='\n' || c=='\t')
			state = OUT;
		else if (state == OUT) {
			state = IN;
			++nw;
		}
	}
	printf("\n\nGeneral statistics : ");
	printf("number of lines, words, and characters : %d %d %d\n", nl, nw, nc);

	rewind(instr);

	lineQ = (short int *) allocate(nl * sizeof(short int));
	for (i=0; i<nl; i++) lineQ[i]=FALSE;

	nw = nrow = ncol = nwxc = 0;
	state = OUT;
	salto = NO;

	i=0;

	while ((c = getc(instr)) != EOF) {

		if(c=='%' || c=='#') {
			while ((c = getc(instr)) != EOF)
				if (c=='\n') break;
			++i;
			continue;
		}

		if (c=='\n' && nw > 0)
			if (salto==NO) {
				++nrow; 
				salto=SI;
				if (ncol != nwxc && nrow>1) {
					printf("\nvalores diferentes : ");
					error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
						nrow, ncol, nwxc);
				}
				ncol = nwxc;
				lineQ[i]=TRUE;
				++i;
				nwxc=0;
			} else {
				++i;
			}

		if (c==' ' || c=='\n' || c=='\t')
			state = OUT;
		else 
			if (state == OUT) {
				state = IN;
				++nw; ++nwxc;
				salto=NO;
			}
	}
	printf("\nValid numbers statistics : ");
	printf("nrow, ncol, nvalues : %d %d %d\n", nrow, ncol, nw);

	if (ncol<3)
		error("\n\nInputData_4c: Error : ncol must be >=4\n");

	rewind(instr);

	npoint=nrow;
	row = (realptr) allocate(ncol*sizeof(real));

	*npts = npoint;
	gd.xval = (real *) allocate(npoint * sizeof(real));
	gd.yval = (real *) allocate(npoint * sizeof(real));
	gd.yminval = (real *) allocate(npoint * sizeof(real));

	ip = 0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			gd.xval[ip] = row[col1-1];
			gd.yval[ip] = row[col2-1];
			gd.yminval[ip] = row[col3-1];
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

    fclose(instr);

	fprintf(stdout,"\n... done.\n");
}

void InputData_4c(string filename, int col1, int col2, int col3, int col4, 
	int *npts)
{
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
	short int *lineQ;

    instr = stropen(filename, "r");

	fprintf(stdout,
		"\nReading columns %d, %d, %d, and %d from file %s... ",
		col1,col2,col3,col4,filename);

	state = OUT;
	nl = nw = nc = 0;
	while ((c = getc(instr)) != EOF) {
		++nc;
		if (c=='\n')
			++nl;
		if (c==' ' || c=='\n' || c=='\t')
			state = OUT;
		else if (state == OUT) {
			state = IN;
			++nw;
		}
	}
	printf("\n\nGeneral statistics : ");
	printf("number of lines, words, and characters : %d %d %d\n", nl, nw, nc);

	rewind(instr);

	lineQ = (short int *) allocate(nl * sizeof(short int));
	for (i=0; i<nl; i++) lineQ[i]=FALSE;

	nw = nrow = ncol = nwxc = 0;
	state = OUT;
	salto = NO;

	i=0;

	while ((c = getc(instr)) != EOF) {

		if(c=='%' || c=='#') {
			while ((c = getc(instr)) != EOF)
				if (c=='\n') break;
			++i;
			continue;
		}

		if (c=='\n' && nw > 0)
			if (salto==NO) {
				++nrow; 
				salto=SI;
				if (ncol != nwxc && nrow>1) {
					printf("\nvalores diferentes : ");
					error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
						nrow, ncol, nwxc);
				}
				ncol = nwxc;
				lineQ[i]=TRUE;
				++i;
				nwxc=0;
			} else {
				++i;
			}

		if (c==' ' || c=='\n' || c=='\t')
			state = OUT;
		else 
			if (state == OUT) {
				state = IN;
				++nw; ++nwxc;
				salto=NO;
			}
	}
	printf("\nValid numbers statistics : ");
	printf("nrow, ncol, nvalues : %d %d %d\n", nrow, ncol, nw);

	if (ncol<4)
		error("\n\nInputData_4c: Error : ncol must be >=4\n");

	rewind(instr);

	npoint=nrow;
	row = (realptr) allocate(ncol*sizeof(real));

	*npts = npoint;
	gd.xval = (real *) allocate(npoint * sizeof(real));
	gd.yval = (real *) allocate(npoint * sizeof(real));
	gd.yminval = (real *) allocate(npoint * sizeof(real));
	gd.ymaxval = (real *) allocate(npoint * sizeof(real));

	ip = 0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			gd.xval[ip] = row[col1-1];
			gd.yval[ip] = row[col2-1];
			gd.yminval[ip] = row[col3-1];
			gd.ymaxval[ip] = row[col4-1];
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

    fclose(instr);

	fprintf(stdout,"\n... done.\n");
}

#undef IN
#undef OUT
#undef SI
#undef NO

void TestData(int *npts)
{
    stream outstr;
	int ip;
	real xmin, xmax, dx;
	real yminerr, ymaxerr;

    outstr = stropen("testdata.dat", "w!");
	fprintf(outstr,"# Sample data to test nplot2d\n");
	fprintf(outstr,"# \n");
	*npts = 50;
	fprintf(outstr,"# The number of points generated is %d\n",*npts);
	fprintf(stdout,
		"\nGenerating test data with %d points... ",*npts);

	gd.xval = (real *) allocate(*npts * sizeof(real));
	gd.yval = (real *) allocate(*npts * sizeof(real));
	xmin=1.0;
	xmax=20.0;
	dx=(xmax-xmin)/((real)*npts-1.0);
	for (ip=0; ip<*npts; ip++) {
		gd.xval[ip] = xmin + (real)(ip)*dx;
		gd.yval[ip] = gd.xval[ip] + 5.0*xrandom(0.0,1.0)*rsin(gd.xval[ip]);
		yminerr = xrandom(0.0,1.0);
		ymaxerr = xrandom(0.0,1.0);
		fprintf(outstr,"%g \t%g \t%g \t%g\n", 
			gd.xval[ip], gd.yval[ip], gd.yval[ip]-yminerr, gd.yval[ip]+ymaxerr);
	}
    fclose(outstr);
	fprintf(stdout,"\n... done.\n");
}

void EndRun(void)
{
	fprintf(stdout,"\n\nTotal running cpu time: %gm\n\n",cputime()-gd.cpuinit);
	fclose(gd.outlog);
}
