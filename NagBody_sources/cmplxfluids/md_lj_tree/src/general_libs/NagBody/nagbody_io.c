/*==============================================================================
	MODULE: nagbody_io.c			[NagBody]
	Written by: M.A. Rodriguez-Meza
	Starting date:	January, 2005
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:	stdinc.h, mathfns.h, vectmath
						vectmath.h, getparam.h
						types.h, stat.h, inout.h
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "../general/stdinc.h"
#include "../math/mathfns.h"
#include "../io/inout.h"
#include "../math/vectdefs.h"
#include "../math/vectmath.h"
#include "nagbody.h"
#include "../physics/physconstants.h"

#include <string.h>								// For Unix
//#include "../../../../General_libs/strings.h"	// For Visual C

#include <sys/stat.h>

// BEGIN NEMO REMOVE :: 2019-05-14
// NEMO IO includes ....
//#include "../nemo/include/filestruct.h"
//#include "../nemo/include/history.h"
//
//#include "../nemo/include/snapshot/snapshot.h"
//#include "../nemo/include/snapshot/body.h"
//#include "../nemo/include/snapshot/put_snap.c"
// END:


//----------------------------- COMIENZA BLOQUE I/O -----------------------------

local void infilefmt_string_to_int(string, int *);

local void inputdata_ascii(char *, int, int *, int *, realptr, bool *, char *);
local void inputdata_ascii_long(char *, int, int *, int *, realptr, bool *, char *);
local void inputdata_bin(char *, int, int *, int *, realptr, bool *, char *);
local void inputdata_pv(char *, int, int *, int *, realptr, bool *, char *);
local void inputdata_gadget11_bin(char *, int, int *, int *, realptr, bool *,
	char *, short);
local void inputdata_gadget11_bin_double(char *, int, int *, int *, realptr,
	bool *, char *, short);
local void inputdata_gadget11_ascii(char *, int, int *, int *, realptr,
	bool *, char *, short);
// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
local void inputdata_gadget11_ascii_long(char *, int, int *, int *, realptr,
	bool *, char *, short);
//
local void inputdata_gadget11_bin_double_reducido(char *, int, int *, int *,
	realptr, bool *, char *, short);
local void inputdata_gadget11_ascii_reducido(char *, int, int *, int *, 
	realptr, bool *, char *, short);
local void inputdata_heitmann_ascii(char *, int, int *, int *, realptr,
									bool *, char *);
local void inputdata_heitmann_ascii_long(char *, int, int *, int *, realptr,
									bool *, char *);
local void inputdata_gadget11_bin_swab(char *, int, int *, int *, realptr,
	bool *, char *, short);
local void readin_header(FILE *);

local void readin_header_reducido(FILE *);

local void outfilefmt_string_to_int(string,int *);
local void outputdata_ascii(char *, int, int, real, io_header_blj *, char *);
// Particle data structure to manipulate I/O
// N > 10^6 purpose...
local void outputdata_ascii_long(char *, int, int, real, io_header_blj *, char *);
//
local void outputpvdata(char *, int, int, real, io_header_blj *, char *);
local void outputbindata(char *, int, int, real, io_header_blj *, char *);
local void outputdata_gadget11_bin_double(char *, int, int, real, io_header_blj *, char *);

// GADGET207 COMIENZO //////////////////////////////////////////////////////////
//local void outputdata_gadget207_bin(char *, int, int, real, io_header_blj *, char *);
// GADGET207 FIN ///////////////////////////////////////////////////////////

// IBERO COMIENZO //////////////////////////////////////////////////////////
//local void outputdata_gadget11_bin_ibero(char *, int, int, real, io_header_blj *, char *);
// IBERO FIN ///////////////////////////////////////////////////////////

local void outputdata_gadget11_ascii(char *, int, int, real, io_header_blj *, char *);
local void outputdata_gadget11_normal_body_ascii_long(char *, int, int, real, io_header_blj *, char *);
local void save_ic_gdgi_format_sph(char *, int, int, real, io_header_blj *, char *);
local void save_ic_gdgi_format_normal_body(char *, int, int, real, io_header_blj *, char *);
local void fprint_header(FILE *);
local void fprint_header_reducido(FILE *);

local void outputdata_gadget11_bin_double_reducido(char *, int, int, real, io_header_blj *, char *);
local void outputdata_gadget11_ascii_reducido(char *, int, int, real, io_header_blj *,
	char *, short);

local void outputdata_tipsy_bin(char *, int, int, real, io_header_blj *, char *);
local void outputdata_powmes_ascii(char *, int, int, real, io_header_blj *, char *);
local void outputdata_powmes_ascii_long(char *, int, int, real, io_header_blj *, char *);
// BEGIN NEMO REMOVE :: 2019-05-14
//local void outputdata_nemo(char *, int, int, real, io_header_blj *, char *);
// END

local void inputdata_blj_ascii(char *, int, int *, int *, realptr, bool *,
	io_header_blj *, char *);
local void inputdata_tlj_ascii(char *, int, int *, int *, realptr, bool *,
	io_header_tlj *, char *);
local void inputdata_blj_bin(char *, int, int *, int *, realptr, bool *,
	io_header_blj *, char *);
local void inputdata_blj_pv(char *, int, int *, int *, realptr, bool *,
	io_header_blj *, char *);
local void input_header_blj(stream, io_header_blj_ptr);
local void input_header_tlj(stream, io_header_tlj_ptr);
local void input_header_blj_pv(stream, io_header_blj_ptr);
local void input_header_blj_bin(stream, io_header_blj_ptr);

local void outputdata_blj(char *, int, int, real, io_header_blj *, char *);
local void outputdata_tlj(char *, int, int, real, io_header_tlj *, char *);
local void outputpvdata_blj(char *, int, int, real, io_header_blj *, char *);
local void outputbindata_blj(char *, int, int, real, io_header_blj *, char *);
local void output_header_blj(stream, io_header_blj);
local void output_header_tlj(stream, io_header_tlj);
local void output_header_blj_pv(stream, io_header_blj);
local void output_header_blj_bin(stream, io_header_blj);

local void HeaderConversion_blj(int, int, void *, io_header_blj *);
local void HeaderConversion_tlj(int, int, void *, io_header_tlj *);

local void scantree(bodyptr, int, bodyptr, int *, int, global_data_tree *);
local void scan_walktree(bodyptr, nodeptr, int);

#define filedump	"dump-data"
#define IN 1
#define OUT 0
#define SI 1
#define NO 0

bool readin_pl_file(string filename, int col1, int col2, int row1, int row2, 
	int *npts)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip, jp;
	short int *lineQ;
	bool flag;

	if (row2==0) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file : Error : row2 must be != 0\n");
		return flag;
	}

	if (row1<1) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file : Error : row1 must be > 0\n");
		return flag;
	}

    if (stat(filename, &buf) != 0) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file : Error : file %s doesnot exist!\n",filename);
		return flag;
    } else {
		flag = TRUE;
        instr = stropen(filename, "r");
	}

	printf("\nReading columns %d and %d from rows %d to %d of the file %s... ",
		col1, col2, row1, row2, filename);

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

	if (row2 < 0 )
		npoint=-row2;
	else
		npoint=row2-row1+1;
	row = (realptr) allocate(ncol*sizeof(real));

	*npts = npoint;
	npltd.xval = (real *) allocate(npoint * sizeof(real));
	npltd.yval = (real *) allocate(npoint * sizeof(real));

	ip = 0; jp=0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			if ((ip+1>=row1 && ip+1<=row2) && row2>=1) {
				npltd.xval[jp] = row[col1-1];
				npltd.yval[jp] = row[col2-1];
				++jp;
			} else
				if (ip+1>nrow+row2 && row2<0) {
					npltd.xval[jp] = row[col1-1];
					npltd.yval[jp] = row[col2-1];
					++jp;
				}
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

	fprintf(stdout,"\nNumber of points read %d",jp);
	if (npoint > jp)
		fprintf(stdout,
			"\nreadin_pl_file : Warning : number of points read is less than requested!");
	if (jp==0) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file : Error : number of points read is zero!");
		return flag;
	}
	*npts = jp;

    fclose(instr);
	return flag;

	printf("\n... done.\n");
}

bool readin_pl_file_3c(string filename, int col1, int col2, int col3, 
	int row1, int row2, int *npts)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip, jp;
	short int *lineQ;
	bool flag;

	if (row2==0) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file_3c : Error : row2 must be != 0\n");
		return flag;
	}

	if (row1<1) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file_3c : Error : row1 must be > 0\n");
		return flag;
	}

    if (stat(filename, &buf) != 0) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file_3c : Error : file %s doesnot exist!\n",filename);
		return flag;
    } else {
		flag = TRUE;
        instr = stropen(filename, "r");
	}

	printf("\nReading columns %d, %d, and %d from rows %d to %d of the file %s... ",
		col1, col2, col3, row1, row2, filename);

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

	if (ncol < 3)
		error("\n\nreadin_pl_file_3c : Error : ncol must be >= 3\n");
	rewind(instr);

	if (row2 < 0 )
		npoint=-row2;
	else
		npoint=row2-row1+1;
	row = (realptr) allocate(ncol*sizeof(real));

	*npts = npoint;
	npltd.xval = (real *) allocate(npoint * sizeof(real));
	npltd.yval = (real *) allocate(npoint * sizeof(real));
	npltd.zval = (real *) allocate(npoint * sizeof(real));

	ip = 0; jp=0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			if ((ip+1>=row1 && ip+1<=row2) && row2>=1) {
				npltd.xval[jp] = row[col1-1];
				npltd.yval[jp] = row[col2-1];
				npltd.zval[jp] = row[col3-1];
				++jp;
			} else
				if (ip+1>nrow+row2 && row2<0) {
					npltd.xval[jp] = row[col1-1];
					npltd.yval[jp] = row[col2-1];
					npltd.zval[jp] = row[col3-1];
					++jp;
				}
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

	fprintf(stdout,"\nNumber of points read %d",jp);
	if (npoint > jp)
		fprintf(stdout,
			"\nreadin_pl_file_3c : Warning : number of points read is less than requested!");
	if (jp==0) {
		flag = FALSE;
		fprintf(stdout,
			"\nreadin_pl_file_3c : Error : number of points read is zero!");
		return flag;
	}
	*npts = jp;

    fclose(instr);
	return flag;

	printf("\n... done.\n");
}

void readin_pl_snap(string filename, int col1, int col2, int *npts)
{
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
	short int *lineQ;
	char firstline[150];

    instr = stropen(filename, "r");

	fprintf(stdout,"\nJumping header ... [2 lines] ... ");
	fgets(firstline,200,instr);
	fgets(firstline,200,instr);

	fprintf(stdout,"\nReading columns %d and %d from file %s... ",
		col1,col2,filename);

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

	fprintf(stdout,"\nJumping header again ... [2 lines] ... ");
	fgets(firstline,200,instr);
	fgets(firstline,200,instr);

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

	fprintf(stdout,"\nJumping header once more ... [2 lines] ... ");
	fgets(firstline,200,instr);
	fgets(firstline,200,instr);

	npoint=nrow;
	row = (realptr) allocate(ncol*sizeof(real));

	*npts = npoint;
	npltd.xval = (real *) allocate(npoint * sizeof(real));
	npltd.yval = (real *) allocate(npoint * sizeof(real));

	ip = 0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			npltd.xval[ip] = row[col1-1];
			npltd.yval[ip] = row[col2-1];
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

    fclose(instr);

	printf("\n... done.\n");
}

void readin_pl_snap_4c(string filename, int col1, int col2, 
	int col3, int col4, int *npts)
{
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
	short int *lineQ;

    instr = stropen(filename, "r");

	fprintf(stdout,"\nReading columns %d and %d from file %s... ",
		col1,col2,filename);

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
	npltd.xval = (real *) allocate(npoint * sizeof(real));
	npltd.yval = (real *) allocate(npoint * sizeof(real));
	npltd.zval = (real *) allocate(npoint * sizeof(real));
	npltd.wval = (real *) allocate(npoint * sizeof(real));

	ip = 0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			npltd.xval[ip] = row[col1-1];
			npltd.yval[ip] = row[col2-1];
			npltd.zval[ip] = row[col3-1];
			npltd.wval[ip] = row[col4-1];
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

    fclose(instr);

	printf("\n... done.\n");
}

void readin_pl_snap_type(string filename, int col1, int col2, int col3, int *npts)
{
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
	short int *lineQ;
	char firstline[150];

    instr = stropen(filename, "r");

	fprintf(stdout,"\nJumping header ... [2 lines] ... ");
	fgets(firstline,200,instr);
	fgets(firstline,200,instr);

	printf("\nReading columns %d, %d, and %d from file %s... ",col1,col2,col3,filename);

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

	fprintf(stdout,"\nJumping header again ... [2 lines] ... ");
	fgets(firstline,200,instr);
	fgets(firstline,200,instr);

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

#ifdef THREEDIM
	if (ncol<7)
		error("\n\nreadin_pl_snap_type : 3D : there is not enough columns in file\n");
#else
	if (ncol<6)
		error("\n\nreadin_pl_snap_type : 2D : there is not enough columns in file\n");
#endif

	rewind(instr);

	fprintf(stdout,"\nJumping header once more ... [2 lines] ... ");
	fgets(firstline,200,instr);
	fgets(firstline,200,instr);

	npoint=nrow;
	row = (realptr) allocate(ncol*sizeof(real));

	*npts = npoint;
	npltd.xval = (real *) allocate(npoint * sizeof(real));
	npltd.yval = (real *) allocate(npoint * sizeof(real));
	npltd.Typeval = (int *) allocate(npoint * sizeof(int));

	ip = 0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			npltd.xval[ip]	= row[col1-1];
			npltd.yval[ip]	= row[col2-1];
			npltd.Typeval[ip] = row[col3-1];
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}

    fclose(instr);

	printf("\n... done.\n");
}

void readin_pl_snap_3d(string filename, int col1, int col2, int col3, int *npts)
{
    stream instr;
    int ncol, nrow;
	real *row;
	int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
	short int *lineQ;

    instr = stropen(filename, "r");

	printf("\nReading columns %d, %d, and %d from file %s... ",col1,col2,col3,filename);

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
	npltd.xval = (real *) allocate(npoint * sizeof(real));
	npltd.yval = (real *) allocate(npoint * sizeof(real));
	npltd.zval = (real *) allocate(npoint * sizeof(real));

	ip = 0;
	for (i=0; i<nl; i++) {
		if (lineQ[i]) {
			in_vector_ndim(instr, row, ncol);
			npltd.xval[ip] = row[col1-1];
			npltd.yval[ip] = row[col2-1];
			npltd.zval[ip] = row[col3-1];
			++ip;
		} else {
			while ((c = getc(instr)) != EOF)		// Reading dummy line ...
				if (c=='\n') break;
		}
	}
    fclose(instr);
	printf("\n... done.\n");
}

#undef filedump
#undef IN
#undef OUT
#undef SI
#undef NO



void PrintSnap(string savesnap, string savesnaptmp, int nbody, char *options,
	char mode[2])
{
    stream outstr_snap;
    bodyptr p;
    char   buf[200];

	fprintf(stdout,"\nPrintSnap: mode : %s\n",mode);
    outstr_snap = stropen(savesnaptmp, mode);

#if (NDIM==3)
	fprintf(outstr_snap,"%1s%10s%14s%15s%15s%15s%15s%7s%5s",
			"#","Posx","Posy","Posz","Velx","Vely","Velz","Id","Type");
	if (scanopt(options, "out-phi"))
		fprintf(outstr_snap,"%12s","Phi");
	if (scanopt(options, "out-acc"))
		fprintf(outstr_snap,"%12s%12s%12s","Accx","Accy","Accz");
	fprintf(outstr_snap,"\n");
	fprintf(outstr_snap,
			"%1s%10s%14s%15s%15s%15s%15s%8s%4s",
			"#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>");
	if (scanopt(options, "out-phi")) {
		fprintf(outstr_snap,"%12s","<9>");
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s%12s","<10>","<11>","<12>");
	} else
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s%12s","<9>","<10>","<11>");
	fprintf(outstr_snap,"\n");
#else
#if (NDIM==2)
	fprintf(outstr_snap,"%1s%10s%14s%15s%15s%7s%5s",
			"#","Posx","Posy","Velx","Vely","Id","Type");
	if (scanopt(options, "out-phi"))
		fprintf(outstr_snap,"%12s","Phi");
	if (scanopt(options, "out-acc"))
		fprintf(outstr_snap,"%12s%12s","Accx","Accy");
	fprintf(outstr_snap,"\n");
	fprintf(outstr_snap,
			"%1s%10s%14s%15s%15s%8s%4s",
			"#","<1>","<2>","<3>","<4>","<5>","<6>");
	if (scanopt(options, "out-phi")) {
		fprintf(outstr_snap,"%12s","<9>");
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s","<10>","<11>");
	} else
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s","<9>","<10>");
	fprintf(outstr_snap,"\n");
#endif
#endif

	DO_BODY(p, bodytab, bodytab+nbody) {
        out_vector_mar(outstr_snap, Pos(p));             
        out_vector_mar(outstr_snap, Vel(p));
        out_int_mar(outstr_snap, Id(p));
        out_int_mar(outstr_snap, Type(p));
        if (scanopt(options, "out-phi"))
            out_real_mar(outstr_snap, Phi(p));
        if (scanopt(options, "out-acc"))
            out_vector_mar(outstr_snap, Acc(p));
        fprintf(outstr_snap,"\n");
    }
    fclose(outstr_snap);
	if (scanopt(mode, "a")) {
		sprintf(buf,"cp %s %s",savesnaptmp,savesnap);
		printf("\nsystem: %s",buf);
	} else {
		sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
		printf("\nsystem: %s",buf);
	}
	system(buf);
}

void PrintBodies(string savesnap, string savesnaptmp, int nbody, 
	int nbodiesID, int *bodyID, char *options, char mode[2])
{
    stream outstr_snap;
    bodyptr p;
	int i, id;
    char   buf[200];

	for (i=0; i<nbodiesID; i++)
		if (bodyID[i]>nbody)
			error("\n\nPrintBodies: bodyID out of range\n");

	fprintf(stdout,"\nPrintBodies: mode : %s\n",mode);
    outstr_snap = stropen(savesnaptmp, mode);

#if (NDIM==3)
	fprintf(outstr_snap,"%1s%10s%14s%15s%15s%15s%15s%7s%5s",
			"#","Posx","Posy","Posz","Velx","Vely","Velz","Id","Type");
	if (scanopt(options, "out-phi"))
		fprintf(outstr_snap,"%12s","Phi");
	if (scanopt(options, "out-acc"))
		fprintf(outstr_snap,"%12s%12s%12s","Accx","Accy","Accz");
	fprintf(outstr_snap,"\n");
	fprintf(outstr_snap,
			"%1s%10s%14s%15s%15s%15s%15s%8s%4s",
			"#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>");
	if (scanopt(options, "out-phi")) {
		fprintf(outstr_snap,"%12s","<9>");
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s%12s","<10>","<11>","<12>");
	} else
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s%12s","<9>","<10>","<11>");
	fprintf(outstr_snap,"\n");
#else
#if (NDIM==2)
	fprintf(outstr_snap,"%1s%10s%14s%15s%15s%7s%5s",
			"#","Posx","Posy","Velx","Vely","Id","Type");
	if (scanopt(options, "out-phi"))
		fprintf(outstr_snap,"%12s","Phi");
	if (scanopt(options, "out-acc"))
		fprintf(outstr_snap,"%12s%12s","Accx","Accy");
	fprintf(outstr_snap,"\n");
	fprintf(outstr_snap,
			"%1s%10s%14s%15s%15s%8s%4s",
			"#","<1>","<2>","<3>","<4>","<5>","<6>");
	if (scanopt(options, "out-phi")) {
		fprintf(outstr_snap,"%12s","<9>");
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s","<10>","<11>");
	} else
		if (scanopt(options, "out-acc"))
			fprintf(outstr_snap,"%12s%12s","<9>","<10>");
	fprintf(outstr_snap,"\n");
#endif
#endif

	DO_BODY(p, bodytab, bodytab+nbody) {
		for (i=0; i<nbodiesID; i++) {	// Si no estan ordenados los cuerpos NO SIRVE!!
			if (Id(p)==bodyID[i]) {
				out_vector_mar(outstr_snap, Pos(p));
				out_vector_mar(outstr_snap, Vel(p));
				out_int_mar(outstr_snap, Id(p));
				out_int_mar(outstr_snap, Type(p));
				if (scanopt(options, "out-phi"))
					out_real_mar(outstr_snap, Phi(p));
				if (scanopt(options, "out-acc"))
					out_vector_mar(outstr_snap, Acc(p));
				fprintf(outstr_snap,"\n");
			}
		}
    }
    fclose(outstr_snap);

	if (scanopt(mode, "a")) {
		sprintf(buf,"cp %s %s",savesnaptmp,savesnap);
		printf("\nsystem: %s",buf);
	} else {
		sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
		printf("\nsystem: %s",buf);
	}
	system(buf);
}

void PrintBodiesSets(string savesnap, string savesnaptmp, int nbody, 
	int nbodiesSets, int *bodyIDMin, int *bodyIDMax, real RMax, char *options)
{
    stream outstr;
    bodyptr p;
	int i, id;
    char   buf[200];
	vector cmpos, tmpv;
	real xi, yi;
	real tmass,r;

	for (i=0; i<nbodiesSets; i++)
		if (bodyIDMin[i]>nbody || bodyIDMax[i]>nbody)
			error("\n\nPrintBodiesSets: IDMin or IDMax out of range\n");

    outstr = stropen(savesnaptmp, "w!");

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+nbody) {
		if ( Id(p)>=bodyIDMin[0] && Id(p)<=bodyIDMax[0] ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nBulge CM: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	xi=cmpos[0]; yi=cmpos[1]; 
	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=bodyIDMin[1] && Id(p)<=bodyIDMax[1]) &&
			 r <= RMax ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nDisk CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	xi=cmpos[0]; yi=cmpos[1];
	DO_BODY(p, bodytab, bodytab+nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if (r<=RMax) {
		for (i=0; i<nbodiesSets; i++) {
			if (Id(p)>=bodyIDMin[i] && Id(p)<=bodyIDMax[i]) {
				out_vector_mar(outstr, Pos(p));
				out_vector_mar(outstr, Vel(p));
				out_int_mar(outstr, Id(p));
				out_int_mar(outstr, Type(p));
				if (scanopt(options, "out-phi"))
					out_real_mar(outstr, Phi(p));
				if (scanopt(options, "out-acc"))
					out_vector_mar(outstr, Acc(p));
				fprintf(outstr,"\n");
			}
		}
		}
    }
    fclose(outstr);
	sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
	printf("\nsystem: %s",buf);
	system(buf);
}


//----------------------------- COMIENZA BLOQUE I/O -----------------------------


void inputdata(char *file, char *filefmt, char *outfmt, int step, 
	int *nbodies, int *ndim, realptr tnow, bool *exist_snap, void *hdr, 
	char *options, char *comment)
{
	char buf[100];
	int infmt_int, outfmt_int;
	short allocate_mode;

	sprintf(buf,"Input from file %s in %s format",file,filefmt);
	strcpy(comment,buf);

	infilefmt_string_to_int(filefmt, &infmt_int);
    outfilefmt_string_to_int(outfmt, &outfmt_int);

//printf("\nAqui voy dentro de inputdata(2)\n");
	switch(infmt_int) {
		case IO_SNAP_FMT:
			inputdata_ascii(file, step, nbodies, ndim, tnow, exist_snap, 
				options);
			break;
		case IO_SNAP_FMT_LONG:
			inputdata_ascii_long(file, step, nbodies, ndim, tnow, exist_snap, 
				options);
			break;
		case IO_NULL_FMT: 
			inputdata_ascii(file, step, nbodies, ndim, tnow, exist_snap,
				options);
			break;
		case IO_PV_FMT:
			inputdata_pv(file, step, nbodies, ndim, tnow, exist_snap,
				options);
			break;
		case IO_SNAP_FMT_BIN: 
			inputdata_bin(file, step, nbodies, ndim, tnow, exist_snap,
				options);
			break;
		case IO_GADGET11_FMT_BIN:
			allocate_mode = 1;
			inputdata_gadget11_bin(file, step, nbodies, ndim, tnow, exist_snap,
				options, allocate_mode);
			break;
		case IO_GADGET11_FMT_ASCII:
			allocate_mode = 1;
			inputdata_gadget11_ascii(file, step, nbodies, ndim, tnow, exist_snap,
				options, allocate_mode);
			break;
// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
		case IO_GADGET11_FMT_ASCII_LONG:
			allocate_mode = 5;
			inputdata_gadget11_ascii_long(file, step, nbodies, ndim, tnow, exist_snap,
				options, allocate_mode);
			break;
//
		case IO_GADGET11_FMT_ASCII_REDUCIDO:
			allocate_mode = 1;
			inputdata_gadget11_ascii_reducido(file, step, nbodies, ndim, tnow,
				exist_snap, options, allocate_mode);
			break;
		case IO_GADGET11_FMT_BIN_DOUBLE:
			allocate_mode = 2;
			inputdata_gadget11_bin_double(file, step, nbodies, ndim, tnow,
			 exist_snap, options, allocate_mode); break;
		case IO_HEITMANN_FMT_ASCII:
			printf("\nInput data heitmann format\n");
			inputdata_heitmann_ascii(file, step, nbodies, ndim, tnow,
									 exist_snap, options);
			break;
		case IO_HEITMANN_FMT_ASCII_LONG:
			printf("\nInput data heitmann long format\n");
			inputdata_heitmann_ascii_long(file, step, nbodies, ndim, tnow,
									 exist_snap, options);
			break;
		case IO_SNAP_BLJ_FMT:
			inputdata_blj_ascii(file, step, nbodies, ndim, tnow, exist_snap,
				hdr, options);
			break;
		case IO_SNAP_BLJ_PV_FMT:
			inputdata_blj_pv(file, step, nbodies, ndim, tnow, exist_snap,
				hdr, options);
			break;
		case IO_SNAP_BLJ_BIN_FMT:
			inputdata_blj_bin(file, step, nbodies, ndim, tnow, exist_snap,
				hdr, options);
			break;
		case IO_SNAP_TLJ_FMT:
			inputdata_tlj_ascii(file, step, nbodies, ndim, tnow, exist_snap,
				hdr, options);
			break;
		default:
			printf("\n\tinput: Unknown input format...");
			printf("\n\tinput in default snap (ascii) format...\n"); 
			inputdata_ascii(file, step, nbodies, ndim, tnow, exist_snap,
				options);
			break;
	}
}

local void infilefmt_string_to_int(string infmt_str,int *infmt_int)
{
    *infmt_int=-1;
    if (strcmp(infmt_str,"snap-ascii") == 0)
		*infmt_int = IO_SNAP_FMT;
    if (strcmp(infmt_str,"snap-ascii-long") == 0)
		*infmt_int = IO_SNAP_FMT_LONG;
    if (strnull(infmt_str))								
		*infmt_int = IO_NULL_FMT;
    if (strcmp(infmt_str,"snap-pv") == 0)				
		*infmt_int = IO_PV_FMT;
    if (strcmp(infmt_str,"snap-bin") == 0)				
		*infmt_int = IO_SNAP_FMT_BIN;
    if (strcmp(infmt_str,"gadget11-bin") == 0)			
		*infmt_int = IO_GADGET11_FMT_BIN;
    if (strcmp(infmt_str,"gadget11-ascii") == 0)		
		*infmt_int = IO_GADGET11_FMT_ASCII;
// Particle data structure to manipulate I/O
// N > 10^6 purpose...
    if (strcmp(infmt_str,"gadget11-ascii-long") == 0)		
		*infmt_int = IO_GADGET11_FMT_ASCII_LONG;
//
    if (strcmp(infmt_str,"gadget11-bin-double") == 0)
		*infmt_int = IO_GADGET11_FMT_BIN_DOUBLE;
    if (strcmp(infmt_str,"heitmann-ascii") == 0)
		*infmt_int = IO_HEITMANN_FMT_ASCII;
    if (strcmp(infmt_str,"heitmann-ascii-long") == 0)
		*infmt_int = IO_HEITMANN_FMT_ASCII_LONG;
    if (strcmp(infmt_str,"gadget11-bin-swab") == 0)
		*infmt_int = IO_GADGET11_FMT_BIN_SWAB;
    if (strcmp(infmt_str,"gadget11-ascii-reducido") == 0)
		*infmt_int = IO_GADGET11_FMT_ASCII_REDUCIDO;
    if (strcmp(infmt_str,"gadget11-bin-double-reducido") == 0)
		*infmt_int = IO_GADGET11_FMT_BIN_DOUBLE_REDUCIDO;
    if (strcmp(infmt_str,"snap-blj-ascii") == 0)
		*infmt_int = IO_SNAP_BLJ_FMT;
    if (strcmp(infmt_str,"snap-blj-pv") == 0)
		*infmt_int = IO_SNAP_BLJ_PV_FMT;
    if (strcmp(infmt_str,"snap-blj-bin") == 0)
		*infmt_int = IO_SNAP_BLJ_BIN_FMT;
    if (strcmp(infmt_str,"snap-tlj-ascii") == 0)
		*infmt_int = IO_SNAP_TLJ_FMT;
}


local void inputdata_ascii(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

    in_int(instr, nbody);
    if (*nbody < 1)
        error("inputdata: nbody = %d is absurd\n", *nbody);
    in_int(instr, ndim);
    if (*ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", *ndim, NDIM);
    in_real(instr, tnow);
    bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+*nbody)
        in_real(instr, &Mass(p));               
	DO_BODY(p, bodytab, bodytab+*nbody)
        in_vector(instr, Pos(p));               
	DO_BODY(p, bodytab, bodytab+*nbody)
        in_vector(instr, Vel(p));               

    if (scanopt(options, "in-phi"))
        for (p = bodytab; p < bodytab+*nbody; p++)
            in_real(instr, &Phi(p));           
    if (scanopt(options, "in-acc"))
        for (p = bodytab; p < bodytab+*nbody; p++)
            in_vector(instr, Acc(p));  

    fclose(instr);
    if (scanopt(options, "reset-time"))         
        *tnow = 0.0;
	DO_BODY(p, bodytab, bodytab+*nbody) {
        Type(p) = BODY;
		Id(p) = p-bodytab+1;
	}
}

// Aun no funciona ...
local void inputdata_ascii_long(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr_long p;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

    in_int(instr, nbody);
    if (*nbody < 1)
        error("inputdata: nbody = %d is absurd\n", *nbody);
    in_int(instr, ndim);
    if (*ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", *ndim, NDIM);
    in_real(instr, tnow);
    bodytab_long = (bodyptr_long) allocate(*nbody * sizeof(body_long));

	DO_BODY(p, bodytab_long, bodytab_long+*nbody)
        in_real(instr, &Mass_long(p));               
	DO_BODY(p, bodytab_long, bodytab_long+*nbody)
        in_vector(instr, Pos_long(p));               
	DO_BODY(p, bodytab_long, bodytab_long+*nbody)
        in_vector(instr, Vel_long(p));               

    fclose(instr);
    if (scanopt(options, "reset-time"))         
        *tnow = 0.0;
	DO_BODY(p, bodytab_long, bodytab_long+*nbody) {
        Type_long(p) = BODY;
		Id_long(p) = p-bodytab_long+1;
	}
}

// (2016-06-15)
local void inputdata_pvm(char *file, int step, int *nbody, int *ndim,
	realptr tnow, bool *exist_snap, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;
	char gato[1], firstline[20];

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

	printf("\n\nReading file in snap-pv format ...\n");
	fgets(firstline,200,instr);
	fscanf(instr,"%1s",gato);
    in_int(instr, nbody);                      
    if (*nbody < 1)
        error("inputdata: nbody = %d is absurd\n", *nbody);
    in_int(instr, ndim);
    if (*ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", *ndim, NDIM);
    in_real(instr, tnow);
    bodytab = (bodyptr) allocate(*nbody * sizeof(body));
	
	printf("nbody ndim tnow : %d %d %g\n", *nbody, *ndim, *tnow);

	DO_BODY(p, bodytab, bodytab+*nbody) {
//		in_int(instr, &Id(p));
        Id(p) = bodytab - p +1;
//        in_real(instr, &Mass(p));
        in_vector(instr, Pos(p));               
        in_vector(instr, Vel(p));
        in_real(instr, &Mass(p));
	}

    fclose(instr);

    if (scanopt(options, "reset-time"))
        *tnow = 0.0;
	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;

	printf("\ndone reading.\n");
}

local void inputdata_pv(char *file, int step, int *nbody, int *ndim,
                        realptr tnow, bool *exist_snap, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;
    char gato[1], firstline[20];
    
    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
        *exist_snap = FALSE;
        return;
    } else {
        *exist_snap = TRUE;
        instr = stropen(namebuf, "r");
    }
    
    printf("\n\nReading file in snap-pv format ...\n");
    fgets(firstline,200,instr);
    fscanf(instr,"%1s",gato);
    in_int(instr, nbody);
    if (*nbody < 1)
        error("inputdata: nbody = %d is absurd\n", *nbody);
    in_int(instr, ndim);
    if (*ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", *ndim, NDIM);
    in_real(instr, tnow);
    bodytab = (bodyptr) allocate(*nbody * sizeof(body));
    
    printf("nbody ndim tnow : %d %d %g\n", *nbody, *ndim, *tnow);
    
    DO_BODY(p, bodytab, bodytab+*nbody) {
        in_int(instr, &Id(p));
        in_real(instr, &Mass(p));
        in_vector(instr, Pos(p));
        in_vector(instr, Vel(p));
    }
    
    fclose(instr);
    
    if (scanopt(options, "reset-time"))
        *tnow = 0.0;
    DO_BODY(p, bodytab, bodytab+*nbody)
    Type(p) = BODY;
    
    printf("\ndone reading.\n");
}


local void inputdata_bin(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

    in_int_bin(instr, nbody);                      
    if (*nbody < 1)
        error("inputdata: nbody = %d is absurd\n", *nbody);
    in_int_bin(instr, ndim);
    if (*ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", *ndim, NDIM);
    in_real_bin(instr, tnow);
    bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+*nbody)
        in_real_bin(instr, &Mass(p));               
	DO_BODY(p, bodytab, bodytab+*nbody)
        in_vector_bin(instr, Pos(p));               
	DO_BODY(p, bodytab, bodytab+*nbody)
        in_vector_bin(instr, Vel(p));               
    fclose(instr);                              
    if (scanopt(options, "reset-time"))
        *tnow = 0.0;
	DO_BODY(p, bodytab, bodytab+*nbody) {
        Type(p) = BODY;                         
		Id(p) = p-bodytab+1;
	}
}

local void inputdata_gadget11_bin(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options, short allocate_mode) 
{
//#define SKIP gdgt_fread(&blklen,sizeof(int4byte),1,fd);

    char namebuf[256];
    bodyptr p;

  FILE *fd;
  int   i,k,massflag,count;
  float dummy[3];
  int   pc,type ;
  int4byte  intdummy, blklen;
  double u_init;

    sprintf(namebuf, file, step);

  if((fd=fopen(namebuf,"r")))
    {
      fprintf(stdout,"Reading file '%s'\n",namebuf); fflush(stdout);

      SKIP; 
      if(blklen!=256)
	{
	  error("\n\nincorrect header format (1)\n");
	  printf("incorrect header format (1)\n");
	  endrun(888);
	}
      gdgt_fread(&header1,sizeof(header1),1,fd);
      SKIP;
      if(blklen!=256)
	{
	  printf("incorrect header format (2)\n");
	  endrun(889);
	}

      All.TotN_gas  = N_gas  = header1.npart[0];
      All.TotN_halo = header1.npart[1];
      All.TotN_disk = header1.npart[2];
      All.TotN_bulge= header1.npart[3];
      All.TotN_stars= header1.npart[4];

      if(RestartFlag==2)		// Incluir su lectura en startrun de cada codigo...
	{
	  All.Time = All.TimeBegin = header1.time;
	}

      for(i=0, massflag=0;i<5;i++)
	{
	  All.MassTable[i]= header1.mass[i];
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    massflag=1;
	}

      printf("\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\n",
	     All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars);
       
      NumPart = All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars;

// Given normally in the parameter file	  
	  All.PartAllocFactor = 1.6;
	  All.TreeAllocFactor =	0.8;

      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;
      All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;

      printf("Numpart=%d\n", NumPart);

      AllocateMemory(allocate_mode);

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(float),3,fd);

	  for(k=0;k<3;k++)
	    P[i].Pos[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(float),3,fd);

	  for(k=0;k<3;k++)
	    P[i].Vel[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&intdummy, sizeof(int4byte), 1, fd);
	  P[i].ID= intdummy;
	}
      SKIP;
      
      if(massflag)
	SKIP;
      for(type=0, count=1; type<5; type++)
	{
	  if(All.MassTable[type]==0 && header1.npart[type]>0)
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  gdgt_fread(&dummy[0],sizeof(float),1,fd);
      
		  P[count++].Mass=dummy[0];
		}
	    }
	  else
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  P[count++].Mass= All.MassTable[type];
		}
	    }
	}
      if(massflag)
	SKIP;

      if(N_gas)
	{
	  SKIP;
	  for(i=1;i<=N_gas;i++)
	    {
	      gdgt_fread(&dummy[0],sizeof(float),1,fd);
	      
	      SphP[i].EgySpec= dummy[0];

	    }
	  SKIP;

	  if(RestartFlag==2)	// Incluir su lectura en startrun de cada codigo...
	    {
	      SKIP;
	      for(i=1;i<=N_gas;i++)
		{
		  gdgt_fread(&dummy[0],sizeof(float),1,fd);
		  
		  SphP[i].Hsml= dummy[0];

		}
	      SKIP;
	    }
	}
      
      fclose(fd);
      fprintf(stdout,"done with reading.\n"); fflush(stdout);

      
      for(type=0, pc=1; type<5; type++)
	for(i=0; i<header1.npart[type]; i++)
	  P[pc++].Type = type;

      if(RestartFlag==0)		// Incluir su lectura en startrun de cada codigo...
	{
	  if(All.InitGasTemp>0)
	    {
	      u_init = (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.InitGasTemp;
	      u_init*= All.UnitMass_in_g/All.UnitEnergy_in_cgs;  
	      
	      for(i=1;i<=N_gas;i++) 
		{
		  if(SphP[i].EgySpec==0)
		    SphP[i].EgySpec= u_init;
		}
	    }
	}

  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);

	*exist_snap = TRUE;


  fprintf(stdout,"Transfiriendo informacion a la estructura bodytab ... ");
  
	*ndim = NDIM;
	*nbody = NumPart;

	bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	*tnow = header1.time;

	for(i=1;i<=NumPart;i++) {
		p=bodytab+i-1;
		Id(p) = P[i].ID;
		Mass(p) = P[i].Mass;
		for(k=0;k<3;k++) {
			Pos(p)[k] = P[i].Pos[k];
			Vel(p)[k] = P[i].Vel[k];
		}
	}

	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;                         

	fprintf(stdout, "done\n");

	}
  else
    {
	*exist_snap = FALSE;
    }

}

local void inputdata_gadget11_bin_swab(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options, short allocate_mode) 
{
//#define SKIP gdgt_fread(&blklen,sizeof(int4byte),1,fd);

    char namebuf[256];
    bodyptr p;

  FILE *fd;
  int   i,k,massflag,count;
  float dummy[3];
  int   pc,type ;
  int4byte  intdummy, blklen;
  double u_init;

	char source[7] = "ABCDEF";
	char destination[7];
	char dest[7], ablklen[10];
	int4byte destino;

    sprintf(namebuf, file, step);
  if((fd=fopen(namebuf,"r")))
    {
      fprintf(stdout,"gadget11_bin_swab : Reading file '%s'\n",namebuf); fflush(stdout);

gdgt_fread(&dest[0],4*sizeof(char),1,fd);

swab(source, destination, 6);
printf("\n\nswab: %s -> %s\n",source,destination);

source[0] = '2';
source[1] = '1';
source[2] = '0';
swab(source, destination, 6);
printf("\n\nswab [2] : %s -> %s %d %s\n\n",source,destination, source[2],source[6]);

swab(&blklen, &destino, 4);
printf("\n\nswab [3] : %d -> %d\n",blklen,destino);

printf("atoi : %d %d\n",atoi(source),atoi(destination));

printf("\nablklen : %s\n",ablklen);

printf("\ndest : %s %d\n",dest, dest[0]);

      if(blklen!=256)
	{
	  error("\n\nincorrect header format (1)\n");
	  printf("incorrect header format (1)\n");
	  endrun(888);
	}

error("\n\nI pass the first test!!!!\n");

      gdgt_fread(&header1,sizeof(header1),1,fd);
      SKIP;
      if(blklen!=256)
	{
	  printf("incorrect header format (2)\n");
	  endrun(889);
	}

      All.TotN_gas  = N_gas  = header1.npart[0];
      All.TotN_halo = header1.npart[1];
      All.TotN_disk = header1.npart[2];
      All.TotN_bulge= header1.npart[3];
      All.TotN_stars= header1.npart[4];

      if(RestartFlag==2)		// Incluir su lectura en startrun de cada codigo...
	{
	  All.Time = All.TimeBegin = header1.time;
	}

      for(i=0, massflag=0;i<5;i++)
	{
	  All.MassTable[i]= header1.mass[i];
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    massflag=1;
	}

      printf("\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\n",
	     All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars);
       
      NumPart = All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars;

	  All.PartAllocFactor = 1.6;
	  All.TreeAllocFactor =	0.8;

      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;
      All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;

      printf("Numpart=%d\n", NumPart);
	  
      AllocateMemory(allocate_mode);

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(float),3,fd);

	  for(k=0;k<3;k++)
	    P[i].Pos[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(float),3,fd);

	  for(k=0;k<3;k++)
	    P[i].Vel[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&intdummy, sizeof(int4byte), 1, fd);
	  P[i].ID= intdummy;
	}
      SKIP;
      
      if(massflag)
	SKIP;
      for(type=0, count=1; type<5; type++)
	{
	  if(All.MassTable[type]==0 && header1.npart[type]>0)
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  gdgt_fread(&dummy[0],sizeof(float),1,fd);
      
		  P[count++].Mass=dummy[0];
		}
	    }
	  else
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  P[count++].Mass= All.MassTable[type];
		}
	    }
	}
      if(massflag)
	SKIP;

      if(N_gas)
	{
	  SKIP;
	  for(i=1;i<=N_gas;i++)
	    {
	      gdgt_fread(&dummy[0],sizeof(float),1,fd);
	      
	      SphP[i].EgySpec= dummy[0];

	    }
	  SKIP;

	  if(RestartFlag==2)		// Incluir su lectura en startrun de cada codigo...
	    {
	      SKIP;
	      for(i=1;i<=N_gas;i++)
		{
		  gdgt_fread(&dummy[0],sizeof(float),1,fd);
		  SphP[i].Hsml= dummy[0];
		}
	      SKIP;
	    }
	}
      
      fclose(fd);
      fprintf(stdout,"done with reading.\n"); fflush(stdout);

      for(type=0, pc=1; type<5; type++)
	for(i=0; i<header1.npart[type]; i++)
	  P[pc++].Type = type;

      if(RestartFlag==0)		// Incluir su lectura en startrun de cada codigo...
	{
	  if(All.InitGasTemp>0)
	    {
	      u_init = (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.InitGasTemp;
	      u_init*= All.UnitMass_in_g/All.UnitEnergy_in_cgs;
	      
	      for(i=1;i<=N_gas;i++) 
		{
		  if(SphP[i].EgySpec==0)
		    SphP[i].EgySpec= u_init;
		}
	    }
	}

  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);

	*exist_snap = TRUE;

  fprintf(stdout,"Transfiriendo informacion a la estructura bodytab ... ");
  
	*ndim = NDIM;
	*nbody = NumPart;

	bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	*tnow = header1.time;

	for(i=1;i<=NumPart;i++) {
		p=bodytab+i-1;
		Id(p) = P[i].ID;
		Mass(p) = P[i].Mass;
		for(k=0;k<3;k++) {
			Pos(p)[k] = P[i].Pos[k];
			Vel(p)[k] = P[i].Vel[k];
		}
	}

	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;                         

	fprintf(stdout, "done\n");

	}
  else
    {
	*exist_snap = FALSE;
    }
}

local void inputdata_gadget11_bin_double(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options, short allocate_mode) 
{
//#define SKIP gdgt_fread(&blklen,sizeof(int4byte),1,fd);

    char namebuf[256];
    bodyptr p;

  FILE *fd;
  int   i,k,massflag,count;
  double dummy[3];
  int   pc,type ;
  int4byte  intdummy, blklen;
  double u_init;


	NumParticleTypes=6;

    sprintf(namebuf, file, step);

  if((fd=fopen(namebuf,"r")))
    {
      fprintf(stdout,"Reading file (double) '%s'\n",namebuf); fflush(stdout);

      SKIP; 
      if(blklen!=256)
	{
	  error("\n\nincorrect header format (1)\n");
	  printf("incorrect header format (1)\n");
	  endrun(888);
	}
      gdgt_fread(&header1,sizeof(header1),1,fd);
      SKIP;
      if(blklen!=256)
	{
	  printf("incorrect header format (2)\n");
	  endrun(889);
	}

      All.TotN_gas  = N_gas  = header1.npart[0];
      All.TotN_halo = header1.npart[1];
      All.TotN_disk = header1.npart[2];
      All.TotN_bulge= header1.npart[3];
      All.TotN_stars= header1.npart[4];
      All.TotN_dm	= header1.npart[5];

      if(RestartFlag==2)		// Incluir su lectura en startrun de cada codigo...
	{
	  All.Time = All.TimeBegin = header1.time;
	}

      for(i=0, massflag=0;i<NumParticleTypes;i++)
	{
	  All.MassTable[i]= header1.mass[i];
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    massflag=1;
	}

      printf("\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\nN_dm: %d\n",
	     All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars,
		 All.TotN_dm);
       
      NumPart = All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars + All.TotN_dm;

// Given normally in the parameter file	
// Estos parametros debe ser dados en la linea de comandos...
	  All.PartAllocFactor = 1.6;
	  All.TreeAllocFactor =	0.8;

      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;
      All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;

      printf("Numpart=%d\n", NumPart);

      AllocateMemory(allocate_mode);

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(double),3,fd);

	  for(k=0;k<3;k++)
	    P_double[i].Pos[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(double),3,fd);

	  for(k=0;k<3;k++)
	    P_double[i].Vel[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&intdummy, sizeof(int4byte), 1, fd);
	  P_double[i].ID= intdummy;
	}
      SKIP;
      
      if(massflag)
	SKIP;
      for(type=0, count=1; type<NumParticleTypes; type++)
	{
	  if(All.MassTable[type]==0 && header1.npart[type]>0)
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  gdgt_fread(&dummy[0],sizeof(double),1,fd);
      
		  P_double[count++].Mass=dummy[0];
		}
	    }
	  else
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  P_double[count++].Mass= All.MassTable[type];
		}
	    }
	}
      if(massflag)
	SKIP;

      if(N_gas)
	{
	  SKIP;
	  for(i=1;i<=N_gas;i++)
	    {
	      gdgt_fread(&dummy[0],sizeof(double),1,fd);
	      SphP_double[i].EgySpec= dummy[0];
	    }
	  SKIP;

	  if(RestartFlag==2)		// Incluir su lectura en startrun de cada codigo...
	    {
	      SKIP;
	      for(i=1;i<=N_gas;i++)
		{
		  gdgt_fread(&dummy[0],sizeof(double),1,fd);
		  
		  SphP_double[i].Hsml= dummy[0];

		}
	      SKIP;
	    }
	}
      
      fclose(fd);
      fprintf(stdout,"done with reading.\n"); fflush(stdout);

	for(type=0, pc=1; type<NumParticleTypes; type++)
	for(i=0; i<header1.npart[type]; i++)
	  P_double[pc++].Type = type;

      if(RestartFlag==0)		// Incluir su lectura en startrun de cada codigo...
	{
	  if(All.InitGasTemp>0)
	    {
	      u_init = (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.InitGasTemp;
	      u_init*= All.UnitMass_in_g/All.UnitEnergy_in_cgs;
	      
	      for(i=1;i<=N_gas;i++) 
		{
		  if(SphP_double[i].EgySpec==0)
		    SphP_double[i].EgySpec= u_init;
		}
	    }
	}

  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);

	*exist_snap = TRUE;

  fprintf(stdout,"Transfiriendo informacion a la estructura bodytab ... ");
  
	*ndim = NDIM;
	*nbody = NumPart;

	bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	*tnow = header1.time;

	for(i=1;i<=NumPart;i++) {
		p=bodytab+i-1;
		Id(p) = P_double[i].ID;
		Mass(p) = P_double[i].Mass;
		for(k=0;k<3;k++) {
			Pos(p)[k] = P_double[i].Pos[k];
			Vel(p)[k] = P_double[i].Vel[k];
		}
	}

	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;                         

	fprintf(stdout, "done\n");

	}
  else
    {
	*exist_snap = FALSE;
    }

}

local void inputdata_gadget11_bin_double_reducido(char *file, int step, 
	int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options, short allocate_mode) 
{
//#define SKIP gdgt_fread(&blklen,sizeof(int4byte),1,fd);

    char namebuf[256];
    bodyptr p;

  FILE *fd;
  int   i,k,massflag,count;
  double dummy[3];
  int   pc,type ;
  int4byte  intdummy, blklen;
  double u_init;

	NumParticleTypes=6;

    sprintf(namebuf, file, step);

  if((fd=fopen(namebuf,"r")))
    {
      fprintf(stdout,"Reading file (double) '%s'\n",namebuf); fflush(stdout);

      SKIP; 
      if(blklen!=256)
	{
	  error("\n\nincorrect header format (1)\n");
	  printf("incorrect header format (1)\n");
	  endrun(888);
	}
      gdgt_fread(&header_reducido,sizeof(header_reducido),1,fd);
      SKIP;
      if(blklen!=256)
	{
	  printf("incorrect header format (2)\n");
	  endrun(889);
	}

N_gas=0;
      All.TotN_p0  = header_reducido.npart[0];
      All.TotN_halo = header_reducido.npart[1];
      All.TotN_disk = header_reducido.npart[2];
      All.TotN_bulge= header_reducido.npart[3];
      All.TotN_stars= header_reducido.npart[4];
      All.TotN_dm	= header_reducido.npart[5];

      if(RestartFlag==2)		// Incluir su lectura en startrun de cada codigo...
	{
	  All.Time = All.TimeBegin = header_reducido.time;
	}

      for(i=0, massflag=0;i<NumParticleTypes;i++)
	{
	  All.MassTable[i]= header_reducido.mass[i];
	  if(All.MassTable[i]==0 && header_reducido.npart[i]>0)
	    massflag=1;
	}

      printf("\nN_p0: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\nN_dm: %d\n",
	     All.TotN_p0, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars,
		 All.TotN_dm);

      NumPart = All.TotNumPart =    All.TotN_p0  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars + All.TotN_dm;

// Given normally in the parameter file	
// Estos parametros debe ser dados en la linea de comandos...
	  All.PartAllocFactor = 1.6;
	  All.TreeAllocFactor =	0.8;

      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;

      printf("Numpart=%d\n", NumPart);

      AllocateMemory(allocate_mode);

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(double),3,fd);

	  for(k=0;k<3;k++)
	    P_double[i].Pos[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&dummy[0],sizeof(double),3,fd);

	  for(k=0;k<3;k++)
	    P_double[i].Vel[k]=dummy[k];
	}
      SKIP;

      SKIP;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fread(&intdummy, sizeof(int4byte), 1, fd);
	  P_double[i].ID= intdummy;
	}
      SKIP;
      
      if(massflag)
	SKIP;
      for(type=0, count=1; type<NumParticleTypes; type++)
	{
	  if(All.MassTable[type]==0 && header_reducido.npart[type]>0)
	    {
	      for(i=1;i<=header_reducido.npart[type];i++)
		{
		  gdgt_fread(&dummy[0],sizeof(double),1,fd);
      
		  P_double[count++].Mass=dummy[0];
		}
	    }
	  else
	    {
	      for(i=1;i<=header_reducido.npart[type];i++)
		{
		  P_double[count++].Mass= All.MassTable[type];
		}
	    }
	}
      if(massflag)
	SKIP;
      
      fclose(fd);
      fprintf(stdout,"done with reading.\n"); fflush(stdout);

      
      for(type=0, pc=1; type<NumParticleTypes; type++)
	for(i=0; i<header_reducido.npart[type]; i++)
	  P_double[pc++].Type = type;


  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);

	*exist_snap = TRUE;

  fprintf(stdout,"Transfiriendo informacion a la estructura bodytab ... ");
  
	*ndim = NDIM;
	*nbody = NumPart;

	bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	*tnow = header_reducido.time;

	for(i=1;i<=NumPart;i++) {
		p=bodytab+i-1;
		Id(p) = P_double[i].ID;
		Mass(p) = P_double[i].Mass;
		for(k=0;k<3;k++) {
			Pos(p)[k] = P_double[i].Pos[k];
			Vel(p)[k] = P_double[i].Vel[k];
		}
	}

	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;                         

	fprintf(stdout, "done\n");

	}
  else
    {
	*exist_snap = FALSE;
    }
}


local void inputdata_gadget11_ascii(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options, short allocate_mode) 
{
    char namebuf[256];
    bodyptr p;

  FILE *fd;
  int   i,k,massflag,count;
  double dummy[3];
  int   pc,type ;
  int4byte  intdummy; 
  int   blklen;
  double u_init;
  vector vectordummy;
  real realdummy;
  
  real tmass;

	NumParticleTypes=6;

    sprintf(namebuf, file, step);

  if((fd=fopen(namebuf,"r")))
    {
      fprintf(stdout,"Reading file '%s'  %d... ",namebuf,sizeof(header1));

	readin_header(fd);

      All.TotN_gas  = N_gas  = header1.npart[0];
      All.TotN_halo = header1.npart[1];
      All.TotN_disk = header1.npart[2];
      All.TotN_bulge= header1.npart[3];
      All.TotN_stars= header1.npart[4];
      All.TotN_dm= header1.npart[5];

      for(i=0, massflag=0;i<NumParticleTypes;i++)
	{
	  All.MassTable[i]= header1.mass[i];
printf("\nMassTable found : %g",All.MassTable[i]);
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    massflag=1;
	}
printf("\n\nmassflag : %d\n",massflag);

      fprintf(stdout,"\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\nN_dm: %d\n",
	     All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars, All.TotN_dm);

      NumPart = All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars + All.TotN_dm;

// Given normally in the parameter file.  
// MODIFICAR PARA SER INCLUIDO EN ARCHIVO DE PARAMETROS O LINEA DE COMANDOS
	  All.PartAllocFactor = 1.6;
	  All.TreeAllocFactor =	0.8;

      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;    
									  
      All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;

      AllocateMemory(allocate_mode);

      for(i=1;i<=NumPart;i++)
	{
		in_vector(fd,vectordummy);

	  for(k=0;k<3;k++) {
	    P[i].Pos[k]=vectordummy[k];
	  }
	}

      for(i=1;i<=NumPart;i++)
	{
		in_vector(fd,vectordummy);

	  for(k=0;k<3;k++)
	    P[i].Vel[k]=vectordummy[k];
	}

	  for(i=1;i<=NumPart;i++)
	{
		in_int(fd,&intdummy);
	  P[i].ID= intdummy;
	}
      
      for(type=0, count=1; type<NumParticleTypes; type++)
	{
	  if(All.MassTable[type]==0 && header1.npart[type]>0)
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		in_real(fd,&realdummy);
      
		  P[count++].Mass=realdummy;
		}
	    }
	  else
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  P[count++].Mass= All.MassTable[type];
		}
	    }
	}

      if(N_gas)
	{
	  for(i=1;i<=N_gas;i++)
	    {
		in_real(fd,&realdummy);
	      
	      SphP[i].EgySpec= realdummy;

	    }

	}
      
      fclose(fd);
      fprintf(stdout,"done.\n"); fflush(stdout);

      for(type=0, pc=1; type<NumParticleTypes; type++)
	for(i=0; i<header1.npart[type]; i++)
	  P[pc++].Type = type;


      if(All.InitGasTemp>0)
	{
#ifdef ADIABATICIDEALEQUATIONOFSTATE			// Procesos adiabaticos
	  u_init = (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.InitGasTemp;
	  u_init*= All.UnitMass_in_g/All.UnitEnergy_in_cgs;
#else
	  u_init = ((double) 1.5) * SOUNDSPEEDSQR;
#endif
	  
	  for(i=1;i<=N_gas;i++) 
	    {
	      if(SphP[i].EgySpec==0)
		SphP[i].EgySpec= u_init;
	    }
	}

  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);
  fflush(stdout);

	*exist_snap = TRUE;

	fprintf(stdout,"Transfiriendo informacion a la estructura bodytab ... ");

	*ndim = NDIM;
	*nbody = NumPart;

	bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	*tnow = header1.time;

tmass=0.;

	for(i=1;i<=NumPart;i++) {
		p=bodytab+i-1;
		Id(p) = P[i].ID;
		Mass(p) = P[i].Mass;

tmass += Mass(p);

		for(k=0;k<3;k++) {
			Pos(p)[k] = P[i].Pos[k];
			Vel(p)[k] = P[i].Vel[k];
		}
	}

printf("\n\nTotal Mass found : %g\n",tmass);

	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;                         

    }
  else
    {
	*exist_snap = FALSE;
    }

}


// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
local void inputdata_gadget11_ascii_long(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options, short allocate_mode) 
{
    char namebuf[256];
    bodyptr_long p;

  FILE *fd;
  int   i,k,massflag,count;
  double dummy[3];
  int   pc,type ;
  int4byte  intdummy; 
  int   blklen;
  double u_init;
  vector vectordummy;
  real realdummy;
  
  real tmass;

	NumParticleTypes=6;

    sprintf(namebuf, file, step);

  if((fd=fopen(namebuf,"r")))
    {
      fprintf(stdout,"Reading file '%s'  %d... ",namebuf,sizeof(header1));

	readin_header(fd);

      All.TotN_gas  = N_gas  = header1.npart[0];
      All.TotN_halo = header1.npart[1];
      All.TotN_disk = header1.npart[2];
      All.TotN_bulge= header1.npart[3];
      All.TotN_stars= header1.npart[4];
      All.TotN_dm= header1.npart[5];

      for(i=0, massflag=0;i<NumParticleTypes;i++)
	{
	  All.MassTable[i]= header1.mass[i];
printf("\nMassTable found : %g",All.MassTable[i]);
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    massflag=1;
	}
printf("\n\nmassflag : %d\n",massflag);

      fprintf(stdout,"\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\nN_dm: %d\n",
	     All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars, All.TotN_dm);

      NumPart = All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars + All.TotN_dm;

	  All.PartAllocFactor = 1.6;
	  All.TreeAllocFactor =	0.8;

      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;    
									  
      All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;

      AllocateMemory(allocate_mode);

      for(i=1;i<=NumPart;i++)
	{
		in_vector(fd,vectordummy);

	  for(k=0;k<3;k++) {
	    P_long[i].Pos[k]=vectordummy[k];
	  }
	}

      for(i=1;i<=NumPart;i++)
	{
		in_vector(fd,vectordummy);

	  for(k=0;k<3;k++)
	    P_long[i].Vel[k]=vectordummy[k];
	}

	  for(i=1;i<=NumPart;i++)
	{
		in_int(fd,&intdummy);
	  P_long[i].ID= intdummy;
	}
      
      for(type=0, count=1; type<NumParticleTypes; type++)
	{
	  if(All.MassTable[type]==0 && header1.npart[type]>0)
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		in_real(fd,&realdummy);

		  P_long[count++].Mass=realdummy;
		}
	    }
	  else
	    {
	      for(i=1;i<=header1.npart[type];i++)
		{
		  P_long[count++].Mass= All.MassTable[type];
		}
	    }
	}

      if(N_gas)
	{
	  for(i=1;i<=N_gas;i++)
	    {
		in_real(fd,&realdummy);
	      
	      SphP[i].EgySpec= realdummy;

	    }

	}
      
      fclose(fd);
      fprintf(stdout,"done.\n"); fflush(stdout);

      for(type=0, pc=1; type<NumParticleTypes; type++)
	for(i=0; i<header1.npart[type]; i++)
	  P_long[pc++].Type = type;

      if(All.InitGasTemp>0)
	{
#ifdef ADIABATICIDEALEQUATIONOFSTATE			// Procesos adiabaticos
	  u_init = (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.InitGasTemp;
	  u_init*= All.UnitMass_in_g/All.UnitEnergy_in_cgs;
#else
	  u_init = ((double) 1.5) * SOUNDSPEEDSQR;
#endif
	  
	  for(i=1;i<=N_gas;i++) 
	    {
	      if(SphP[i].EgySpec==0)
		SphP[i].EgySpec= u_init;
	    }
	}

  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);
  fflush(stdout);

	*exist_snap = TRUE;

	fprintf(stdout,"Transfiriendo informacion a la estructura bodytab_long ... ");

	*ndim = NDIM;
	*nbody = NumPart;

	bodytab_long = (bodyptr_long) allocate(*nbody * sizeof(body_long));

	*tnow = header1.time;

tmass=0.;

	for(i=1;i<=NumPart;i++) {
		p=bodytab_long+i-1;
		Id_long(p) = P_long[i].ID;
		Mass_long(p) = P_long[i].Mass;

tmass += Mass_long(p);

		for(k=0;k<3;k++) {
			Pos_long(p)[k] = P_long[i].Pos[k];
			Vel_long(p)[k] = P_long[i].Vel[k];
		}
	}

printf("\n\nTotal Mass and bodies found : %g %ld\n",tmass,*nbody);

	DO_BODY(p, bodytab_long, bodytab_long+*nbody)
        Type_long(p) = BODY;                         

    }
  else
    {
	*exist_snap = FALSE;
    }

printf("\nDone with gadget11_ascii_long ...\n");

}
// End inputdata_gadget11_ascii_long


local void inputdata_gadget11_ascii_reducido(char *file, int step, 
	int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options, short allocate_mode) 
{
    char namebuf[256];
    bodyptr p;

  FILE *fd;
  int   i,k,massflag,count;
  double dummy[3];
  int   pc,type ;
  int4byte  intdummy; 
  int   blklen;
  double u_init;
  vector vectordummy;
  real realdummy;

	NumParticleTypes=6;

    sprintf(namebuf, file, step);

  if((fd=fopen(namebuf,"r")))
    {
      fprintf(stdout,"Reading file '%s'  %d... ",namebuf,sizeof(header_reducido));

	readin_header_reducido(fd);

N_gas=0;
      All.TotN_gas  = N_gas;
      All.TotN_p0  = header_reducido.npart[0];
      All.TotN_halo = header_reducido.npart[1];
      All.TotN_disk = header_reducido.npart[2];
      All.TotN_bulge= header_reducido.npart[3];
      All.TotN_stars= header_reducido.npart[4];
      All.TotN_dm= header_reducido.npart[5];

      for(i=0, massflag=0;i<NumParticleTypes;i++)
	{
	  All.MassTable[i]= header_reducido.mass[i];
	  if(All.MassTable[i]==0 && header_reducido.npart[i]>0)
	    massflag=1;
	}

      fprintf(stdout,"\nN_sph: %d\nN_p0: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\nN_dm: %d\n",
	     All.TotN_gas, All.TotN_p0, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars, All.TotN_dm);

      NumPart = All.TotNumPart =   All.TotN_gas + All.TotN_p0  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars + All.TotN_dm;

	  All.PartAllocFactor = 1.6;
	  All.TreeAllocFactor =	0.8;

      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;    
									  
      All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;

      AllocateMemory(allocate_mode);

      for(i=1;i<=NumPart;i++)
	{
		in_vector(fd,vectordummy);

	  for(k=0;k<3;k++)
	    P[i].Pos[k]=vectordummy[k];
	}

      for(i=1;i<=NumPart;i++)
	{
		in_vector(fd,vectordummy);

	  for(k=0;k<3;k++)
	    P[i].Vel[k]=vectordummy[k];
	}

	  for(i=1;i<=NumPart;i++)
	{
		in_int(fd,&intdummy);
	  P[i].ID= intdummy;
	}
      
      if(massflag)
      for(type=0, count=1; type<NumParticleTypes; type++)
	{
	  if(All.MassTable[type]==0 && header_reducido.npart[type]>0)
	    {
	      for(i=1;i<=header_reducido.npart[type];i++)
		{
		in_real(fd,&realdummy);
      
		  P[count++].Mass=realdummy;
		}
	    }
	  else
	    {
	      for(i=1;i<=header_reducido.npart[type];i++)
		{
		  P[count++].Mass= All.MassTable[type];
		}
	    }
	}

    if (scanopt(options, "out-phi"))
		for(i=1;i<=NumPart;i++) {
			in_real(fd,&realdummy);
			P[i].Potential=realdummy;
		}
    if (scanopt(options, "out-acc"))
		for(i=1;i<=NumPart;i++) {
			in_vector(fd,vectordummy);
			for(k=0;k<3;k++)
				P[i].Accel[k]=vectordummy[k];
		}


      fclose(fd);
      fprintf(stdout,"done.\n"); fflush(stdout);

      for(type=0, pc=1; type<NumParticleTypes; type++)
	for(i=0; i<header_reducido.npart[type]; i++)
	  P[pc++].Type = type;

  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);
  fflush(stdout);

	*exist_snap = TRUE;

	fprintf(stdout,"Transfiriendo informacion a la estructura bodytab ... ");

	*ndim = NDIM;
	*nbody = NumPart;

	bodytab = (bodyptr) allocate(*nbody * sizeof(body));

	*tnow = header_reducido.time;

	for(i=1;i<=NumPart;i++) {
		p=bodytab+i-1;
		Id(p) = P[i].ID;
		Mass(p) = P[i].Mass;
		for(k=0;k<3;k++) {
			Pos(p)[k] = P[i].Pos[k];
			Vel(p)[k] = P[i].Vel[k];
			if (scanopt(options, "out-acc"))
				Acc(p)[k] = P[i].Accel[k];
		}
		if (scanopt(options, "out-phi"))
			Phi(p) = P[i].Potential;
	}

	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;                         

    }
  else
    {
	*exist_snap = FALSE;
    }
}


local void readin_header(FILE *fd)
{
	in_int(fd, &header1.npart[0]);
	in_int(fd, &header1.npart[1]);	// Si las particulas estan
									// aqui hay errores en 'energy.txt'
    in_int(fd, &header1.npart[2]);
    in_int(fd, &header1.npart[3]);
    in_int(fd, &header1.npart[4]);	// Aqui si funciona!!!!
    in_int(fd, &header1.npart[5]);

    in_real(fd,&header1.mass[0]);
    in_real(fd,&header1.mass[1]);
    in_real(fd,&header1.mass[2]);
    in_real(fd,&header1.mass[3]);
    in_real(fd,&header1.mass[4]);
    in_real(fd,&header1.mass[5]);

    in_real(fd,&header1.time);
    in_real(fd,&header1.redshift);
    in_int(fd, &header1.flag_sfr);				// Necesaria esta variable?
    in_int(fd, &header1.flag_feedback);

    in_int(fd, &header1.npartTotal[0]);
    in_int(fd, &header1.npartTotal[1]);
    in_int(fd, &header1.npartTotal[2]);
    in_int(fd, &header1.npartTotal[3]);
    in_int(fd, &header1.npartTotal[4]);
    in_int(fd, &header1.npartTotal[5]);

    in_int(fd,&header1.flag_cooling);	// Se queda para tener compatibilidad en el
										// formato de IC-files
    in_int(fd, &header1.num_files);
    in_real(fd,&header1.BoxSize);
    in_real(fd,&header1.Omega0);
    in_real(fd,&header1.OmegaLambda);
    in_real(fd,&header1.HubbleParam);
}


local void readin_header_reducido(FILE *fd)
{
	in_int(fd, &header_reducido.npart[0]);
	in_int(fd, &header_reducido.npart[1]);	// Si las particulas estan
									// aqui hay errores en 'energy.txt'
    in_int(fd, &header_reducido.npart[2]);
    in_int(fd, &header_reducido.npart[3]);
    in_int(fd, &header_reducido.npart[4]);	// Aqui si funciona!!!!
    in_int(fd, &header_reducido.npart[5]);

    in_real(fd,&header_reducido.mass[0]);
    in_real(fd,&header_reducido.mass[1]);
    in_real(fd,&header_reducido.mass[2]);
    in_real(fd,&header_reducido.mass[3]);
    in_real(fd,&header_reducido.mass[4]);
    in_real(fd,&header_reducido.mass[5]);

    in_real(fd,&header_reducido.time);

    in_int(fd, &header_reducido.npartTotal[0]);
    in_int(fd, &header_reducido.npartTotal[1]);
    in_int(fd, &header_reducido.npartTotal[2]);
    in_int(fd, &header_reducido.npartTotal[3]);
    in_int(fd, &header_reducido.npartTotal[4]);
    in_int(fd, &header_reducido.npartTotal[5]);

    in_int(fd, &header_reducido.num_files);
    in_real(fd,&header_reducido.BoxSize);
}


local void inputdata_heitmann_ascii(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;
	char gato[1], firstline[20];
	real tmp;
	long int tmplongmem;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

	fgets(firstline,200,instr);
	fscanf(instr,"%1s",gato);
    in_int(instr, nbody);                      
    if (*nbody < 1)
        error("inputdata_heitmann_ascii: nbody = %d is absurd\n", *nbody);
    in_int(instr, ndim);                       
    if (*ndim != NDIM)
        error("inputdata_heitmann_ascii: ndim = %d; expected %d\n", *ndim, NDIM);
    in_real(instr, tnow);


							// Tener cuidado con pasarse del limite
							// del entero long mas grande (+- 2,147,483,648)
							// que se maneja en una maquina con 32 bits...
	tmplongmem = *nbody * sizeof(body);
fprintf(stdout,"\nSize of body struct %ld; and total bytes reserved %ld\n",
sizeof(body),tmplongmem);
	printf("\n\nAllocating %ld bytes of memory for %d particles\n",
		tmplongmem, *nbody);
    bodytab = (bodyptr) allocate(tmplongmem);

fprintf(stdout,"\nbodytab, bodytab+nbody: %ld %ld\n", bodytab, bodytab + *nbody);

	DO_BODY(p, bodytab, bodytab+*nbody) {

		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii : input conversion error\n");
		else
			Pos(p)[0] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii : input conversion error\n");
		else
			Vel(p)[0] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii : input conversion error\n");
		else
			Pos(p)[1] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii : input conversion error\n");
		else
			Vel(p)[1] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii : input conversion error\n");
		else
			Pos(p)[2] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii : input conversion error\n");
		else
			Vel(p)[2] = tmp;

        in_real(instr, &Mass(p));               
		in_int(instr, &Id(p));

	}
fprintf(stdout,"\nRead %d total particles, last Id %d\n",p-bodytab, Id(p-1));

    fclose(instr);

	DO_BODY(p, bodytab, bodytab+*nbody)
        Type(p) = BODY;
}

local void inputdata_heitmann_ascii_long(char *file, int step, int *nbody, int *ndim, 
									realptr tnow, bool *exist_snap, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
//    bodyptr p;
    bodyptr_long p;
	char gato[1], firstline[20];
	real tmp;
	long int tmplongmem;
	
    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}
	
	fgets(firstline,200,instr);
	fscanf(instr,"%1s",gato);
    in_int(instr, nbody);
    if (*nbody < 1)
        error("inputdata_heitmann_ascii_long: nbody = %d is absurd\n", *nbody);
    in_int(instr, ndim);
    if (*ndim != NDIM)
        error("inputdata_heitmann_ascii_long: ndim = %d; expected %d\n", *ndim, NDIM);
    in_real(instr, tnow);
	
	
	// Tener cuidado con pasarse del limite
	// del entero long mas grande (+- 2,147,483,648)
	// que se maneja en una maquina con 32 bits...
	tmplongmem = *nbody * sizeof(body_long);
	fprintf(stdout,"\nSize of body struct %ld; and total bytes reserved %ld\n",
			sizeof(body_long),tmplongmem);
	printf("\n\nAllocating %ld bytes of memory for %d particles\n",
		   tmplongmem, *nbody);
    bodytab_long = (bodyptr_long) allocate(tmplongmem);
	
	fprintf(stdout,"\nbodytab, bodytab+nbody: %ld %ld\n", bodytab_long, bodytab_long + *nbody);
	
	DO_BODY(p, bodytab_long, bodytab_long+*nbody) {
		
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii_long: input conversion error\n");
		else
			Pos_long(p)[0] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii_long: input conversion error\n");
		else
			Vel_long(p)[0] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii_long: input conversion error\n");
		else
			Pos_long(p)[1] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii_long: input conversion error\n");
		else
			Vel_long(p)[1] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii_long: input conversion error\n");
		else
			Pos_long(p)[2] = tmp;
		if (fscanf(instr, "%lf", &tmp) != 1)
			error("inputdata_heitmann_ascii_long: input conversion error\n");
		else
			Vel_long(p)[2] = tmp;
		
        in_real(instr, &Mass_long(p));               
		in_int(instr, &Id_long(p));
		
	}
	fprintf(stdout,"\nRead %d total particles, last Id %d\n",p-bodytab_long, Id_long(p-1));

    fclose(instr);

	DO_BODY(p, bodytab_long, bodytab_long+*nbody)
	Type_long(p) = BODY;
}

// -------------COMIENZAN I/O ROUTINES FOR BLJ DATA FORMAT----------------------

local void inputdata_tlj_ascii(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, io_header_tlj_ptr hdr, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;
	real tmass;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

	input_header_tlj(instr, hdr);

	*nbody=hdr->nbody;
	*ndim = hdr->ndim;

    bodytab = (bodyptr) allocate(hdr->nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_int(instr, &Id(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_short(instr, &Type(p));

	tmass=0.0;
	DO_BODY(p, bodytab, bodytab+hdr->nbody) {
        in_real(instr, &Mass(p));
		tmass += Mass(p);
	}
	if (tmass != hdr->mass1*hdr->nbody1+hdr->mass2*hdr->nbody2+hdr->mass3*hdr->nbody3)
        error("inputdata_tlj_ascii: inconsistent values of the masses : %g %g %g\n",
			tmass, hdr->mass1*hdr->nbody1,hdr->mass2*hdr->nbody2); 

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_vector(instr, Pos(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_vector(instr, Vel(p));
    fclose(instr);
    if (scanopt(options, "reset-time"))
        hdr->tnow = 0.0;
	*tnow = hdr->tnow;
}

local void inputdata_blj_ascii(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, io_header_blj_ptr hdr, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;
	real tmass;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

	input_header_blj(instr, hdr);

	*nbody=hdr->nbody;
	*ndim = hdr->ndim;

    bodytab = (bodyptr) allocate(hdr->nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_int(instr, &Id(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_short(instr, &Type(p));

	tmass=0.0;
	DO_BODY(p, bodytab, bodytab+hdr->nbody) {
        in_real(instr, &Mass(p));
		tmass += Mass(p);
	}
	if (rabs(tmass - hdr->mass1*hdr->nbody1-hdr->mass2*hdr->nbody2)>1.0e-5)
        error("inputdata_blj_ascii: inconsistent values of the masses : %g %g %g %g\n",
			tmass, hdr->mass1*hdr->nbody1,hdr->mass2*hdr->nbody2,
			rabs(tmass - hdr->mass1*hdr->nbody1-hdr->mass2*hdr->nbody2)); 

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_vector(instr, Pos(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_vector(instr, Vel(p));
    fclose(instr);
    if (scanopt(options, "reset-time"))
        hdr->tnow = 0.0;
	*tnow = hdr->tnow;
}

local void inputdata_blj_pv(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, io_header_blj_ptr hdr, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;
	real tmass;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

	input_header_blj_pv(instr, hdr);

	*nbody=hdr->nbody;
	*ndim = hdr->ndim;

    bodytab = (bodyptr) allocate(hdr->nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+hdr->nbody) {
		in_int(instr, &Id(p));
		in_short(instr, &Type(p));
        in_real(instr, &Mass(p));               
        in_vector(instr, Pos(p));               
        in_vector(instr, Vel(p));
	}
    fclose(instr);

	tmass=0.0;
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
		tmass += Mass(p);
	if (tmass != hdr->mass1*hdr->nbody1+hdr->mass2*hdr->nbody2)
        error("inputdata_ascii: inconsistent values of the masses : %g %g %g\n",
			tmass, hdr->mass1*hdr->nbody1,hdr->mass2*hdr->nbody2); 

    if (scanopt(options, "reset-time"))
        hdr->tnow = 0.0;
	*tnow = hdr->tnow;
}

local void inputdata_blj_bin(char *file, int step, int *nbody, int *ndim, 
	realptr tnow, bool *exist_snap, io_header_blj_ptr hdr, char *options)
{
    char namebuf[256];
    struct stat buf;
    stream instr;
    bodyptr p;
	real tmass;

    sprintf(namebuf, file, step);
    if (stat(namebuf, &buf) != 0) {
		*exist_snap = FALSE;
		return;
    } else {
		*exist_snap = TRUE;
        instr = stropen(namebuf, "r");
	}

	input_header_blj_bin(instr, hdr);

	*nbody=hdr->nbody;
	*ndim = hdr->ndim;

    bodytab = (bodyptr) allocate(hdr->nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_int_bin(instr, &Id(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_short_bin(instr, &Type(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_real_bin(instr, &Mass(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_vector_bin(instr, Pos(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        in_vector_bin(instr, Vel(p));               
    fclose(instr);                              

	tmass=0.0;
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
		tmass += Mass(p);
	if (tmass != hdr->mass1*hdr->nbody1+hdr->mass2*hdr->nbody2)
        error("inputdata_ascii: inconsistent values of the masses : %g %g %g\n",
			tmass, hdr->mass1*hdr->nbody1,hdr->mass2*hdr->nbody2); 

    if (scanopt(options, "reset-time"))         
        hdr->tnow = 0.0;
	*tnow = hdr->tnow;
}

local void input_header_tlj(stream instr, io_header_tlj_ptr hdr)
{
    int ndim;

    in_int(instr, &hdr->nbody);
    if (hdr->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", hdr->nbody);
    in_int(instr, &hdr->nbody1);
    if (hdr->nbody1 < 0)
        error("inputdata: nbody1 = %d is absurd\n", hdr->nbody1);
    in_int(instr, &hdr->nbody2);
    if (hdr->nbody2 < 0)
        error("inputdata: nbody2 = %d is absurd\n", hdr->nbody2);
    in_int(instr, &hdr->nbody3);
    if (hdr->nbody3 < 0)
        error("inputdata: nbody3 = %d is absurd\n", hdr->nbody3);
	if (hdr->nbody != hdr->nbody1+hdr->nbody2+hdr->nbody3)
        error("input_header: inconsistent values of nbody, nbody1, nbody2 and nbody3 : %d %d %d %d\n", 
			hdr->nbody, hdr->nbody1, hdr->nbody2, hdr->nbody3);

    in_int(instr, &hdr->ndim);
    if (hdr->ndim != NDIM)
        error("input_header: ndim = %d; expected %d\n", hdr->ndim, NDIM);
    in_real(instr, &hdr->tnow);
    in_real(instr, &hdr->temperature);
    in_real(instr, &hdr->density);
    in_real(instr, &hdr->mass1);
    in_real(instr, &hdr->mass2);
    in_real(instr, &hdr->mass3);
    in_real(instr, &hdr->Lx);
    in_real(instr, &hdr->Ly);
#ifdef THREEDIM
    in_real(instr, &hdr->Lz);
#endif
    in_real(instr, &hdr->eps11);
    in_real(instr, &hdr->eps12);
    in_real(instr, &hdr->eps13);
    in_real(instr, &hdr->eps22);
    in_real(instr, &hdr->eps23);
    in_real(instr, &hdr->eps33);
    in_real(instr, &hdr->sigma11);
    in_real(instr, &hdr->sigma12);
    in_real(instr, &hdr->sigma13);
    in_real(instr, &hdr->sigma22);
    in_real(instr, &hdr->sigma23);
    in_real(instr, &hdr->sigma33);
    in_real(instr, &hdr->Rcut11);
    in_real(instr, &hdr->Rcut12);
    in_real(instr, &hdr->Rcut13);
    in_real(instr, &hdr->Rcut22);
    in_real(instr, &hdr->Rcut23);
    in_real(instr, &hdr->Rcut33);
}

local void input_header_blj(stream instr, io_header_blj_ptr hdr)
{
    int ndim;

    in_int(instr, &hdr->nbody);
    if (hdr->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", hdr->nbody);
    in_int(instr, &hdr->nbody1);
    if (hdr->nbody1 < 0)
        error("inputdata: nbody1 = %d is absurd\n", hdr->nbody1);
    in_int(instr, &hdr->nbody2);
    if (hdr->nbody2 < 0)
        error("inputdata: nbody2 = %d is absurd\n", hdr->nbody2);
	if (hdr->nbody != hdr->nbody1+hdr->nbody2)
        error("input_header: inconsistent values of nbody, nbody1 and nbody2 : %d %d %d\n", 
			hdr->nbody, hdr->nbody1, hdr->nbody2);

    in_int(instr, &hdr->ndim);
    if (hdr->ndim != NDIM)
        error("input_header: ndim = %d; expected %d\n", hdr->ndim, NDIM);
    in_real(instr, &hdr->tnow);
    in_real(instr, &hdr->temperature);
    in_real(instr, &hdr->density);
    in_real(instr, &hdr->mass1);
    in_real(instr, &hdr->mass2);
    in_real(instr, &hdr->Lx);
    in_real(instr, &hdr->Ly);
#ifdef THREEDIM
    in_real(instr, &hdr->Lz);
#endif
    in_real(instr, &hdr->eps11);
    in_real(instr, &hdr->eps12);
    in_real(instr, &hdr->eps22);
    in_real(instr, &hdr->sigma11);
    in_real(instr, &hdr->sigma12);
    in_real(instr, &hdr->sigma22);
    in_real(instr, &hdr->Rcut11);
    in_real(instr, &hdr->Rcut12);
    in_real(instr, &hdr->Rcut22);
}

local void input_header_blj_pv(stream instr, io_header_blj_ptr hdr)
{
    int ndim;
	char gato[1], firstline[20];

	fgets(firstline,200,instr);

	fscanf(instr,"%1s",gato);
    in_int(instr, &hdr->nbody);
    if (hdr->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", hdr->nbody);

	fscanf(instr,"%1s",gato);
    in_int(instr, &hdr->nbody1);
    if (hdr->nbody1 < 0)
        error("inputdata: nbody1 = %d is absurd\n", hdr->nbody1);

	fscanf(instr,"%1s",gato);
    in_int(instr, &hdr->nbody2);
    if (hdr->nbody2 < 0)
        error("inputdata: nbody2 = %d is absurd\n", hdr->nbody2);
	if (hdr->nbody != hdr->nbody1+hdr->nbody2)
        error("input_header: inconsistent values of nbody, nbody1 and nbody2 : %d %d %d\n", 
			hdr->nbody, hdr->nbody1, hdr->nbody2);

	fscanf(instr,"%1s",gato);
    in_int(instr, &hdr->ndim);
    if (hdr->ndim != NDIM)
        error("input_header: ndim = %d; expected %d\n", hdr->ndim, NDIM);

	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->tnow);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->temperature);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->density);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->mass1);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->mass2);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->Lx);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->Ly);
#ifdef THREEDIM
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->Lz);
#endif
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->eps11);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->eps12);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->eps22);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->sigma11);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->sigma12);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->sigma22);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->Rcut11);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->Rcut12);
	fscanf(instr,"%1s",gato);
    in_real(instr, &hdr->Rcut22);
}

local void input_header_blj_bin(stream instr, io_header_blj_ptr hdr)
{
    int ndim;

    in_int_bin(instr, &hdr->nbody);
    if (hdr->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", hdr->nbody);
    in_int_bin(instr, &hdr->nbody1);
    if (hdr->nbody1 < 0)
        error("inputdata: nbody1 = %d is absurd\n", hdr->nbody1);
    in_int_bin(instr, &hdr->nbody2);
    if (hdr->nbody2 < 0)
        error("inputdata: nbody2 = %d is absurd\n", hdr->nbody2);
	if (hdr->nbody != hdr->nbody1+hdr->nbody2)
        error("input_header: inconsistent values of nbody, nbody1 and nbody2 : %d %d %d\n", 
			hdr->nbody, hdr->nbody1, hdr->nbody2);

    in_int_bin(instr, &hdr->ndim);
    if (hdr->ndim != NDIM)
        error("input_header: ndim = %d; expected %d\n", hdr->ndim, NDIM);
    in_real_bin(instr, &hdr->tnow);
    in_real_bin(instr, &hdr->temperature);
    in_real_bin(instr, &hdr->density);
    in_real_bin(instr, &hdr->mass1);
    in_real_bin(instr, &hdr->mass2);
    in_real_bin(instr, &hdr->Lx);
    in_real_bin(instr, &hdr->Ly);
#ifdef THREEDIM
    in_real_bin(instr, &hdr->Lz);
#endif
    in_real_bin(instr, &hdr->eps11);
    in_real_bin(instr, &hdr->eps12);
    in_real_bin(instr, &hdr->eps22);
    in_real_bin(instr, &hdr->sigma11);
    in_real_bin(instr, &hdr->sigma12);
    in_real_bin(instr, &hdr->sigma22);
    in_real_bin(instr, &hdr->Rcut11);
    in_real_bin(instr, &hdr->Rcut12);
    in_real_bin(instr, &hdr->Rcut22);
}

// --------------TERMINAN I/O ROUTINES FOR BLJ DATA FORMAT----------------------

// -------------TERMINA INPUT ROUTINES -----------------------------------------

void SnapLessNBody(int *nbody, int snapcount, real reductionFac, char *options, 
	global_data_tree *gdtree, global_data_tree_bljforcecalc *gdforce)
{
    bodyptr p, bodytabsmall;
	int nbodysmall, nbodiesgroup, nbodysave;
	real tmass;

	nbodysave = *nbody;
	nbodysmall = *nbody*reductionFac;
	nbodiesgroup = 1.0/reductionFac;
	fprintf(stdout,
		"\nnbody inital-nbodysmall inital-reduction nbodiesgroup : %d %d %g %d",
		*nbody, nbodysmall, reductionFac, nbodiesgroup);
	gdtree->rsize = 1.01;
	DO_BODY(p,bodytab,bodytab+*nbody)		// Set num of bodies in each node type body
        NBodies(p) = 1;

	maketree(bodytab, *nbody, options, gdtree, gdforce);
	printf("\nncells tdepth rsize : %d %d %g\n",gdtree->ncell, gdtree->tdepth, gdtree->rsize);
	printf("\nRoot NBodies Mass Pos : %d %g %g %g %g\n",
		NBodies(gdtree->root), Mass(gdtree->root), Pos(gdtree->root)[0],Pos(gdtree->root)[1],Pos(gdtree->root)[2]);

    bodytabsmall = (bodyptr) allocate(*nbody * sizeof(body));

	scantree(bodytab, *nbody, bodytabsmall, &nbodysmall, nbodiesgroup, gdtree);

	printf("\n\nUpdating body table : %d\n", nbodysmall);
	free(bodytab);
	bodytab = bodytabsmall;
	*nbody=nbodysmall;

	tmass=0.;
	DO_BODY(p,bodytab,bodytab+*nbody)
        tmass += Mass(p);
	printf("\n\nFinal reduction factor and total mass: %g %g\n", 
		(real) nbodysmall/ (real) nbodysave, tmass);
}

local int ncellcounter;
local int nbodycounter;
local int ntotalcounter;

local void scantree(bodyptr btab, int nbod, bodyptr btabsmall, int *nbodsmall, int npgroup, global_data_tree *gdtree)
{
	bodyptr p;

	nbodycounter=0;
	ncellcounter=0;
	ntotalcounter=0;

	p = btabsmall;
	scan_walktree(p, ((nodeptr) gdtree->root), npgroup);

	printf("\nTotal-nbody npgroup nbodies ncells ntotal (bodies and cells) : %d %d %d %d %d\n",
		nbod, npgroup, nbodycounter, ncellcounter, ntotalcounter);

	*nbodsmall = nbodycounter + ncellcounter;
}

local void scan_walktree(bodyptr p, nodeptr q, int npg)
{
    nodeptr l;

	if (Type(q) == CELL) {
		if (NBodies(q) > npg) {
			for (l = More(q); l != Next(q); l = Next(l)) {
				scan_walktree(p, l, npg);
			}
		} else {
			Mass(p+ntotalcounter) = Mass(q);
			SETV(Pos(p+ntotalcounter),Pos(q));
			SETV(Vel(p+ntotalcounter),Vel(q));
			Type(p+ntotalcounter) = BODY2;
			++ncellcounter;
			++ntotalcounter;
		}
	} else {
		Mass(p+ntotalcounter) = Mass(q);
		SETV(Pos(p+ntotalcounter),Pos(q));
		SETV(Vel(p+ntotalcounter),Vel(q));
		Type(p+ntotalcounter) = BODY1;
		++nbodycounter;
		++ntotalcounter;
	}
}

void outputdata(char *file, char *outfmt, char *infmt, int snapcount, 
				int nbody, real tnow, void *hdr, char *options)
{
	int infmt_int, outfilefmt_int;
	void *hdrout;
	short allocate_mode;

	infilefmt_string_to_int(infmt, &infmt_int);
    outfilefmt_string_to_int(outfmt, &outfilefmt_int);

    if (! strnull(file) ) { 
        switch(outfilefmt_int) {
            case IO_SNAP_FMT: 
                printf("\n\tsnap-ascii format output"); 
				outputdata_ascii(file, snapcount, nbody, tnow, hdr, options); break;
// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
            case IO_SNAP_FMT_LONG:
                printf("\n\tsnap-ascii-long format output"); 
				outputdata_ascii_long(file, snapcount, nbody, tnow, hdr, options); break;
//

            case IO_NULL_FMT: 
                printf("\n\tsnap-ascii format output"); 
				outputdata_ascii(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_PV_FMT: 
                printf("\n\tpv format output"); 
				outputpvdata(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_SNAP_FMT_BIN: 
                printf("\n\tsnap-bin format output"); 
				outputbindata(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_GADGET11_FMT_NORMAL_BODY_ASCII: 
                printf("\n\tgadget-normal-body-ascii format output\n"); 
				save_ic_gdgi_format_normal_body(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_GADGET11_FMT_NORMAL_BODY_ASCII_LONG: 
                printf("\n\tgadget-normal-body-ascii-long format output\n"); 
				outputdata_gadget11_normal_body_ascii_long(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_GADGET11_FMT_SPH_ASCII: 
                printf("\n\tgadget-sph-ascii format output"); 
				save_ic_gdgi_format_sph(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_GADGET11_FMT_BIN_DOUBLE:
                printf("\n\tgadget11-bin-double format output"); 
				outputdata_gadget11_bin_double(file, snapcount, nbody, tnow, hdr, options); break;
// GADGET207 COMIENZO //////////////////////////////////////////////////////////
//            case IO_GADGET207_FMT_BIN:
//                printf("\n\tgadget207-bin format output"); 
//				outputdata_gadget207_bin(file, snapcount, nbody, tnow, hdr, options); break;
// GADGET207 FIN ///////////////////////////////////////////////////////////
                
// IBERO COMIENZO //////////////////////////////////////////////////////////
//            case IO_GADGET11_FMT_BIN_IBERO:
//                printf("\n\tgadget11-bin-ibero format output"); 
//				outputdata_gadget11_bin_ibero(file, snapcount, nbody, tnow, hdr, options); break;
// IBERO FIN ///////////////////////////////////////////////////////////
                
            case IO_GADGET11_FMT_ASCII:
                printf("\n\tgadget11-ascii format output"); 
				outputdata_gadget11_ascii(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_GADGET11_FMT_BIN_DOUBLE_REDUCIDO:
                printf("\n\tgadget11-bin-double-reducido format output"); 
				outputdata_gadget11_bin_double_reducido(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_GADGET11_FMT_ASCII_REDUCIDO:
                printf("\n\tgadget11-ascii-reducido format output");
				allocate_mode = 4;
				outputdata_gadget11_ascii_reducido(file, snapcount, nbody, 
					tnow, hdr, options, allocate_mode); break;
            case IO_TIPSY_FMT_BIN:
                printf("\n\ttipsy-bin format output"); 
				outputdata_tipsy_bin(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_POWMES_FMT_ASCII:
                printf("\n\tpowmes-ascii format output"); 
				outputdata_powmes_ascii(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_POWMES_FMT_ASCII_LONG:
                fprintf(stdout,"\n\tpowmes-ascii-long format output");
				outputdata_powmes_ascii_long(file, snapcount, nbody, tnow, hdr, options); break;
// BEGIN NEMO REMOVE :: 2019-05-14
//            case IO_NEMO_FMT:
//                printf("\n\tnemo format output");
//				outputdata_nemo(file, snapcount, nbody, tnow, hdr, options); break;
// END
            case IO_SNAP_BLJ_FMT:
                printf("\n\tsnap-blj-ascii format output"); 
				if (infmt_int == outfilefmt_int)
					outputdata_blj(file, snapcount, nbody, tnow, hdr, options); 
				else {
					hdrout = (io_header_blj *) allocate(sizeof(io_header_blj));
					HeaderConversion_blj(infmt_int,outfilefmt_int,hdr,hdrout);
					outputdata_blj(file,snapcount,nbody,tnow,hdrout,options); 
				}
				break;
            case IO_SNAP_BLJ_PV_FMT:
                printf("\n\tsnap-blj-pv format output"); 
				outputpvdata_blj(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_SNAP_BLJ_BIN_FMT:
                printf("\n\tsnap-blj-bin format output"); 
				outputbindata_blj(file, snapcount, nbody, tnow, hdr, options); break;
            case IO_SNAP_TLJ_FMT:
                printf("\n\tsnap-tlj-ascii format output");
				if (infmt_int == outfilefmt_int)
					outputdata_tlj(file,snapcount,nbody,tnow,hdr,options);
				else {
					hdrout = (io_header_tlj *) allocate(sizeof(io_header_tlj));
					HeaderConversion_tlj(infmt_int,outfilefmt_int,hdr,hdrout);
					outputdata_tlj(file,snapcount,nbody,tnow,hdrout,options);
				}
				break;
            default:
                printf("\n\toutput: Unknown output format...");
                printf("\n\tprinting in default snap-ascii format..."); 
				outputdata_ascii(file, snapcount, nbody, tnow, hdr, options); break;
        }
    } else
		error("\n\noutputdata : You should give an output file name\n");
}


local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"snap-ascii") == 0)
		*outfmt_int = IO_SNAP_FMT;
// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
    if (strcmp(outfmt_str,"snap-ascii-long") == 0)
		*outfmt_int = IO_SNAP_FMT_LONG;
//
    if (strnull(outfmt_str))
		*outfmt_int = IO_NULL_FMT;
    if (strcmp(outfmt_str,"snap-pv") == 0)
		*outfmt_int = IO_PV_FMT;
    if (strcmp(outfmt_str,"snap-bin") == 0)
		*outfmt_int = IO_SNAP_FMT_BIN;
    if (strcmp(outfmt_str,"gadget11-normal-body-ascii") == 0)
		*outfmt_int = IO_GADGET11_FMT_NORMAL_BODY_ASCII;
    if (strcmp(outfmt_str,"gadget11-normal-body-ascii-long") == 0)
		*outfmt_int = IO_GADGET11_FMT_NORMAL_BODY_ASCII_LONG;
    if (strcmp(outfmt_str,"gadget11-sph-ascii") == 0)
		*outfmt_int = IO_GADGET11_FMT_SPH_ASCII;
    if (strcmp(outfmt_str,"gadget11-bin-double") == 0)
		*outfmt_int = IO_GADGET11_FMT_BIN_DOUBLE;

// GADGET207 COMIENZO //////////////////////////////////////////////////////////
//    if (strcmp(outfmt_str,"gadget207-bin") == 0)
//		*outfmt_int = IO_GADGET207_FMT_BIN;
// GADGET207 FIN ///////////////////////////////////////////////////////////
    
// IBERO COMIENZO //////////////////////////////////////////////////////////
//    if (strcmp(outfmt_str,"gadget11-bin-ibero") == 0)
//		*outfmt_int = IO_GADGET11_FMT_BIN_IBERO;
// IBERO FIN ///////////////////////////////////////////////////////////
    
    if (strcmp(outfmt_str,"gadget11-ascii") == 0)
		*outfmt_int = IO_GADGET11_FMT_ASCII;
    if (strcmp(outfmt_str,"tipsy-bin") == 0)
		*outfmt_int = IO_TIPSY_FMT_BIN;
    if (strcmp(outfmt_str,"powmes-ascii") == 0)
		*outfmt_int = IO_POWMES_FMT_ASCII;
    if (strcmp(outfmt_str,"powmes-ascii-long") == 0)
		*outfmt_int = IO_POWMES_FMT_ASCII_LONG;
// BEGIN NEMO REMOVE :: 2019-05-14
//    if (strcmp(outfmt_str,"nemo") == 0)
//		*outfmt_int = IO_NEMO_FMT;
//
    if (strcmp(outfmt_str,"gadget11-bin-double-reducido") == 0)
		*outfmt_int = IO_GADGET11_FMT_BIN_DOUBLE_REDUCIDO;
    if (strcmp(outfmt_str,"gadget11-ascii-reducido") == 0)
		*outfmt_int = IO_GADGET11_FMT_ASCII_REDUCIDO;
    if (strcmp(outfmt_str,"snap-blj-ascii") == 0)
		*outfmt_int = IO_SNAP_BLJ_FMT;
    if (strcmp(outfmt_str,"snap-blj-pv") == 0)
		*outfmt_int = IO_SNAP_BLJ_PV_FMT;
    if (strcmp(outfmt_str,"snap-blj-bin") == 0)
		*outfmt_int = IO_SNAP_BLJ_BIN_FMT;
    if (strcmp(outfmt_str,"snap-tlj-ascii") == 0) 
		*outfmt_int = IO_SNAP_TLJ_FMT;
}


#define SNAPTMP			"snap.tmp"

local void outputdata_ascii(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, file, snapcount);           
    outstr = stropen(namebuf, "w!");
    out_int(outstr, nbody);
    out_int(outstr, NDIM);
    out_real(outstr, tnow);

fprintf(stdout,"\nWriting %d particle masses ... ",nbody);
    for (p = bodytab; p < bodytab+nbody; p++)
        out_real(outstr, Mass(p));
fprintf(stdout,"done. Counter %ld",p-bodytab+1);

fprintf(stdout,"\nWriting %d particle positions ... ",nbody);
    for (p = bodytab; p < bodytab+nbody; p++)
        out_vector(outstr, Pos(p));
fprintf(stdout,"done. Counter %d",p-bodytab+1);

fprintf(stdout,"\nWriting %d particle velocities ... ",nbody);
    for (p = bodytab; p < bodytab+nbody; p++)
        out_vector(outstr, Vel(p));
fprintf(stdout,"done. Counter %d",p-bodytab+1);

    if (scanopt(options, "out-phi"))            
        for (p = bodytab; p < bodytab+nbody; p++)
            out_real(outstr, Phi(p));           
    if (scanopt(options, "out-acc"))            
        for (p = bodytab; p < bodytab+nbody; p++)
            out_vector(outstr, Acc(p));  
    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
}

// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
local void outputdata_ascii_long(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr_long p;

    sprintf(namebuf, file, snapcount);           
    outstr = stropen(namebuf, "w!");
    out_int(outstr, nbody);
    out_int(outstr, NDIM);
    out_real(outstr, tnow);

fprintf(stdout,"\nWriting %d particle masses ... ",nbody);
    for (p = bodytab_long; p < bodytab_long+nbody; p++) //{
        out_real(outstr, Mass_long(p));
fprintf(stdout,"done. Counter %ld",p-bodytab_long+1);

fprintf(stdout,"\nWriting %d particle positions ... ",nbody);
    for (p = bodytab_long; p < bodytab_long+nbody; p++)
        out_vector(outstr, Pos_long(p));
fprintf(stdout,"done. Counter %d",p-bodytab_long+1);

fprintf(stdout,"\nWriting %d particle velocities ... ",nbody);
    for (p = bodytab_long; p < bodytab_long+nbody; p++)
        out_vector(outstr, Vel_long(p));
fprintf(stdout,"done. Counter %d",p-bodytab_long+1);

    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
}
//


local void outputbindata(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, file, snapcount);
    outstr = stropen(namebuf, "w!");         
    out_int_bin(outstr, nbody);
    out_int_bin(outstr, NDIM);                      
    out_real_bin(outstr, tnow);                     
    for (p = bodytab; p < bodytab+nbody; p++)   
        out_real_bin(outstr, Mass(p));              
    for (p = bodytab; p < bodytab+nbody; p++)
        out_vector_bin(outstr, Pos(p));             
    for (p = bodytab; p < bodytab+nbody; p++)
        out_vector_bin(outstr, Vel(p));             
    if (scanopt(options, "out-phi"))           
        for (p = bodytab; p < bodytab+nbody; p++)
            out_real_bin(outstr, Phi(p));           
    if (scanopt(options, "out-acc"))            
        for (p = bodytab; p < bodytab+nbody; p++)
            out_vector_bin(outstr, Acc(p));         
    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
}

local void outputpvdata(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;


    sprintf(namebuf, file, snapcount);
    outstr = stropen(namebuf, "w!");         
    for (p = bodytab; p < bodytab+nbody; p++) {
		out_int_mar(outstr, (p-bodytab+1));

        out_vector_mar(outstr, Pos(p));
        out_vector_mar(outstr, Vel(p));

        if (scanopt(options, "out-phi"))
            out_real_mar(outstr, Phi(p));
        if (scanopt(options, "out-acc"))
            out_vector_mar(outstr, Acc(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);                         
    printf("\tpos-vel data output to file %s at time %f\n\n",namebuf,tnow);
}

local void outputdata_tipsy_bin(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    int ndim ;
    int nbodies ;
    int ngas ;
    int ndark ;
    int nstar ;

    gas_particle_ptr gp;
    dark_particle_ptr dp;
    star_particle_ptr sp;

    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, file, snapcount);           
    outstr = stropen(namebuf, "w!");

	header.time = tnow;
	header.nbodies = nbody;
	header.ndim = NDIM;
	header.nsph = 0;
	header.nstar = 0;

	ndim=header.ndim;
	nbodies=header.nbodies;
	ngas=header.nsph;
	nstar = header.nstar ;
	ndark = header.ndark = nbodies - nstar - ngas ;

	if(gas_particles != NULL) free(gas_particles);
	if(ngas != 0) {
	    gas_particles =
		(gas_particle_ptr) malloc(ngas*sizeof(gas_particle));
	    if(gas_particles == NULL)
			error("<sorry, no memory for gas particles, master>\n");
	}
	else
	  gas_particles = NULL;
	if(dark_particles != NULL) free(dark_particles);
	if(ndark != 0) {
	    dark_particles =
		(dark_particle_ptr) malloc(ndark*sizeof(dark_particle));
	    if(dark_particles == NULL)
			error("<sorry, no memory for dark particles, master>\n");
	}
	else
	  dark_particles = NULL;
	if(star_particles != NULL) free(star_particles);
	if(nstar != 0) {
	    star_particles =
		 (star_particle_ptr)malloc(nstar*sizeof(star_particle));
	    if(star_particles == NULL)
			error("<sorry, no memory for star particles, master>\n");
	}
	else
	  star_particles = NULL;

    for (p = bodytab, dp=dark_particles; p < bodytab+nbody; p++, dp++) {
		dp->mass = Mass(p);
		dp->pos[0] = Pos(p)[0];
		dp->pos[1] = Pos(p)[1];
		dp->pos[2] = Pos(p)[2];
		dp->vel[0] = Vel(p)[0];
		dp->vel[1] = Vel(p)[1];
		dp->vel[2] = Vel(p)[2];
		dp->eps = 0.;
		dp->phi = 0.;
	}

	fwrite((char *)&header,sizeof(header),1, outstr);
	fwrite((char *)gas_particles,sizeof(gas_particle),ngas, outstr);
	fwrite((char *)dark_particles,sizeof(dark_particle),ndark, outstr);
	fwrite((char *)star_particles,sizeof(star_particle),nstar, outstr);

    fclose(outstr);
    printf("\n\tdata output (tipsy bin format) to file %s at time %f\n", namebuf, tnow);
}

local void outputdata_powmes_ascii(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;
	int k;

// Normalizacion de las coordenadas ... esta bajo desarrollo ...
// Solo si se conoce el L (BoxSize) podremos hacerla...
    printf("\n\tSize of the box is (from gadget file) %f\n", header1.BoxSize);

    sprintf(namebuf, file, snapcount);
    outstr = stropen(namebuf, "w!");
	fprintf(outstr,"%d\n",nbody);
    for (p = bodytab; p < bodytab+nbody; p++) {
		for (k=0; k<NDIM; k++)
			Pos(p)[k] /= header1.BoxSize;
        out_vector_mar(outstr, Pos(p));
		out_real_mar(outstr, Mass(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);
    printf("\n\tdata output (powmes ascii format) to file %s at time %f\n", namebuf, tnow);
}

local void outputdata_powmes_ascii_long(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr_long p;
	int k;

    fprintf(stdout,"\n\tdata output (powmes-ascii-long format) ... begining\n");
	fflush(stdout);

// Normalizacion de las coordenadas ... esta bajo desarrollo ...
// Solo si se conoce el L (BoxSize) podremos hacerla...
    printf("\n\tSize of the box is (from gadget file) %f\n", header1.BoxSize);


    sprintf(namebuf, file, snapcount);
    outstr = stropen(namebuf, "w!");
	fprintf(outstr,"%ld\n",nbody);
    for (p = bodytab_long; p < bodytab_long+nbody; p++) {
		for (k=0; k<NDIM; k++)
			Pos_long(p)[k] /= header1.BoxSize;
        out_vector_mar(outstr, Pos_long(p));
		out_real_mar(outstr, Mass_long(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);                         
    printf("\n\tdata output (powmes-ascii-long format ) to file %s at time %f\n", namebuf, tnow);
}

// BEGIN NEMO REMOVE :: 2019-05-14
/*
local void outputdata_nemo(char *file, int snapcount,
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
	int bits;
	int seed=0;
    char    hisline[80];
    char namebuf[256];
    stream outstr;
	string headline;
	Body *btab;
	bodyptr p;
	Body *dp;

    sprintf(namebuf, file, snapcount);           
    outstr = stropen(namebuf, "w!");

    btab = (Body *) allocate (nbody * sizeof(Body));

    for (p = bodytab, dp=btab; p < bodytab+nbody; p++, dp++) {
		NBMass(dp) = Mass(p);
		NBPos(dp)[0] = Pos(p)[0];
		NBPos(dp)[1] = Pos(p)[1];
		NBPos(dp)[2] = Pos(p)[2];
		NBVel(dp)[0] = Vel(p)[0];
		NBVel(dp)[1] = Vel(p)[1];
		NBVel(dp)[2] = Vel(p)[2];
	}

	sprintf(hisline,"init_xrandom: seed used %d",seed);
	bits = (MassBit | PhaseSpaceBit | TimeBit);
	put_string(outstr, HeadlineTag, hisline);
	put_history (outstr);           // update history
	if (*headline)
		put_string (outstr, HeadlineTag, headline);
	put_snap (outstr, &btab, &nbody, &tnow, &bits);

    fclose(outstr);
}
*/
// END

local void outputdata_gadget11_bin_double(char *file, int num, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
  FILE *fd;
  char buf[100];
  double dummy[3];
  int i,k;
  int   blklen,masscount;
  double a3inv;

    char namebuf[256];

//#define BLKLEN gdgt_fwrite(&blklen, sizeof(blklen), 1, fd);

  if(All.ComovingIntegrationOn)
    a3inv=  1/(All.Time*All.Time*All.Time);
  else
    a3inv=  1.0;

    sprintf(namebuf, file, num);
  sprintf(buf,"%s%s",All.OutputDir,namebuf);

  if((fd=fopen(buf,"w")))
    {
      header1.npart[0]= header1.npartTotal[0]= All.TotN_gas;
      header1.npart[1]= header1.npartTotal[1]= All.TotN_halo;
      header1.npart[2]= header1.npartTotal[2]= All.TotN_disk;
      header1.npart[3]= header1.npartTotal[3]= All.TotN_bulge;
      header1.npart[4]= header1.npartTotal[4]= All.TotN_stars;
      header1.npart[5]= header1.npartTotal[5]= All.TotN_dm;

      for(i=0;i<NumParticleTypes;i++)
	header1.mass[i]=0;

      for(i=0, masscount=0; i<NumParticleTypes; i++)
	{
	  header1.mass[i]= All.MassTable[i];
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    masscount+= header1.npart[i];
	}

      header1.time= All.Time;

      if(All.ComovingIntegrationOn)
	header1.redshift=1.0/All.Time - 1.0;
      else
	header1.redshift=0;  
      
      header1.flag_sfr=0;
      header1.flag_feedback=0;
      header1.flag_cooling= 0;			// Hacer que funciones con cooling off
      header1.num_files= 1;
      header1.BoxSize= All.BoxSize;
      header1.Omega0=  All.Omega0;
      header1.OmegaLambda= All.OmegaLambda;
      header1.HubbleParam= All.HubbleParam;
      
      blklen=sizeof(header1);
      BLKLEN;
      gdgt_fwrite(&header1, sizeof(header1), 1, fd);
      BLKLEN;

      blklen=NumPart*3*sizeof(double);	

      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P[i].PosPred[k];
	  gdgt_fwrite(dummy,sizeof(double),3,fd);	
	}
      BLKLEN;

      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P[i].VelPred[k];
	  gdgt_fwrite(dummy,sizeof(double),3,fd);
	}
      BLKLEN;
 
      blklen=NumPart*sizeof(int);
      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fwrite(&P[i].ID,sizeof(int),1,fd);
	}
      BLKLEN;

      blklen=masscount*sizeof(double);
      if(masscount)
	BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  dummy[0]= P[i].Mass;
	  if(All.MassTable[P[i].Type]==0)
	    gdgt_fwrite(dummy,sizeof(double),1,fd);		
	}
      if(masscount)
	BLKLEN;

      if(N_gas)
	{
	  blklen=N_gas*sizeof(double);
	  BLKLEN;
	  for(i=1;i<=N_gas;i++)
	    {
	      dummy[0]=SphP[i].EgySpecPred;
	      gdgt_fwrite(dummy,sizeof(double),1,fd);	
	    }
	  BLKLEN;

	  blklen=N_gas*sizeof(double);  
	  BLKLEN;
	  for(i=1;i<=N_gas;i++)
	    {
	      dummy[0]=SphP[i].DensityPred;
	      gdgt_fwrite(dummy,sizeof(double),1,fd);	
	    }
	  BLKLEN;
	}
      fclose(fd);
    }
  else
    {
      fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
      endrun(10);
    }
}



// GADGET207 COMIENZO //////////////////////////////////////////////////////////
/*
local void outputdata_gadget207_bin(char *file, int num, 
                                         int nbody, real tnow, io_header_blj *hdr, char *options)
{
    FILE *fd;
    char buf[100];
    double dummy[3];
    int i,k;
    int   blklen,masscount;
    double a3inv;

    char namebuf[256];
    
    //#define BLKLEN gdgt_fwrite(&blklen, sizeof(blklen), 1, fd);
    
    if(All_GADGET.ComovingIntegrationOn)
        a3inv=  1/(All_GADGET.Time*All_GADGET.Time*All_GADGET.Time);
    else
        a3inv=  1.0;
    
    sprintf(namebuf, file, num);
    sprintf(buf,"%s%s",All_GADGET.OutputDir,namebuf);
    
    if((fd=fopen(buf,"w")))
    {
        header_GADGET.npart[0]= header_GADGET.npartTotal[0]= All_GADGET.TotN_gas;
        header_GADGET.npart[1]= header_GADGET.npartTotal[1]= All_GADGET.TotN_halo;
        header_GADGET.npart[2]= header_GADGET.npartTotal[2]= All_GADGET.TotN_disk;
        header_GADGET.npart[3]= header_GADGET.npartTotal[3]= All_GADGET.TotN_bulge;
        header_GADGET.npart[4]= header_GADGET.npartTotal[4]= All_GADGET.TotN_stars;
        header_GADGET.npart[5]= header_GADGET.npartTotal[5]= All_GADGET.TotN_dm;
        
        for(i=0;i<NumParticleTypes;i++)
            header_GADGET.mass[i]=0;
        
        for(i=0, masscount=0; i<NumParticleTypes; i++)
        {
            header_GADGET.mass[i]= All_GADGET.MassTable[i];
            if(All_GADGET.MassTable[i]==0 && header_GADGET.npart[i]>0)
                masscount+= header_GADGET.npart[i];
        }
        
        header_GADGET.time= All_GADGET.Time;
        
        if(All_GADGET.ComovingIntegrationOn)
            header_GADGET.redshift=1.0/All_GADGET.Time - 1.0;
        else
            header_GADGET.redshift=0;  
        
        header_GADGET.flag_sfr=0;
        header_GADGET.flag_feedback=0;
        header_GADGET.flag_cooling= 0;			// Hacer que funciones con cooling off
        header_GADGET.num_files= 1;
        header_GADGET.BoxSize= All_GADGET.BoxSize;
        header_GADGET.Omega0=  All_GADGET.Omega0;
        header_GADGET.OmegaLambda= All_GADGET.OmegaLambda;
        header_GADGET.HubbleParam= All_GADGET.HubbleParam;
        
        blklen=sizeof(header_GADGET);
        BLKLEN;
        gdgt_fwrite(&header_GADGET, sizeof(header_GADGET), 1, fd);
        BLKLEN;
        
        blklen=NumPart*3*sizeof(double);	
        
        BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            for(k=0;k<3;k++)
                dummy[k]=P_GADGET[i].PosPred[k];
            gdgt_fwrite(dummy,sizeof(double),3,fd);	
        }
        BLKLEN;
        
        BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            for(k=0;k<3;k++)
                dummy[k]=P_GADGET[i].VelPred[k];
            gdgt_fwrite(dummy,sizeof(float),3,fd);
        }
        BLKLEN;
        
        blklen=NumPart*sizeof(int);
        BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            gdgt_fwrite(&P_GADGET[i].ID,sizeof(int),1,fd);
        }
        BLKLEN;
        
        blklen=masscount*sizeof(double);
        if(masscount)
            BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            dummy[0]= P_GADGET[i].Mass;
            if(All_GADGET.MassTable[P_GADGET[i].Type]==0)
                gdgt_fwrite(dummy,sizeof(float),1,fd);		
        }
        if(masscount)
            BLKLEN;
        
        if(N_gas)
        {
            blklen=N_gas*sizeof(double);
            BLKLEN;
            for(i=1;i<=N_gas;i++)
            {
                dummy[0]=SphP[i].EgySpecPred;
                gdgt_fwrite(dummy,sizeof(double),1,fd);	
            }
            BLKLEN;
            
            blklen=N_gas*sizeof(double);  
            BLKLEN;
            for(i=1;i<=N_gas;i++)
            {
                dummy[0]=SphP[i].DensityPred;
                gdgt_fwrite(dummy,sizeof(double),1,fd);	
            }
            BLKLEN;
        }
        fclose(fd);
    }
    else
    {
        fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
        endrun(10);
    }
}
*/
// GADGET207 FIN ///////////////////////////////////////////////////////////


// IBERO COMIENZO //////////////////////////////////////////////////////////
/*
local void outputdata_gadget11_bin_ibero(char *file, int num, 
                                          int nbody, real tnow, io_header_blj *hdr, char *options)
{
    FILE *fd;
    char buf[100];
    double dummy[3];
    int i,k;
    int   blklen,masscount;
    double a3inv;
    
    char namebuf[256];
    
    //#define BLKLEN gdgt_fwrite(&blklen, sizeof(blklen), 1, fd);
    
    if(All.ComovingIntegrationOn)
        a3inv=  1/(All.Time*All.Time*All.Time);
    else
        a3inv=  1.0;
    
    sprintf(namebuf, file, num);
    sprintf(buf,"%s%s",All.OutputDir,namebuf);
    
    if((fd=fopen(buf,"w")))
    {
        header1_IBERO.npart[0]= header1_IBERO.npartTotal[0]= All.TotN_gas;
        header1_IBERO.npart[1]= header1_IBERO.npartTotal[1]= All.TotN_halo;
        header1_IBERO.npart[2]= header1_IBERO.npartTotal[2]= All.TotN_disk;
        header1_IBERO.npart[3]= header1_IBERO.npartTotal[3]= All.TotN_bulge;
        header1_IBERO.npart[4]= header1_IBERO.npartTotal[4]= All.TotN_stars;
        header1_IBERO.npart[5]= header1_IBERO.npartTotal[5]= All.TotN_dm;
        
        for(i=0;i<NumParticleTypes;i++)
            header1_IBERO.mass[i]=0;
        
        for(i=0, masscount=0; i<NumParticleTypes; i++)
        {
            header1_IBERO.mass[i]= All.MassTable[i];
            if(All.MassTable[i]==0 && header1_IBERO.npart[i]>0)
                masscount+= header1_IBERO.npart[i];
        }
        
        header1_IBERO.time= All.Time;
        
        if(All.ComovingIntegrationOn)
            header1_IBERO.redshift=1.0/All.Time - 1.0;
        else
            header1_IBERO.redshift=0;  
        
        header1_IBERO.flag_sfr=0;
        header1_IBERO.flag_feedback=0;
        header1_IBERO.flag_cooling= 0;			// Hacer que funciones con cooling off
        header1_IBERO.num_files= 1;
        header1_IBERO.BoxSize= All.BoxSize;
        header1_IBERO.Omega0=  All.Omega0;
        header1_IBERO.OmegaLambda= All.OmegaLambda;
        header1_IBERO.HubbleParam= All.HubbleParam;
        
        blklen=sizeof(header1_IBERO);
        BLKLEN;
        gdgt_fwrite(&header1_IBERO, sizeof(header1_IBERO), 1, fd);
        BLKLEN;
        
        blklen=NumPart*3*sizeof(double);	
        
        BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            for(k=0;k<3;k++)
                dummy[k]=P_IBERO[i].PosPred[k];
            gdgt_fwrite(dummy,sizeof(double),3,fd);	
        }
        BLKLEN;
        
        BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            for(k=0;k<3;k++)
                dummy[k]=P_IBERO[i].VelPred[k];
            gdgt_fwrite(dummy,sizeof(float),3,fd);
        }
        BLKLEN;
        
        blklen=NumPart*sizeof(int);
        BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            gdgt_fwrite(&P_IBERO[i].ID,sizeof(int),1,fd);
        }
        BLKLEN;
        
        blklen=masscount*sizeof(double);
        if(masscount)
            BLKLEN;
        for(i=1;i<=NumPart;i++)
        {
            dummy[0]= P_IBERO[i].Mass;
            if(All.MassTable[P_IBERO[i].Type]==0)
                gdgt_fwrite(dummy,sizeof(float),1,fd);		
        }
        if(masscount)
            BLKLEN;
        
        if(N_gas)
        {
            blklen=N_gas*sizeof(double);
            BLKLEN;
            for(i=1;i<=N_gas;i++)
            {
                dummy[0]=SphP[i].EgySpecPred;
                gdgt_fwrite(dummy,sizeof(double),1,fd);	
            }
            BLKLEN;
            
            blklen=N_gas*sizeof(double);  
            BLKLEN;
            for(i=1;i<=N_gas;i++)
            {
                dummy[0]=SphP[i].DensityPred;
                gdgt_fwrite(dummy,sizeof(double),1,fd);	
            }
            BLKLEN;
        }
        fclose(fd);
    }
    else
    {
        fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
        endrun(10);
    }
}
*/
// IBERO FIN ///////////////////////////////////////////////////////////



local void outputdata_gadget11_ascii(char *file, int num, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
  FILE *fd;
  char buf[100];
  double dummy[3];
  int i,k;
  int   blklen,masscount;
  double a3inv;
  vector dummyvec;
  double realdummy;

    char namebuf[256];


  if(All.ComovingIntegrationOn)
    a3inv=  1/(All.Time*All.Time*All.Time);
  else
    a3inv=  1.0;

    sprintf(namebuf, file, num);
  sprintf(buf,"%s%s",All.OutputDir,namebuf);

  if((fd=fopen(buf,"w")))
    {
      header1.npart[0]= header1.npartTotal[0]= All.TotN_gas;
      header1.npart[1]= header1.npartTotal[1]= All.TotN_halo;
      header1.npart[2]= header1.npartTotal[2]= All.TotN_disk;
      header1.npart[3]= header1.npartTotal[3]= All.TotN_bulge;
      header1.npart[4]= header1.npartTotal[4]= All.TotN_stars;
      header1.npart[5]= header1.npartTotal[5]= 0;

      for(i=0;i<6;i++)			
	header1.mass[i]=0;

      for(i=0, masscount=0; i<5; i++)	
	{
	  header1.mass[i]= All.MassTable[i];
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    masscount+= header1.npart[i];
	}

      header1.time= All.Time;

      if(All.ComovingIntegrationOn)
	header1.redshift=1.0/All.Time - 1.0;
      else
	header1.redshift=0;  

      
      header1.flag_sfr=0;
      header1.flag_feedback=0;
      header1.flag_cooling= 0;		// Hacer que funcione con cooling off
      header1.num_files= 1;
      header1.BoxSize= All.BoxSize;
      header1.Omega0=  All.Omega0;
      header1.OmegaLambda= All.OmegaLambda;
      header1.HubbleParam= All.HubbleParam;
      
	fprint_header(fd);

      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummyvec[k]=P[i].Pos[k];
        out_vector(fd, dummyvec);
	}

      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummyvec[k]=P[i].Vel[k];
        out_vector(fd, dummyvec);
	}

      for(i=1;i<=NumPart;i++)
	{
		out_int(fd,P[i].ID);
	}

printf("\n\nmasscount & N_gas : %d %d\n",masscount, N_gas);

      for(i=1;i<=NumPart;i++)
	{
	  realdummy= P[i].Mass;
	  if(All.MassTable[P[i].Type]==0)
        out_real(fd, realdummy);
	}

      if(N_gas)
	{
	  for(i=1;i<=N_gas;i++)
	    {
	      dummy[0]=SphP[i].EgySpec;
		  out_real(fd,dummy[0]);
	    }

	}
      fclose(fd);
    }
  else
    {
      fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
      endrun(10);
    }
}

#define IFMT  " %d"
#define RFMT  " %14.7E" 

local void save_ic_gdgi_format_normal_body(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
	bodyptr p;

    sprintf(namebuf, file, snapcount);
    outstr = stropen(namebuf, "w!");

	header1.npart[0]=0;
	header1.npart[1]=nbody;
	header1.npart[2]=0;
	header1.npart[3]=0;
	header1.npart[4]=0;
	header1.npart[5]=0;
	header1.mass[0]=0.0;
	header1.mass[1]=0.0;
	header1.mass[2]=0.0;
	header1.mass[3]=0.0;
	header1.mass[4]=0.0;
	header1.mass[5]=0.0;
	header1.time=0.0;
	header1.redshift=0.0;
	header1.flag_sfr=0;
	header1.flag_feedback=0;
	header1.npartTotal[0]=0;
	header1.npartTotal[1]=nbody;
	header1.npartTotal[2]=0;
	header1.npartTotal[3]=0;
	header1.npartTotal[4]=0;
	header1.npartTotal[5]=0;
	header1.flag_cooling=0;
	header1.num_files=1;
	header1.BoxSize=0.0;
	header1.Omega0=0.0;
	header1.OmegaLambda=0.0;
	header1.HubbleParam=0.0;
	fprint_header(outstr);
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr,RFMT RFMT RFMT "\n",Pos(p)[0],Pos(p)[1],Pos(p)[2]);
	}
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr,RFMT RFMT RFMT "\n",Vel(p)[0],Vel(p)[1],Vel(p)[2]);
	}
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr,IFMT "\n",p-bodytab+1);
	}
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr, RFMT "\n",Mass(p));
	}
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
	fclose(outstr);
}

local void outputdata_gadget11_normal_body_ascii_long(char *file, int snapcount, 
							int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
	bodyptr_long p;

    sprintf(namebuf, file, snapcount);
    outstr = stropen(namebuf, "w!");
	
	header1.npart[0]=0;
	header1.npart[1]=nbody;
	header1.npart[2]=0;
	header1.npart[3]=0;
	header1.npart[4]=0;
	header1.npart[5]=0;
	header1.mass[0]=0.0;
	header1.mass[1]=0.0;
	header1.mass[2]=0.0;
	header1.mass[3]=0.0;
	header1.mass[4]=0.0;
	header1.mass[5]=0.0;
	header1.time=0.0;
	header1.redshift=0.0;
	header1.flag_sfr=0;
	header1.flag_feedback=0;
	header1.npartTotal[0]=0;
	header1.npartTotal[1]=nbody;
	header1.npartTotal[2]=0;
	header1.npartTotal[3]=0;
	header1.npartTotal[4]=0;
	header1.npartTotal[5]=0;
	header1.flag_cooling=0;
	header1.num_files=1;
	header1.BoxSize=0.0;
	header1.Omega0=0.0;
	header1.OmegaLambda=0.0;
	header1.HubbleParam=0.0;
	fprint_header(outstr);
	for (p = bodytab_long; p < bodytab_long+nbody; p++) { 
		fprintf(outstr,RFMT RFMT RFMT "\n",Pos_long(p)[0],Pos_long(p)[1],Pos_long(p)[2]);
	}
	for (p = bodytab_long; p < bodytab_long+nbody; p++) { 
		fprintf(outstr,RFMT RFMT RFMT "\n",Vel_long(p)[0],Vel_long(p)[1],Vel_long(p)[2]);
	}
	for (p = bodytab_long; p < bodytab_long+nbody; p++) { 
		fprintf(outstr,IFMT "\n",p-bodytab_long+1);
	}
	for (p = bodytab_long; p < bodytab_long+nbody; p++) { 
		fprintf(outstr, RFMT "\n",Mass_long(p));
	}
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
	fclose(outstr);
}

local void save_ic_gdgi_format_sph(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;
	real InternalEnergy;
	real SoundSpeed;

    sprintf(namebuf, file, snapcount);
    outstr = stropen(namebuf, "w!");

	SoundSpeed = 0.369331351;	// Checar la consistencia de este valor con la
								//	condicion inicial en general

	header1.npart[0]=nbody;
	header1.npart[1]=0;
	header1.npart[2]=0;
	header1.npart[3]=0;
	header1.npart[4]=0;
	header1.npart[5]=0;
	header1.mass[0]=0.0;
	header1.mass[1]=0.0;
	header1.mass[2]=0.0;
	header1.mass[3]=0.0;
	header1.mass[4]=0.0;
	header1.mass[5]=0.0;
	header1.time=0.0;
	header1.redshift=0.0;
	header1.flag_sfr=0;
	header1.flag_feedback=0;
	header1.npartTotal[0]=nbody;
	header1.npartTotal[1]=0;
	header1.npartTotal[2]=0;
	header1.npartTotal[3]=0;
	header1.npartTotal[4]=0;
	header1.npartTotal[5]=0;
	header1.flag_cooling=0;
	header1.num_files=1;
	header1.BoxSize=0.0;
	header1.Omega0=0.0;
	header1.OmegaLambda=0.0;
	header1.HubbleParam=0.0;
	fprint_header(outstr);
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr,RFMT RFMT RFMT "\n",Pos(p)[0],Pos(p)[1],Pos(p)[2]);
	}
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr,RFMT RFMT RFMT "\n",Vel(p)[0],Vel(p)[1],Vel(p)[2]);
	}
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr,IFMT "\n",p-bodytab+1);
	}
	for (p = bodytab; p < bodytab+nbody; p++) { 
		fprintf(outstr, RFMT "\n",Mass(p));
	}
	for (p = bodytab; p < bodytab+nbody; p++) { 
//		InternalEnergy = ((double) 1.5) * SOUNDSPEEDSQR;
//		InternalEnergy = ((double) 1.5) * rsqr(SoundSpeed);	// Standar model
		InternalEnergy = 0.05;	// Value to test (Evrard, Ver Articulo de Hernquist & Katz)
		fprintf(outstr, RFMT "\n",InternalEnergy);
	}
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
	fclose(outstr);
}

#undef IMFT
#undef RMFT


local void fprint_header(FILE *fd)
{
  fprintf(fd,"%d\n", header1.npart[0]);
  fprintf(fd,"%d\n", header1.npart[1]);
  fprintf(fd,"%d\n", header1.npart[2]);
  fprintf(fd,"%d\n", header1.npart[3]);
  fprintf(fd,"%d\n", header1.npart[4]);
  fprintf(fd,"%d\n", header1.npart[5]);
  fprintf(fd,"%e\n", header1.mass[0]);
  fprintf(fd,"%e\n", header1.mass[1]);
  fprintf(fd,"%e\n", header1.mass[2]);
  fprintf(fd,"%e\n", header1.mass[3]);
  fprintf(fd,"%e\n", header1.mass[4]);
  fprintf(fd,"%e\n", header1.mass[5]);
  fprintf(fd,"%e\n", header1.time);
  fprintf(fd,"%e\n", header1.redshift);
  fprintf(fd,"%d\n", header1.flag_sfr);
  fprintf(fd,"%d\n", header1.flag_feedback);
  fprintf(fd,"%d\n", header1.npartTotal[0]);
  fprintf(fd,"%d\n", header1.npartTotal[1]);
  fprintf(fd,"%d\n", header1.npartTotal[2]);
  fprintf(fd,"%d\n", header1.npartTotal[3]);
  fprintf(fd,"%d\n", header1.npartTotal[4]);
  fprintf(fd,"%d\n", header1.npartTotal[5]);
  fprintf(fd,"%d\n", header1.flag_cooling);
  fprintf(fd,"%d\n", header1.num_files);
  fprintf(fd,"%e\n", header1.BoxSize);
  fprintf(fd,"%e\n", header1.Omega0);
  fprintf(fd,"%e\n", header1.OmegaLambda);
  fprintf(fd,"%e\n", header1.HubbleParam);
}

local void fprint_header_reducido(FILE *fd)
{
  fprintf(fd,"%d\n", header_reducido.npart[0]);
  fprintf(fd,"%d\n", header_reducido.npart[1]);
  fprintf(fd,"%d\n", header_reducido.npart[2]);
  fprintf(fd,"%d\n", header_reducido.npart[3]);
  fprintf(fd,"%d\n", header_reducido.npart[4]);
  fprintf(fd,"%d\n", header_reducido.npart[5]);
  fprintf(fd,"%e\n", header_reducido.mass[0]);
  fprintf(fd,"%e\n", header_reducido.mass[1]);
  fprintf(fd,"%e\n", header_reducido.mass[2]);
  fprintf(fd,"%e\n", header_reducido.mass[3]);
  fprintf(fd,"%e\n", header_reducido.mass[4]);
  fprintf(fd,"%e\n", header_reducido.mass[5]);
  fprintf(fd,"%e\n", header_reducido.time);
  fprintf(fd,"%d\n", header_reducido.npartTotal[0]);
  fprintf(fd,"%d\n", header_reducido.npartTotal[1]);
  fprintf(fd,"%d\n", header_reducido.npartTotal[2]);
  fprintf(fd,"%d\n", header_reducido.npartTotal[3]);
  fprintf(fd,"%d\n", header_reducido.npartTotal[4]);
  fprintf(fd,"%d\n", header_reducido.npartTotal[5]);
  fprintf(fd,"%d\n", header_reducido.num_files);
  fprintf(fd,"%e\n", header_reducido.BoxSize);
}

// Print all the bodies like halo type... and BoxSize=0...
local void outputdata_gadget11_bin_double_reducido(char *file, int num, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
  FILE *fd;
  char buf[100];
  double dummy[3];
  int i,k;
  int   blklen,masscount;
  double a3inv;
  bodyptr p;

    char namebuf[256];

//#define BLKLEN gdgt_fwrite(&blklen, sizeof(blklen), 1, fd);

    a3inv=  1.0;

    sprintf(namebuf, file, num);
  sprintf(buf,"%s%s",All.OutputDir,namebuf);

	All.TotN_p0=All.TotN_disk=All.TotN_bulge=All.TotN_stars=All.TotN_dm=0;
	NumPart=All.TotN_halo=nbody;
	NumParticleTypes=6;
	for (i=0; i<NumParticleTypes; i++)
		All.MassTable[i]=0;
	All.BoxSize=0.;

fprintf(stdout,"\nTransfering bodytab info to Particle structure ...");
	for (p = bodytab; p < bodytab+nbody; p++) {
		i = p-bodytab+1;
		P[i].Type = 1;								// Halo particles...
		P[i].Mass=Mass(p);
		P[i].ID=Id(p);
		for (k=0; k<NDIM; k++) {
			P[i].PosPred[k]=Pos(p)[k];
			P[i].VelPred[k]=Vel(p)[k];
		}
	}
fprintf(stdout,"done ... \n%d bodies were transfered\n",i);

  if((fd=fopen(buf,"w")))
    {
      header_reducido.npart[0]= header_reducido.npartTotal[0]= All.TotN_p0;
      header_reducido.npart[1]= header_reducido.npartTotal[1]= All.TotN_halo;
      header_reducido.npart[2]= header_reducido.npartTotal[2]= All.TotN_disk;
      header_reducido.npart[3]= header_reducido.npartTotal[3]= All.TotN_bulge;
      header_reducido.npart[4]= header_reducido.npartTotal[4]= All.TotN_stars;
      header_reducido.npart[5]= header_reducido.npartTotal[5]= All.TotN_dm;

      for(i=0;i<NumParticleTypes;i++)
	header_reducido.mass[i]=0;

      for(i=0, masscount=0; i<NumParticleTypes; i++)
	{
	  header_reducido.mass[i]= All.MassTable[i];
	  if(All.MassTable[i]==0 && header_reducido.npart[i]>0)
	    masscount+= header_reducido.npart[i];
	}

      header_reducido.time= All.Time;

      header_reducido.num_files= 1;
      header_reducido.BoxSize= All.BoxSize;
      
      blklen=sizeof(header_reducido);
      BLKLEN;
      gdgt_fwrite(&header_reducido, sizeof(header_reducido), 1, fd);
      BLKLEN;

      blklen=NumPart*3*sizeof(double);

      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P[i].PosPred[k];
	  gdgt_fwrite(dummy,sizeof(double),3,fd);	
	}
      BLKLEN;

      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P[i].VelPred[k];
	  gdgt_fwrite(dummy,sizeof(double),3,fd);	
	}
      BLKLEN;
 
      blklen=NumPart*sizeof(int);
      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  gdgt_fwrite(&P[i].ID,sizeof(int),1,fd);
	}
      BLKLEN;

      blklen=masscount*sizeof(double);
      if(masscount)
	BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  dummy[0]= P[i].Mass;
	  if(All.MassTable[P[i].Type]==0)
	    gdgt_fwrite(dummy,sizeof(double),1,fd);		
	}
      if(masscount)
	BLKLEN;

      fclose(fd);
    }
  else
    {
      fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
      endrun(10);
    }
}


// Print all the bodies like halo type... and BoxSize=0...
local void outputdata_gadget11_ascii_reducido(char *file, int num, 
		int nbody, real tnow, io_header_blj *hdr, char *options, short allocate_mode)
{
  FILE *fd;
  char buf[100];
  double dummy[3];
  int i,k;
  int   blklen,masscount;
  double a3inv;
  vector dummyvec;
  double realdummy;
  bodyptr p;

    char namebuf[256];

    a3inv=  1.0;

    sprintf(namebuf, file, num);
  sprintf(buf,"%s%s",All.OutputDir,namebuf);

	All.TotN_p0=All.TotN_disk=All.TotN_bulge=All.TotN_stars=All.TotN_dm=0;
	NumPart=All.TotN_halo=nbody;
	NumParticleTypes=6;
	for (i=0; i<NumParticleTypes; i++)
		All.MassTable[i]=0;
	All.BoxSize=0.;

	All.TotNumPart = NumPart;
	All.TotN_gas=0;
	All.PartAllocFactor = 1.6;
	All.TreeAllocFactor =	0.8;
	All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;    
	All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;
	AllocateMemory(allocate_mode);

fprintf(stdout,"\nTransfering bodytab info to Particle structure ...");
	for (p = bodytab; p < bodytab+nbody; p++) {
		i = p-bodytab+1;
		P[i].Type = 1;								// Halo particles...
		P[i].Mass=Mass(p);
		P[i].ID=Id(p);
		for (k=0; k<NDIM; k++) {
			P[i].PosPred[k]=Pos(p)[k];
			P[i].VelPred[k]=Vel(p)[k];
		}
	}
fprintf(stdout,"done ... \n%d bodies were transfered\n",i);

  if((fd=fopen(buf,"w")))
    {
      header_reducido.npart[0]= header_reducido.npartTotal[0]= All.TotN_p0;
      header_reducido.npart[1]= header_reducido.npartTotal[1]= All.TotN_halo;
      header_reducido.npart[2]= header_reducido.npartTotal[2]= All.TotN_disk;
      header_reducido.npart[3]= header_reducido.npartTotal[3]= All.TotN_bulge;
      header_reducido.npart[4]= header_reducido.npartTotal[4]= All.TotN_stars;
      header_reducido.npart[5]= header_reducido.npartTotal[5]= All.TotN_dm;

      for(i=0;i<NumParticleTypes;i++)
	header_reducido.mass[i]=0;

      for(i=0, masscount=0; i<NumParticleTypes; i++)
	{
	  header1.mass[i]= All.MassTable[i];
	  if(All.MassTable[i]==0 && header_reducido.npart[i]>0)
	    masscount+= header_reducido.npart[i];
	}

      header_reducido.time= All.Time;


      
      header_reducido.num_files= 1;
      header_reducido.BoxSize= All.BoxSize;

	fprint_header_reducido(fd);

      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummyvec[k]=P[i].PosPred[k];
        out_vector(fd, dummyvec);
	}

      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummyvec[k]=P[i].VelPred[k];
        out_vector(fd, dummyvec);
	}

      for(i=1;i<=NumPart;i++)
	{
		out_int(fd,P[i].ID);
	}


		for(i=1;i<=NumPart;i++) {
			realdummy= P[i].Mass;
			if(All.MassTable[P[i].Type]==0)
			out_real(fd, realdummy);
		}

		if (scanopt(options, "out-phi"))
			for(i=1;i<=NumPart;i++) {
				realdummy= P[i].Potential;
				out_real(fd, realdummy);
			}
		if (scanopt(options, "out-acc"))
			for(i=1;i<=NumPart;i++) {
				for(k=0;k<3;k++)
					dummyvec[k]=P[i].Accel[k];
				out_vector(fd, dummyvec);
			}

      fclose(fd);
    }
  else
    {
      fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
      endrun(10);
    }
}

// -------------COMIENZAN OUTPUT ROUTINES FOR BLJ DATA FORMAT-------------------


local void outputdata_tlj(char *file, int snapcount, 
	int nbody, real tnow, io_header_tlj *hdr, char *options)
{
    char namebuf[256], buf[256];
    stream outstr;
    bodyptr p;

    outstr = stropen(SNAPTMP, "w!");

	output_header_tlj(outstr, *hdr);

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_int(outstr, Id(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_short(outstr, Type(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_real(outstr, Mass(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_vector(outstr, Pos(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_vector(outstr, Vel(p));
    if (scanopt(options, "out-phi"))
		DO_BODY(p, bodytab, bodytab+hdr->nbody)
            out_real(outstr, Phi(p));
    if (scanopt(options, "out-acc"))
		DO_BODY(p, bodytab, bodytab+hdr->nbody)
            out_vector(outstr, Acc(p));
    fclose(outstr);

    sprintf(namebuf, file, snapcount);           
	sprintf(buf,"mv %s %s",SNAPTMP,namebuf);
	system(buf);

    printf("\n\tdata output to file %s at time %f\n", namebuf, hdr->tnow);
}

local void outputdata_blj(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256], buf[256];
    stream outstr;
    bodyptr p;

    outstr = stropen(SNAPTMP, "w!");

	output_header_blj(outstr, *hdr);

// Id
// Type
// Mass
// Pos
// Vel
// Phi
// Acc

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_int(outstr, Id(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_short(outstr, Type(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_real(outstr, Mass(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_vector(outstr, Pos(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_vector(outstr, Vel(p));
    if (scanopt(options, "out-phi"))
		DO_BODY(p, bodytab, bodytab+hdr->nbody)
            out_real(outstr, Phi(p));
    if (scanopt(options, "out-acc"))
		DO_BODY(p, bodytab, bodytab+hdr->nbody)
            out_vector(outstr, Acc(p));
    fclose(outstr);

    sprintf(namebuf, file, snapcount);           
	sprintf(buf,"mv %s %s",SNAPTMP,namebuf);
	system(buf);

    printf("\n\tdata output to file %s at time %f\n", namebuf, hdr->tnow);
}

local void outputbindata_blj(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, file, snapcount);           
    outstr = stropen(namebuf, "w!");

	output_header_blj_bin(outstr, *hdr);

	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_int_bin(outstr, Id(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_short_bin(outstr, Type(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_real_bin(outstr, Mass(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_vector_bin(outstr, Pos(p));
	DO_BODY(p, bodytab, bodytab+hdr->nbody)
        out_vector_bin(outstr, Vel(p));
    if (scanopt(options, "out-phi"))
		DO_BODY(p, bodytab, bodytab+hdr->nbody)
            out_real_bin(outstr, Phi(p));
    if (scanopt(options, "out-acc"))
		DO_BODY(p, bodytab, bodytab+hdr->nbody)
            out_vector_bin(outstr, Acc(p));
	fclose(outstr);
    printf("\n\tdata output to file %s at time %f\n", namebuf, hdr->tnow);
}

local void outputpvdata_blj(char *file, int snapcount, 
	int nbody, real tnow, io_header_blj *hdr, char *options)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, file, snapcount);           
    outstr = stropen(namebuf, "w!");

	output_header_blj_pv(outstr, *hdr);

	DO_BODY(p, bodytab, bodytab+hdr->nbody) {
		out_int_mar(outstr, Id(p));
		out_short_mar(outstr, Type(p));
		out_real_mar(outstr, Mass(p));
        out_vector_mar(outstr, Pos(p));
        out_vector_mar(outstr, Vel(p));
        if (scanopt(options, "out-phi"))
            out_real_mar(outstr, Phi(p));
        if (scanopt(options, "out-acc"))
            out_vector_mar(outstr, Acc(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);
    printf("\n\tdata output to file %s at time %f\n", namebuf, hdr->tnow);
}

local void output_header_tlj(stream outstr, io_header_tlj hdr)
{
	hdr.ndim=NDIM;

    out_int(outstr, hdr.nbody);
    out_int(outstr, hdr.nbody1);
    out_int(outstr, hdr.nbody2);
    out_int(outstr, hdr.nbody3);
    out_int(outstr, hdr.ndim);
    out_real(outstr, hdr.tnow);
    out_real(outstr, hdr.temperature);
    out_real(outstr, hdr.density);
    out_real(outstr, hdr.mass1);
    out_real(outstr, hdr.mass2);
    out_real(outstr, hdr.mass3);
    out_real(outstr, hdr.Lx);
    out_real(outstr, hdr.Ly);
#ifdef THREEDIM
    out_real(outstr, hdr.Lz);
#endif
    out_real(outstr, hdr.eps11);
    out_real(outstr, hdr.eps12);
    out_real(outstr, hdr.eps13);
    out_real(outstr, hdr.eps22);
    out_real(outstr, hdr.eps23);
    out_real(outstr, hdr.eps33);
    out_real(outstr, hdr.sigma11);
    out_real(outstr, hdr.sigma12);
    out_real(outstr, hdr.sigma13);
    out_real(outstr, hdr.sigma22);
    out_real(outstr, hdr.sigma23);
    out_real(outstr, hdr.sigma33);
    out_real(outstr, hdr.Rcut11);
    out_real(outstr, hdr.Rcut12);
    out_real(outstr, hdr.Rcut13);
    out_real(outstr, hdr.Rcut22);
    out_real(outstr, hdr.Rcut23);
    out_real(outstr, hdr.Rcut33);
}

local void output_header_blj(stream outstr, io_header_blj hdr)
{
	hdr.ndim=NDIM;
// Header number of lines: 21 (3D); 20 (2D)

    out_int(outstr, hdr.nbody);         // 1
    out_int(outstr, hdr.nbody1);        // 2
    out_int(outstr, hdr.nbody2);        // 3
    out_int(outstr, hdr.ndim);          // 4
    out_real(outstr, hdr.tnow);         // 5
    out_real(outstr, hdr.temperature);  // 6
    out_real(outstr, hdr.density);      // 7
    out_real(outstr, hdr.mass1);        // 8
    out_real(outstr, hdr.mass2);        // 9
    out_real(outstr, hdr.Lx);           // 10
    out_real(outstr, hdr.Ly);           // 11
#ifdef THREEDIM
    out_real(outstr, hdr.Lz);           // 12
#endif
    out_real(outstr, hdr.eps11);        // 13
    out_real(outstr, hdr.eps12);        // 14
    out_real(outstr, hdr.eps22);        // 15
    out_real(outstr, hdr.sigma11);      // 16
    out_real(outstr, hdr.sigma12);      // 17
    out_real(outstr, hdr.sigma22);      // 18
    out_real(outstr, hdr.Rcut11);       // 19
    out_real(outstr, hdr.Rcut12);       // 20
    out_real(outstr, hdr.Rcut22);       // 21
}

local void output_header_blj_pv(stream outstr, io_header_blj hdr)
{
	hdr.ndim=NDIM;

	fprintf(outstr,
		"# nbody nbody1 nbody2 NDIM time temperature density masses Ls epss sigmas Rcuts\n");
	fprintf(outstr,"# ");
    out_int(outstr, hdr.nbody);
	fprintf(outstr,"# ");
    out_int(outstr, hdr.nbody1);
	fprintf(outstr,"# ");
    out_int(outstr, hdr.nbody2);
	fprintf(outstr,"# ");
    out_int(outstr, hdr.ndim);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.tnow);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.temperature);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.density);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.mass1);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.mass2);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.Lx);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.Ly);
#ifdef THREEDIM
	fprintf(outstr,"# ");
    out_real(outstr, hdr.Lz);
#endif
	fprintf(outstr,"# ");
    out_real(outstr, hdr.eps11);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.eps12);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.eps22);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.sigma11);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.sigma12);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.sigma22);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.Rcut11);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.Rcut12);
	fprintf(outstr,"# ");
    out_real(outstr, hdr.Rcut22);
}

local void output_header_blj_bin(stream outstr, io_header_blj hdr)
{
	hdr.ndim=NDIM;

    out_int_bin(outstr, hdr.nbody);
    out_int_bin(outstr, hdr.nbody1);
    out_int_bin(outstr, hdr.nbody2);
    out_int_bin(outstr, hdr.ndim);
    out_real_bin(outstr, hdr.tnow);
    out_real_bin(outstr, hdr.temperature);
    out_real_bin(outstr, hdr.density);
    out_real_bin(outstr, hdr.mass1);
    out_real_bin(outstr, hdr.mass2);
    out_real_bin(outstr, hdr.Lx);
    out_real_bin(outstr, hdr.Ly);
#ifdef THREEDIM
    out_real_bin(outstr, hdr.Lz);
#endif
    out_real_bin(outstr, hdr.eps11);
    out_real_bin(outstr, hdr.eps12);
    out_real_bin(outstr, hdr.eps22);
    out_real_bin(outstr, hdr.sigma11);
    out_real_bin(outstr, hdr.sigma12);
    out_real_bin(outstr, hdr.sigma22);
    out_real_bin(outstr, hdr.Rcut11);
    out_real_bin(outstr, hdr.Rcut12);
    out_real_bin(outstr, hdr.Rcut22);
}

#undef SNAPTMP

// --------------TERMINAN OUTPUT ROUTINES FOR BLJ DATA FORMAT-------------------


// --------------COMIENZAN ROUTINES FOR HEADER CONVERSION DATA FORMAT-----------
local void HeaderConversion_blj(int infmt, int outfmt, 
								void *hdrin, io_header_blj *hdrout)
{
	io_header_tlj *hdrtljin;

	switch(infmt) {
		case IO_SNAP_TLJ_FMT:
			hdrtljin = (io_header_tlj *) hdrin;
			hdrout->nbody		= hdrtljin->nbody;
			hdrout->nbody1		= hdrtljin->nbody1;
			hdrout->nbody2		= hdrtljin->nbody2;
			hdrout->tnow		= hdrtljin->tnow;
			hdrout->temperature	= hdrtljin->temperature;
			hdrout->density		= hdrtljin->density;
			hdrout->mass1		= hdrtljin->mass1;
			hdrout->mass2		= hdrtljin->mass2;
			hdrout->Lx			= hdrtljin->Lx;
			hdrout->Ly			= hdrtljin->Ly;
#ifdef THREEDIM
			hdrout->Lz			= hdrtljin->Lz;
#endif
			hdrout->eps11		= hdrtljin->eps11;
			hdrout->eps12		= hdrtljin->eps12;
			hdrout->eps22		= hdrtljin->eps22;
			hdrout->sigma11		= hdrtljin->sigma11;
			hdrout->sigma12		= hdrtljin->sigma12;
			hdrout->sigma22		= hdrtljin->sigma22;
			hdrout->Rcut11		= hdrtljin->Rcut11;
			hdrout->Rcut12		= hdrtljin->Rcut12;
			hdrout->Rcut22		= hdrtljin->Rcut22;
			break;
		default:
			error("\n\tinput: Unknown input format...");
			break;
	}
}

local void HeaderConversion_tlj(int infmt, int outfmt, 
								void *hdrin, io_header_tlj *hdrout)
{
	io_header_blj *hdrbljin;

	switch(infmt) {
		case IO_SNAP_BLJ_FMT:
			hdrbljin = (io_header_blj *) hdrin;
			hdrout->nbody		= hdrbljin->nbody;
			hdrout->nbody1		= hdrbljin->nbody1;
			hdrout->nbody2		= hdrbljin->nbody2;
			hdrout->tnow		= hdrbljin->tnow;
			hdrout->temperature	= hdrbljin->temperature;
			hdrout->density		= hdrbljin->density;
			hdrout->mass1		= hdrbljin->mass1;
			hdrout->mass2		= hdrbljin->mass2;
			hdrout->Lx			= hdrbljin->Lx;
			hdrout->Ly			= hdrbljin->Ly;
#ifdef THREEDIM
			hdrout->Lz			= hdrbljin->Lz;
#endif
			hdrout->eps11		= hdrbljin->eps11;
			hdrout->eps12		= hdrbljin->eps12;
			hdrout->eps22		= hdrbljin->eps22;
			hdrout->sigma11		= hdrbljin->sigma11;
			hdrout->sigma12		= hdrbljin->sigma12;
			hdrout->sigma22		= hdrbljin->sigma22;
			hdrout->Rcut11		= hdrbljin->Rcut11;
			hdrout->Rcut12		= hdrbljin->Rcut12;
			hdrout->Rcut22		= hdrbljin->Rcut22;
			break;
		default:
			error("\n\tinput: Unknown input format...");
			break;
	}
}
// --------------TERMINAN ROUTINES FOR HEADER CONVERSION DATA FORMAT-----------

//----------------------------- TERMINA BLOQUE I/O -----------------------------

