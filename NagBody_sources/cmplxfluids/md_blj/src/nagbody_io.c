/*==============================================================================
	MODULE: nagbody_io.c			[NagBody]
	Written by: Mario A. Rodriguez-Meza
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
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"
#include "nagbody.h"

#include <string.h>								// For Unix
//#include "../../../../General_libs/strings.h"	// For Visual C

#include <sys/stat.h>


//----------------------------- COMIENZA BLOQUE I/O -----------------------------

local void infilefmt_string_to_int(string, int *);

local void outfilefmt_string_to_int(string,int *);


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

	switch(infmt_int) {
		case IO_NULL_FMT:
			inputdata_blj_ascii(file, step, nbodies, ndim, tnow, exist_snap,
                                hdr, options);
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
/*		case IO_SNAP_TLJ_FMT:
			inputdata_tlj_ascii(file, step, nbodies, ndim, tnow, exist_snap,
				hdr, options);
			break; */
		default:
			printf("\n\tinput: Unknown input format...");
			printf("\n\tinput in default snap (ascii) format...\n"); 
			inputdata_blj_ascii(file, step, nbodies, ndim, tnow, exist_snap,
                                hdr, options);
			break;
	}
}

local void infilefmt_string_to_int(string infmt_str,int *infmt_int)
{
    *infmt_int=-1;
    if (strnull(infmt_str))
		*infmt_int = IO_NULL_FMT;
    if (strcmp(infmt_str,"snap-blj-ascii") == 0)
		*infmt_int = IO_SNAP_BLJ_FMT;
    if (strcmp(infmt_str,"snap-blj-pv") == 0)
		*infmt_int = IO_SNAP_BLJ_PV_FMT;
    if (strcmp(infmt_str,"snap-blj-bin") == 0)
		*infmt_int = IO_SNAP_BLJ_BIN_FMT;
/*    if (strcmp(infmt_str,"snap-tlj-ascii") == 0)
		*infmt_int = IO_SNAP_TLJ_FMT; */
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
            case IO_NULL_FMT: 
                printf("\n\tsnap-blj-ascii format output");
				if (infmt_int == outfilefmt_int)
					outputdata_blj(file, snapcount, nbody, tnow, hdr, options);
				else {
					hdrout = (io_header_blj *) allocate(sizeof(io_header_blj));
					HeaderConversion_blj(infmt_int,outfilefmt_int,hdr,hdrout);
					outputdata_blj(file,snapcount,nbody,tnow,hdrout,options);
				}
				break;
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
/*
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
*/
            default:
                printf("\n\toutput: Unknown output format...");
                printf("\n\tprinting in default snap-ascii format..."); 
				outputdata_blj(file, snapcount, nbody, tnow, hdr, options); break;
        }
    } else
		error("\n\noutputdata : You should give an output file name\n");
}


local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strnull(outfmt_str))
		*outfmt_int = IO_NULL_FMT;
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

