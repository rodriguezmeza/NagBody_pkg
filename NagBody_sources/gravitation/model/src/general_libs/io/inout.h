/*==============================================================================
	HEADER: inout.h				[General_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date:	January, 2005
	Purpose: Headers of utilities for input and output data
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

	Major revisions:
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#ifndef _inout_h
#define _inout_h

#include "vectdefs.h"

// ------------[	inout normal definitions	 	]------------

void in_int(stream, int *);
void in_short(stream, short *);
void in_real(stream, real *);
void in_vector(stream, vector);
void out_int(stream, int);
void out_short(stream, short);
void out_real(stream, real);
void out_vector(stream, vector);


// ------------[	inout mar definitions	 		]------------

void out_int_mar(stream, int);
void out_short_mar(stream, short);
void out_bool_mar(stream, bool);
void out_real_mar(stream, real);          
void out_vector_mar(stream, vector);      


// ------------[	inout binary definitions	 	]------------

void in_int_bin(stream, int *);               
void in_short_bin(stream, short *);
void in_real_bin(stream, real *);             
void in_vector_bin(stream, vector);           
void out_int_bin(stream, int);
void out_short_bin(stream, short);
void out_real_bin(stream, real);
void out_vector_bin(stream, vector);
void out_bool_mar_bin(stream, bool);


// ------------[	inout other definitions	 		]------------

void in_vector_ruben(stream, vector);           

void out_vector_ndim(stream, double *, int);
void in_vector_ndim(stream, real *, int);

// Macros for binary in/out
#define safewrite(ptr,len,str)                  \
    if (fwrite((void *) ptr, len, 1, str) != 1) \
        error("safewrite: fwrite failed\n")

#define saferead(ptr,len,str)                  \
    if (fread((void *) ptr, len, 1, str) != 1) \
        error("saferead: fread failed\n")


// Macros and routines for binary gdgt in/out:: Uncomment if necessary
//size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
//size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
// This are the same as above. Just to be explicit their use in gadget type of codes.
size_t gdgt_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t gdgt_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

#define SKIP gdgt_fread(&blklen,sizeof(int4byte),1,fd);
#define BLKLEN gdgt_fwrite(&blklen, sizeof(blklen), 1, fd);

// IN/OUT ROUTINES TO HANDLE STRINGS ...
void ReadInString(stream, char *);
void ReadInLineString(stream, char *);


// InputData (like the one used by nplot2d)
void inout_InputData(string, int, int, int *);

#endif	// ! _inout_h
