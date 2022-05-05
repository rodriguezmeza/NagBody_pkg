/*==============================================================================
	HEADER: globaldefs.h			[galaxy_hernquist]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: main, direct_gravcalc,
					nbody_n2_io, startrun, timestep
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 2007; November 2008;
	Copyright: (c) 2005-2014 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h

// =================================================================
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>


//#ifdef HAVE_HDF5
//#include <hdf5.h>
//#endif


#ifndef NOGNU
#include "./general_libs/general/stdinc.h"
#include "./general_libs/general/getparam.h"
#else
#include "stdinc.h"
#include "getparam.h"
#endif

// =================================================================

#include "data_struc_defs.h"

typedef struct {
	string options;                  

    int seed;

    string paramfile;
	string snapoutfile;
	string snapoutfilefmt;

    int Ngas;
    int Nhalo;
    int Ndisk;
    int Nbulge;

    real Mgas;
    real Mhalo;
    real Mdisk;
    real Mbulge;

    real masscut;

    real ag;
    real ab;
    real ad;
    real ah;
    
    real zg;
    real zd;
    real gammah;
    real gammab;
    
    real rmaxg;
    real rmaxb;
    real rmaxd;
    real rmaxh;

} cmdline_data, *cmdline_data_ptr;

typedef struct {
	real cpuinit;

	char model_comment[100];

	string headline0;
	string headline1;
	string headline2;
	string headline3;

	FILE *outlog;

	char mode[4];

} global_data, *global_data_ptr;

global global_data gd;
global cmdline_data cmd;

extern double vcirs[];

// STATIC problem: gcc version 11
// From inout.h
global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_wval;

// STATIC problem: gcc version 11
// From diffeqs.h
global double dxsav,*xp,**yp;
global int kmax,kount;
global int nrhs;

// STATIC problem: gcc version 11
// From stdinc.h
global long idum;                // seed for random generators


#endif // ! _globaldefs_h

