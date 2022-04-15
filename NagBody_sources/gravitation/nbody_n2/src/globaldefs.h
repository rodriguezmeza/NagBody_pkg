/*==============================================================================
	HEADER: globaldefs.h			[nbody_n2]
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
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.																		!
==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h

//===============================================
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "vectmath.h"
#include "getparam.h"
#include "machines.h"
#include "inout.h"
#include "strings.h"

#include "data_struc_defs.h"
#include "protodefs.h"

typedef struct {
	string options;                  

	real eps;

	int seed;

	string icfile;
	string icfilefmt;
	string paramfile;
	string snapoutfile;
	string snapoutfilefmt;
	string statefile;
	int stepState;
	string restorefile;

	string dtimestr;
	string dtoutstr;
	string dtoutinfostr;
	real tstop;

	int nbody;             

} cmdline_data, *cmdline_data_ptr;

typedef struct {
	int nbbcalc;                     

	real cpuforce;
	real cpuinit;

	real dtime;
	real dtout;
	real dtoutinfo;

	char model_comment[100];

	string headline0;
	string headline1;
	string headline2;
	string headline3;

	real tnow;
	real tout;
	real toutinfo;
	int nstep;

	FILE *outlog;

	int stopflag;

	char mode[2];

} global_data, *global_data_ptr;

global global_data gd;
global cmdline_data cmd;
global bodyptr bodytab;      

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

