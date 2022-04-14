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

#ifndef NOGNU
#include "./general_libs/stdinc.h"
#include "./general_libs/mathfns.h"
#include "./general_libs/vectdefs.h"
#include "./general_libs/vectmath.h"
#include "./general_libs/getparam.h"
#include "./general_libs/machines.h"
#include "./general_libs/inout.h"
// #include <strings.h>							// For unix
#include "./general_libs/strings.h"	// For Visual c
#else
#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "vectmath.h"
#include "getparam.h"
#include "machines.h"
#include "inout.h"
// #include <strings.h>							// For unix
#include "strings.h"	// For Visual c
#endif

#include "data_struc_defs.h"
//===============================================
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

#endif // ! _globaldefs_h

