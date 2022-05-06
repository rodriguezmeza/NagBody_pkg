/* ==============================================================================
!	HEADER: global_defs.h														!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date: February 2005												!
!	Purpose: Definitions of global variables and parameters						!
!	Language: C																	!
!	Use: '#include "...."														!
!	Use in routines and functions: md_lj_tree (main), start_run, time_step,		!
!					tree_ljforcecalc, md_lj_tree_io								!
!	External headers: data_struct_defs											!
!	Comments and notes:															!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.							!
!		e-mail: marioalberto.rodriguez@inin.gob.mx
!		http://www.astro.inin.mx/mar											!
!																				!
!	Major revisions:															!
!	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.							!
!===============================================================================
!	Legal matters:																!
!	The author does not warrant that the program and routines it contains		!
!	listed below are free from error or suitable for particular applications,	!
!	and he disclaims all liability from any consequences arising from their		!
!	use.																		!
!==============================================================================*/

#ifndef _global_defs_h
#define _global_defs_h

#include "stdinc.h"
#include "vectdefs.h"

#include "data_struc_defs.h"

#define DIAGNOSTICS				// Definition to unbug the program
#undef DIAGNOSTICS

// Block of command line definitions

global string paramfile;
global string forcecalc_method;

global real theta;
global bool usequad;

global real temperature;
global real density;

global int stepEquil;

global int stepAvgVel;
global int stepVel;
global int sizeHistVel;
global real rangeVel;

global int stepAvgRdf;
global int stepRdf;
global int sizeHistRdf;
global real rangeRdf;

global int nbody;             
global real dtime;
global real tstop;

global string icfile;
global string icfilefmt;
global string snapoutfile;                 
global string snapoutfilefmt;
global real dtout;
global real dtoutinfo;
global string statefile;
global string restorefile;

global string options;                  
global int seed;

//------------------------------------------------------------------------------

global int forcemethod_int;
global int icfilefmt_int;

global int nbbcalc;                     

global real cpuforce;
global real cputotforce;
global real cpuinit;
global real cputotout;
global real cputotal;

global string headline0;
global string headline1;
global string headline2;
global string headline3;

global string model_comment;			

global real tnow;
global real tout;
global real toutinfo;
global int nstep;

global bodyptr bodytab;

global FILE *outlog;

global int stopflag;

// Specific to tree construction

global cellptr root;                    
global real rsize;                      
global int ncell;                       
global int tdepth;                      
global real cputree;                    

global int actmax;                      
global int nbccalc;                     


// Specific MD parameters and variables

global real eps;			// LJ Effective energy
global real sigma;			// LJ Effective length
global real mass;			// Mass of each molecule
global real Rcut;			// Cutoff radius

global real virSum;
global real Scale;
global real RcutSq;
global vector Box;
global real kinEnergy;
global real potEnergy;
global real totEnergy;
global real pressure;
global real vSum;
global real vvSum;
global real uSum;

global real sKinEnergy;
global real ssKinEnergy;
global real sTotEnergy;
global real ssTotEnergy;
global real sPressure;
global real ssPressure;

global int countVel;
global realptr histVel;
global int countRdf;

global realptr histRdf;

global int stepSnapInit;
global int stepSnap;
global int stepAvg;

global real kB;				// Boltzmann constant (=1.0)
global real vMag;			// to re-scale vels

global real lAvg;
global real dtcrit;
global real fAvg;

global real fphi;			// Needed in ljforcecalc
global real ssq;			// Needed in ljforcecalc
global real fa;				// Needed in ljforcecalc

global vectorI cells;
global int *cellList;

global char mode[2];		// thermo file mode

#define DO_CELL(j, m)  for (j = cellList[m]; j >= 0; j = cellList[j])

/*
#define VWrap(v, t)											\
	v[t] = v[t] - ( (real) (nint(v[t]/Box[t])) )*Box[t]
*/
// Alternative definition for VWrap

#define VWrap(v, t)											\
   if (v[t] >= 0.5 * Box[t])      v[t] -= Box[t];			\
   else if (v[t] < -0.5 * Box[t]) v[t] += Box[t]


#define VCellWrap(t)                                        \
   if (m2v[t] >= cells[t]) {								\
     m2v[t] = 0;											\
     shift[t] = Box[t];										\
   } else if (m2v[t] < 0) {									\
     m2v[t] = cells[t] - 1;                                 \
     shift[t] = - Box[t];									\
   }

#if NDIM == 2
#define VWrapAll(v)                                         \
   {VWrap (v, 0);                                           \
   VWrap (v, 1);}
#define VCellWrapAll()                                      \
   {VCellWrap (0);                                          \
   VCellWrap (1);}
#define OFFSET_VALS                                         \
   {{0,0}, {1,0}, {1,1}, {0,1}, {-1,1}}
#define N_OFFSET  5
#endif

#if NDIM == 3
#define VWrapAll(v)                                         \
   {VWrap (v, 0);                                           \
   VWrap (v, 1);                                            \
   VWrap (v, 2);}
#define VCellWrapAll()                                      \
   {VCellWrap (0);                                          \
   VCellWrap (1);                                           \
   VCellWrap (2);}
#define OFFSET_VALS                                         \
   {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0},           \
	{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, \
	{-1,-1,1}, {0,-1,1}, {1,-1,1}                           \
   }

#define N_OFFSET  14
#endif

// Definitions to add a parameter in the scheme of parameterfile

#define IPName(param)										\
  {strcpy(tag[nt],"param");									\
  addr[nt]=&(param);										\
  id[nt++]=INT;}

#define RPName(param)										\
  {strcpy(tag[nt],"param");									\
  addr[nt]=&param;											\
  id[nt++]=DOUBLE;}

#define BPName(param)										\
  {strcpy(tag[nt],"param");									\
  addr[nt]=&param;											\
  id[nt++]=BOOLEAN;}

#define SPName(param,n)										\
  {strcpy(tag[nt],"param");									\
  param=(string) malloc(n);									\
  addr[nt]=param;											\
  id[nt++]=STRING;}


// STATIC problem: gcc version 11
// From inout.h
//global real *inout_xval;
//global real *inout_yval;
//global real *inout_zval;
//global real *inout_wval;

// STATIC problem: gcc version 11
// From diffeqs.h
//global double dxsav,*xp,**yp;
//global int kmax,kount;
//global int nrhs;

// STATIC problem: gcc version 11
// From stdinc.h
global long idum;                // seed for random generators


#endif /* ! _global_defs_h */

