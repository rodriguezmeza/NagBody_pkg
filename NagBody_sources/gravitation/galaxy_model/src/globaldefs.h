/*==============================================================================
	HEADER: globaldefs.h		[galaxy_models]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "global_defs.h"
	Use in routines and functions:
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: June 6, 2008;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved.
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h


//===============================================
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include "stdinc.h"
#include "vectdefs.h"
#include "vectmath.h"
#include "nagbody.h"
#include "constant.h"
#include "mathfns.h"
#include "numrec.h"
#include "inout.h"
#include "getparam.h"

//===============================================

#define nbodsmax	1500000
#define maxtabh		30000
#define maxtab		210				// maximum number of radial bins in tables


global realptr surfd;				// disk gas composition. array [1:nbodsmax]

// Particle composition ==== All arrays [1:nbodsmax]
global realptr x,
			y,
			z,
			vx,
			vy,
			vz,
            pmass,
			radcyl,
			radsph,
			ax,
			ay,
			az,
			aradcyl,
			kappa,
			sigr,						// radial velocity dispersion
			sigphi,						// phi velocity dispersion
			sigz,						// z velocity dispersion
			sigt,						// critical Toomre dispersion (xQ)
			pot,
			rotcirc,					// circular equilibrium velocity
            rotmean;					// mean circular velocity

// Halo composition ==== Arrays [1:maxtabh]
global realptr rhalo,
			xmhalo,
			uhalo;

global realptr dadrtab;					// radial acceleration gradient table. 
										// Array [1:maxtab]


typedef struct {
	string paramfile;

	string filename;
	string filenamefmt;
	string statfile;
	long seed;

	bool outpteps;

// Disk's parameters ======
	bool	usedisk;
	bool	usegas,
			selfggas;
	int		ndstars,
			ndgas;
	real	z0,								// disk scale height
			qsolar,							// Toomre Q parameter and radius 
											// in the solar neighborhood
			epsdisk,
			zmax,
			rmax;

	string	rsolarstr;						// Radius in the solar neighborhood

// Disk gas's parameters ==
	real	gasmass,
            gastemp,
			z0gas,
			zmaxgas,
			rmaxgas,
			rmingas;

// Bulge's parameters =====
	bool	usebulge,
			selfgbul,
			axibulge,
			bulgerot;
	int		nbulge,
			nsimpson;
	real	bulgmass,
			abulge,
			epsbulge,
			rmaxbulg,
			zmaxbulg,
			brotfrac,
			cbulge;

// Halo's parameters ======
	bool	usehalo,
			selfghal;
	int		nhalo;
	real	halomass,
			ahalo,
			gamhalo,
			rthalo,
			rmaxhalo,
			epshalo;
	string	halotype;

// Satellite's parameters =
	bool	usesat,
			selfgsat;
	int		nsat;
	real	satmass,
			asat,
			xsat,
			ysat,
			zsat,
            vxsat,
			vysat,
			vzsat,
			rmaxsat,
            epssat;

// Add models composition's parameters
	bool	addmods;
	real	thetmod1,
			phimod1,
			thetmod2,
			phimod2,
			rp,
			rsep;

} cmdline_data, *cmdline_data_ptr;


typedef struct {

	real cpuinit;
	string comment;			
	FILE *outlog;

	string headline0;
	string headline1;
	string headline2;
	string headline3;

	real	h,								// Units system
			diskmass,
			G;

	int		nbodies;						// total number of particles	

// Disk's parameters ======
	int		ndisk;
	real	epsdisk2,
			rsolar;							// Radius in the solar neighborhood

// Disk gas's parameters ==
	real	xgasmass,
			surfd0,
			sigr0,
			acorr,
			acorrgas;

// Halo's parameters ======
	int		ntabhalo;

// Satellite's parameters =
	real	axsat,
			aysat,
			azsat,
			potsat,
			radsat;

// Add models composition's parameters
	real	xmod1,
			ymod1,
            zmod1,
			vxmod1,
			vymod1,
			vzmod1,
			xmod2,
			ymod2,
			zmod2,
			vxmod2,
            vymod2,
			vzmod2;

} global_data, *global_data_ptr;

global global_data gd;
global cmdline_data cmd;
global global_data_treegrav gdtreegrav;

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
static long idum;                // seed for random generators

#endif /* ! _globaldefs_h */

