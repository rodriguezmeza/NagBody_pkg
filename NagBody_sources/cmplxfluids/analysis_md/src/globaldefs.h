/* ==============================================================================
	HEADER: globaldefs.h			[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: analysis_md (main), start_run, time_step,
					forcecalc, md_lj_tree_io
	External headers: data_struct_defs
	Comments and notes:
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
===============================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _global_defs_h
#define _global_defs_h

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

#ifndef NOGNU
#include "./general_libs/stdinc.h"
#include "./general_libs/vectdefs.h"
#include "./general_libs/nagbody.h"
#include "./general_libs/physconstants.h"
#include "./general_libs/mathfns.h"
#include "./general_libs/vectmath.h"
#include "./general_libs/getparam.h"
#include "./general_libs/inout.h"
#include "./general_libs/constant.h"
/* #include <strings.h> */							// For unix
#include "./general_libs/strings.h"	// For Visual C
#include "./general_libs/pldefs.h"
#include "./general_libs/lic.h"
#include "./general_libs/machines.h"
#else
#include "stdinc.h"
#include "vectdefs.h"
#include "nagbody.h"
#include "physconstants.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "inout.h"
#include "constant.h"
/* #include <strings.h> */							// For unix
#include "strings.h"	// For Visual C
#include "pldefs.h"
#include "lic.h"
#include "machines.h"
#endif

#include "protodefs.h"


typedef struct {
	string paramfile;

	real temperature;
	real density;

//	string nbodyprop;
//	string massprop;
//#ifdef THREEDIM
//	string LxLyprop;
//#else
//	string Lx;
//#endif

	real eps11;
	real eps12;
	real eps22;
	real sigma11;
	real sigma12;
	real sigma22;
	real Rcut11;
	real Rcut12;
	real Rcut22;

//	int stepEquil;

	int stepAvg;
	int sizeHist;
	real rangeVal;

//	int stepAvgRhoAxes;
//	int sizeHistRhoAxes;

//	int stepAvgNFrecAxes;
//	int sizeHistNFrecAxes;

//	int stepAvgVel;
//	int sizeHistVel;
//	real rangeVel;

//	int stepAvgRdf;
//	int sizeHistRdf;
//	real rangeRdf;

//	int stepAvgPressAxes;
//	int sizeHistPressAxes;

//	int stepAvgChemPot;
//	int sizeHistChemPot;
//	int numTestBodies;

//	int limitAcfAv;
//	int nBuffAcf;
//	int nValAcf;

//	int nbody;
//	real dtime;
//	real dtout;                

	string in;
	string infmt;
	string out;                 
	string outfmt;
	string basedir;

	int isnap;
	int fsnap;

	string unitsset;

	string options;

	string data_analysis_type;

	real reductionFac;

	string bodiesID;

	string xrange;
	string yrange;
#ifdef THREEDIM
	string zrange;

	real xmin;
	real xmax;
	real ymin;
	real ymax;
	real zmin;
	real zmax;
#endif

	string usingcolumns;
	string usingrows;

	string xlabel;
	string ylabel;
	real labelfontsize;
	string plotlabel;
	bool plotjoined;
	bool withdots;
	bool withsymbols;
	real symbolsize;

	int labelfontweight;
	real nlsize;
	int linewidth;
	int axeswidth;
	int symbolweight;
	int symbolcolor;

	string pl_a;
	string pl_dev;
	string pl_geo;
	string pl_ori;
	string pl_bg;
	string pl_ncol0;
	string pl_ncol1;
} cmdline_data, *cmdline_data_ptr;

#define MAXBODIESID		20

typedef struct {
	int nbody;
	int nbody1;
	int nbody2;
	real mass1;
	real mass2;
	real Lx;
	real Ly;
#ifdef THREEDIM
	real Lz;
#endif
	string headline0;
	string headline1;
	string headline2;
	string headline3;
//	string model_comment;
	char model_comment[100];

	string headerfmt;

	real tnow;
	real tout;
	int nstep;
	real cpuinit;
	FILE *outlog;
	real mass;
	vector Box;

	int count;

//	int countRhoAxes;
	realptr histRhoX;
	realptr histRhoY;

	realptr histRhoX1;
	realptr histRhoY1;

	realptr histRhoX2;
	realptr histRhoY2;

#ifdef THREEDIM
	realptr histRhoZ;
	realptr histRhoZ1;
	realptr histRhoZ2;
#endif

//	int countNFrecAxes;
	int *histNFrecX;
	int *histNFrecY;

	int *histNFrecX1;
	int *histNFrecY1;

	int *histNFrecX2;
	int *histNFrecY2;

	int *histNFrecXD;
	int *histNFrecYD;

#ifdef THREEDIM
	int *histNFrecZ;
	int *histNFrecZ1;
	int *histNFrecZ2;
	int *histNFrecZD;
#endif

//	int countVel;
	realptr histVel;

//	int countRdf;
	realptr histRdf;

//	int countPressAxes;
//	realptr histPressX;
//	realptr histPressY;
//#ifdef THREEDIM
//	realptr histPressZ;
//#endif

//	int countChemPot;
//	realptr histChemPotX;
//	realptr histChemPotY;
//#ifdef THREEDIM
//	realptr histChemPotZ;
//#endif

//	real ChemPot;
//	real sChemPot;
//	real ssChemPot;
	
//	int countAcfAv;
//	realptr avAcfVel;
//	real intAcfVel;
//	realptr avAcfTherm;
//	realptr avAcfVisc;
//	real intAcfTherm;
//	real intAcfVisc;

//	string snapfilename;

	real unitLength;
	real unitEnergy;
	real unitMass;

	int nbodiesID;
	int bodyID[MAXBODIESID];

	bool x_autoscale;
	bool y_autoscale;
#ifdef THREEDIM
	bool z_autoscale;
#endif
	real xmin;
	real xmax;
	real ymin;
	real ymax;
#ifdef THREEDIM
	real zmin;
	real zmax;
#endif

	int column1;
	int column2;
	int column3;
	int row1;
	int row2;

	int Rdf_flag;
	int Vel_flag;
	int RhoAxes_flag;
	int NFrecAxes_flag;
} global_data, *global_data_ptr;

#undef MAXBODIESID

typedef struct {
  vector *orgVel;
  real *acfVel;
  int count;
  vector orgTherm, orgVisc;
  real *acfTherm, *acfVisc;
} TBuf;

global global_data gd;
global cmdline_data cmd;
global io_header_blj hdr;
// gcc11
global global_data_tree gdtree;
global global_data_tree_bljforcecalc gdforce;
global TBuf *tBuf;
//

#endif /* !_global_defs_h */

