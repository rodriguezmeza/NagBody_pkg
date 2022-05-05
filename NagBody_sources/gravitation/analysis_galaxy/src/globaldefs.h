/* ==============================================================================
	HEADER: global_defs.h			[analysis_galaxy]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main), start_run, time_step,
					forcecalc, md_lj_tree_io
	External headers: data_struct_defs
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
===============================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _global_defs_h
#define _global_defs_h

//====================================================
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/types.h>

#include <string.h>

#ifndef NOGNU
#include "./general_libs/general/stdinc.h"
#include "./general_libs/math/vectdefs.h"
#include "./general_libs/NagBody/nagbody.h"
#include "./general_libs/math/mathfns.h"
#include "./general_libs/math/mathutil.h"
#include "./general_libs/math/vectmath.h"
#include "./general_libs/general/getparam.h"
#include "./general_libs/io/inout.h"
#include "./general_libs/general/constant.h"
/* #include <strings.h> */							// For unix
#include "./general_libs/general/strings.h"	// For Visual C
#include "./general_libs/visual/pldefs.h"
#include "./general_libs/general/lic.h"
#else
#include "stdinc.h"
#include "vectdefs.h"
#include "nagbody.h"
#include "mathfns.h"
#include "mathutil.h"
#include "vectmath.h"
#include "getparam.h"
#include "inout.h"
#include "constant.h"
/* #include <strings.h> */							// For unix
#include "strings.h"	// For Visual C
#include "pldefs.h"
#include "lic.h"
#endif

#include "protodefs.h"
//====================================================
//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/vectdefs.h"
//#include "../../../General_libs/NagBody/nagbody.h"


typedef struct {
	string paramfile;

	int stepAvgRhoTheta;
	int sizeHistRhoTheta;
	real RhoDeltaZ;
	real RhoR;
	real RMax;
	real RhoDeltaR;
	real ThetaMin;
	real ThetaMax;

	int stepAvgVcR;
	int sizeHistVcR;
	real rangeR;

	int limitAcfAv;
	int nBuffAcf;
	int nValAcf;

	string cmpos1;
	string cmvel1;
	string cmpos2;
	string cmvel2;

	string in;
	string infmt;
	string out;                 
	string outfmt;
	string basedir;

	int isnap;
	int fsnap;

	string options;

	string data_analysis_type;

	string bodiesID;
	string bodiesSets;

	string xrange;
	string yrange;
	string zrange;

	string usingcolumns;
	string usingrows;

	string xlabel;
	string ylabel;
	real labelfontsize;
	string plotlabel;
	bool withdots;
	bool withsymbols;
	int symboltype;
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

	bool usequad;
	real eps;
	real theta;

} cmdline_data, *cmdline_data_ptr;

#define MAXBODIESID		20
#define MAXBODIESSETS	20
#define MAXYRANGESETS	4

typedef struct {
	int nbody;

	string headline0;
	string headline1;
	string headline2;
	string headline3;
	char model_comment[100];

	string headerfmt;

	char *filenames[40];
	char *filenamesfmt[40];
	vector cmpos1;
	vector cmvel1;
	vector cmpos2;
	vector cmvel2;

	real tnow;
	real tout;
	int nstep;
	real cpuinit;
	FILE *outlog;
	real mass;

	int countRhoTheta;
	realptr histRho0Theta;
	realptr histRho1Theta;
	realptr histRho2Theta;

	realptr histRho0ThetaSave;				// De hecho no son necesarios
	realptr histRho1ThetaSave;				// A la hora de usarlos se puede
	realptr histRho2ThetaSave;				// intercambiar el orden de guardar...

	realptr histAccRTheta;
	realptr histAccTTheta;
	realptr histVelTTheta;

	int countVcR;
	realptr histVcR;

	int countAcfAv;
	realptr avAcfVel;
	real intAcfVel;
	realptr avAcfTherm;
	realptr avAcfVisc;
	real intAcfTherm;
	real intAcfVisc;

	int nbodiesID;
	int bodyID[MAXBODIESID];

	int nbodiesSets;
	int bodyIDMin[MAXBODIESSETS];
	int bodyIDMax[MAXBODIESSETS];

	int nyrangeSets;
	real yrangeMin[MAXYRANGESETS];
	real yrangeMax[MAXYRANGESETS];

	bool x_autoscale;
	bool y_autoscale;
	bool z_autoscale;
	real xmin;
	real xmax;
	real ymin;
	real ymax;
	real zmin;
	real zmax;

	int column1;
	int column2;
	int column3;
	int row1;
	int row2;

	int RhoTheta_flag;
	int VcR_flag;
	int Acf_flag;
} global_data, *global_data_ptr;

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
// STATIC problem: gcc version 11
global global_data_tree gdtree;
global global_data_tree_bljforcecalc gdforce;
global global_data_treegrav gdtreegrav;
global TBuf *tBuf;
//

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


#endif /* !_global_defs_h */

