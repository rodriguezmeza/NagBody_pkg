/* ==============================================================================
	HEADER: global_defs.h			[analysis_grav]
	Written by: Mario A. Rodriguez-Meza
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
		e-mail: marioalberto.rodriguez@inin.gob.mx
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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/types.h>

#include <string.h>

#ifndef NOGNU
#include "./general_libs/stdinc.h"
#include "./general_libs/vectdefs.h"
#include "./general_libs/nagbody.h"
#include "./general_libs/mathfns.h"
#include "./general_libs/mathutil.h"
#include "./general_libs/vectmath.h"
#include "./general_libs/getparam.h"
#include "./general_libs/inout.h"
#include "./general_libs/constant.h"
/* #include <strings.h> */							// For unix
#include "./general_libs/strings.h"	// For Visual C
#include "./general_libs/pldefs.h"
#include "./general_libs/lic.h"
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


typedef struct {
	string paramfile;

//	real temperature;
//	real density;

//	string nbodyprop;
//	string massprop;
//	string LxLyprop;

//	real eps11;
//	real eps12;
//	real eps22;
//	real sigma11;
//	real sigma12;
//	real sigma22;
//	real Rcut11;
//	real Rcut12;
//	real Rcut22;

//	int stepEquil;

/*
	int stepAvgRhoAxes;
	int sizeHistRhoAxes;

	int stepAvgNFrecAxes;
	int sizeHistNFrecAxes;

	int stepAvgVel;
	int sizeHistVel;
	real rangeVel;

	int stepAvgRdf;
	int sizeHistRdf;
	real rangeRdf;
*/

	int stepAvgRhoTheta;
	int sizeHistRhoTheta;
	real RhoDeltaZ;
	real RhoR;
	real RMax;
	real RhoDeltaR;
	real ThetaMin;
	real ThetaMax;

	real xmin;
	real xmax;
	real ymin;
	real ymax;
	real zmin;
	real zmax;

	int stepAvgRho;
	int sizeHistRho;

	int stepAvgVcR;
	int sizeHistVcR;
	real rangeR;

//	bool computeTransport;
//	int stepAcf;
	int limitAcfAv;
	int nBuffAcf;
	int nValAcf;

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

	string options;

	string data_analysis_type;

	string bodiesID;
	string bodiesSets;

	real reductionFac;
	string foffile;
	real gravityConstant;

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

// Seran activados cuando se necesite calcular la fuerza entre particulas.
//	bool usequad;
//	real eps;
//	real theta;

} cmdline_data, *cmdline_data_ptr;

#define MAXBODIESID		20
#define MAXBODIESSETS	20
#define MAXYRANGESETS	4

typedef struct {
	int nbody;
/*
	int nbody1;
	int nbody2;
	real mass1;
	real mass2;
	real Lx;
	real Ly;
#ifdef THREEDIM
	real Lz;
#endif
*/
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
//	vector Box;

/*
	int countRhoAxes;
	realptr histRhoX;
	realptr histRhoY;
	realptr histRhoZ;

	int countNFrecAxes;
	int *histNFrecX;
	int *histNFrecY;
	int *histNFrecZ;

	int countVel;
	realptr histVel;

	int countRdf;
	realptr histRdf;
*/

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

	int countRho;
	realptr histRho;

	int countVcR;
	realptr histVcR;
	realptr histMass;
	realptr histForce;
	realptr histNPart;
	realptr histRVelDisp;
	realptr histVrR;

	int countAcfAv;
	realptr avAcfVel;
	real intAcfVel;
	realptr avAcfTherm;
	realptr avAcfVisc;
	real intAcfTherm;
	real intAcfVisc;

//	string snapfilename;

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

/*
	int Rdf_flag;
	int Vel_flag;
	int RhoAxes_flag;
	int NFrecAxes_flag;
*/

	int RhoTheta_flag;
	int Rho_flag;
	int VcR_flag;
	int Acf_flag;
/*
	int nbbcalc;
	int nbccalc;
	cellptr root;
	real rsize;
	double cpuforce;
	double cputree;
	int tdepth;
	int ncell;
*/

	int in_long_fmt;
//	int out_long_fmt;			// por el momento no es necesario ....

} global_data, *global_data_ptr;

typedef struct {
  vector *orgVel;
  real *acfVel;
  int count;
  vector orgTherm, orgVisc;
  real *acfTherm, *acfVisc;
} TBuf;

//#undef MAXBODIESID
//#undef MAXBODIESSETS

//global bodyptr bodytab;
global global_data gd;
global cmdline_data cmd;
global io_header_blj hdr;
// gcc11
global global_data_tree gdtree;
global global_data_tree_bljforcecalc gdforce;
global global_data_treegrav gdtreegrav;
global TBuf *tBuf;
//

#endif /* !_global_defs_h */

