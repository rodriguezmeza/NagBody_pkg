/*==============================================================================
	HEADER: globaldefs.h			[md_blj_n2]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: md_lj_n2 (main), start_run, time_step,
					tree_ljforcecalc, md_lj_io
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:  November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.																		!
==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h

#include "./general_libs/general/stdinc.h"
#include "./general_libs/math/vectdefs.h"
#include "./general_libs/NagBody/nagbody.h"

#include <stdio.h>

typedef struct {										// CHECK 2D --- OK!!!
	string paramfile;
	string forcecalc_method;
	int potType;
	string fnamePot;

	real temperature;
	real density;
	bool adjustTemperature;
	int stepAdjustTemperature;
	bool adjustCenterOfMass;

	int stepEquil;
	int stepAvg;
//	bool printSnap;
//	int stepSnap;

	bool computeRhoAxes;
	int stepAvgRhoAxes;
	int stepRhoAxes;
	int sizeHistRhoAxes;

	bool computeNFrecAxes;
	int stepAvgNFrecAxes;
	int stepNFrecAxes;
	int sizeHistNFrecAxes;

	bool computeVelDist;
	int stepAvgVel;
	int stepVel;
	int sizeHistVel;
	real rangeVel;

	bool computeRdf;
	int stepAvgRdf;
	int stepRdf;
	int sizeHistRdf;
	real rangeRdf;

	bool computePressAxes;
	int stepAvgPressAxes;
	int stepPressAxes;
	int sizeHistPressAxes;

	bool computeChemPot;
	int stepAvgChemPot;
	int stepChemPot;
	int sizeHistChemPot;
	int numTestBodies;

	bool computeDiffusion;
	int stepAvgDiffuse;
	int nBuffDiffuse;
	int nValDiffuse;
	int stepDiffuse;

	bool computeVelAcf;
//	int limitAcfAv;
//	int nBuffAcf; 
//	int nValAcf; 
//	int stepAcf;

	bool computeBulkViscosity;

	bool computeTransport;
	int stepAcf;
	int stepAvgAcf;
	int nBuffAcf;
	int nValAcf;

	bool computeSTCorr;
	int stepCorr;
	int stepAvgCorr;
	int nBuffCorr;
	int nValCorr;
	int nFunCorr;

	real lattCorr_kx;
	real lattCorr_ky;
#if (NDIM==3)
	real lattCorr_kz;
#endif

	string nbodyprop;
	string massprop;
#if (NDIM==3)
	string LxLyprop;
#else
	string Lx;
#endif

	real eps11;
	real eps12;
	real eps22;
	real sigma11;
	real sigma12;
	real sigma22;
	real Rcut11;
	real Rcut12;
	real Rcut22;

	string dtimestr;
	real tstop;
	int icModel;
	int intMethod;
	int seed;

	string icfile;
	string icfilefmt;
	string snapoutfile;                 
	string snapoutfilefmt;
	string dtoutstr;
	string dtoutinfostr;
	string statefile;
	int stepState;
	string restorefile;

	string options;            

	string unitCells;

} cmdline_data, *cmdline_data_ptr;

typedef struct {										// CHECK 2D --- OK!!!
	int nbody;

	real dtime;
	real dtout;
	real dtoutinfo;

	int forcemethod_int;

	real cputotforce;
	real cpuinit;
	real cputotout;
	real cputotal;
	real cpuchempot;

	string headline0;
	string headline1;
	string headline2;
	string headline3;

	char model_comment[100];

	string headerfmt;

	real tnow;
	real tout;
	real toutinfo;
	int nstep;

	FILE *outlog;

	int stopflag;

	real eps;
	real sigma;
	real mass;
	real kB;

	real kinEnergySave;

	real kinEnergy;
	real potEnergy;
	real totEnergy;
	real pressure;
	real vSum;
	real vvSum;

	real sKinEnergy;
	real ssKinEnergy;
	real sPotEnergy;
	real ssPotEnergy;
	real sTotEnergy;
	real ssTotEnergy;
	real sPressure;
	real ssPressure;

	int countRhoAxes;
	realptr histRhoX;
	realptr histRhoY;

	realptr histRhoX1;
	realptr histRhoY1;

	realptr histRhoX2;
	realptr histRhoY2;

#if (NDIM==3)
	realptr histRhoZ;
	realptr histRhoZ1;
	realptr histRhoZ2;
#endif

	int countNFrecAxes;
	int *histNFrecX;
	int *histNFrecY;

	int *histNFrecX1;
	int *histNFrecY1;

	int *histNFrecX2;
	int *histNFrecY2;

	int *histNFrecXD;
	int *histNFrecYD;

#if (NDIM==3)
	int *histNFrecZ;
	int *histNFrecZ1;
	int *histNFrecZ2;
	int *histNFrecZD;
#endif

	int countVel;
	realptr histVel;

	int countRdf;
	realptr histRdf;
	realptr histRdf11;
	realptr histRdf12;
	realptr histRdf22;

	int countPressAxes;
	realptr histPressX;
	realptr histPressY;
#if (NDIM==3)
	realptr histPressZ;
#endif

	int countChemPot;
	realptr histChemPotX;
	realptr histChemPotY;
#if (NDIM==3)
	realptr histChemPotZ;
#endif
	real ChemPot;
	real sChemPot;
	real ssChemPot;

	real *rrDiffuseAv;
	int countDiffuseAv;

//	real *avAcfVel
//	real intAcfVel;
//	int countAcfAv;

	real *avAcfdP;
	real intAcfdP;
	real sPressureSave;

	int countAcfAv;
	realptr avAcfVel;
	real intAcfVel;
	realptr avAcfTherm;
// - Bulk Viscosity
	real PV_K;
	real PV_P;
	real sPV_K;
	real sPV_P;
	real ssPV_K;
	real ssPV_P;
	real PV_KAvg;
	real PV_KAvgSave;
	real PV_PAvg;
	real PV_PAvgSave;
	realptr avAcfBVisc_KK;
	realptr avAcfBVisc_KP;
	realptr avAcfBVisc_PP;
	real intAcfBVisc_KK;
	real intAcfBVisc_KP;
	real intAcfBVisc_PP;
// - Shear Viscosity
	realptr avAcfVisc_KK;
	realptr avAcfVisc_KP;
	realptr avAcfVisc_PP;
	real intAcfVisc_KK;
	real intAcfVisc_KP;
	real intAcfVisc_PP;
//
	realptr avAcfVisc;

	real intAcfTherm;
//
//
	real intAcfVisc;

	int countCorrAv;
	real **avAcfST;
	real *valST;

//	int stepSnapInit;	// borrar cuando ya este seguro de que no sirven...
//	int stepSnap;
//	int stepAvg;

	int stepAvgOld;				// In order to restore work ...
	int nstepOld;				// In order to restore work ...
	int nstepNew;				// In order to restore work ...
	bool stepAvgFlag;			// In order to restore work ...
	bool stepAvgStatus;			// In order to restore work ...

	real vMag;

	real lAvg;
	real dtcrit;
	real fAvg;
	real phiAvg;

	char mode[2];

	int nbody1;
	int nbody2;

	real Lx;
	real Ly;
#ifdef THREEDIM
	real Lz;
#endif

	real Rcut11Min;
	real Rcut22Min;
	real RcutAllMax;
	real RcutAllMin;

	real RcutSq11Max;
	real RcutSq22Max;
	real RcutSq11Min;
	real RcutSq22Min;
	real RcutSqAllMax;
	real RcutSqAllMin;

	real mass1;
	real mass2;

	real LBoxMin;
	real VelMax;
	real TimeMin;

	vectorI numUcell;
	int nbc1;
	int nbc2;
	real latticeCorr;

} global_data, *global_data_ptr;

typedef struct {
  vector *orgR, *rTrue;
  real *rrDiffuse;
  int count;
} TBufDiffusion;

typedef struct {
  vector *orgVel;
  real *acfVel;
  int count;
} TBufVelAcf;

typedef struct {
  real *orgdP;
  real *acfdP;
  int count;
} TBufdPressAcf;

typedef struct {										// CHECK 2D --- OK!!!
  vector *orgVel;
  real *acfVel;
  int count;
  vector orgTherm, orgVisc;
  real *acfTherm, *acfVisc;
// - Bulk Visocisty
  vector orgBVisc_K;
  vector orgBVisc_P;
  real *acfBVisc_KK;
  real *acfBVisc_KP;
  real *acfBVisc_PP;
// - Shear Viscosity
  vector orgVisc_K;
  vector orgVisc_P;
  real *acfVisc_KK;
  real *acfVisc_KP;
  real *acfVisc_PP;
} TBuf;

typedef struct {										// CHECK 2D --- OK!!!
  real **acfST, *orgST;
  int count;
} TBufCorr;

global global_data gd;
global cmdline_data cmd;
global io_header_blj hdr;

// gcc11
global global_data_tree gdtree;
global global_data_tree_bljforcecalc gdforce;

global TBufDiffusion *tBufD;
global TBufVelAcf *tBufVAcf;
global TBufdPressAcf *tBufdPAcf;
global TBuf *tBuf;
global TBufCorr *tBufC;
//

#endif /* ! _globaldefs_h */

