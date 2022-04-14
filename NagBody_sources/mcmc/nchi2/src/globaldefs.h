/*==============================================================================
    HEADER: globaldefs.h		[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date: January 2018
    Purpose: Definitions of global variables and parameters
    Language: C
    Use: '#include "global_defs.h"
    Use in routines and functions:
    External headers: stdinc.h, data_struc_defs.h
    Comments and notes:
    Info: Mario A. Rodriguez-Meza
        Depto. de Fisica, ININ
        Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
        e-mail: marioalberto.rodriguez@inin.gob.mx
        https://github.com/rodriguezmeza

    Major revisions:
    Copyright: (c) 2005-2020 Mar.  All Rights Reserved
================================================================================
    Legal matters:
    The author does not warrant that the program and routines it contains
    listed below are free from error or suitable for particular applications,
    and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h


//========================================================
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef NOGNU
#include "general_libs/general/stdinc.h"
#include "general_libs/math/vectmath.h"
#include "general_libs/general/getparam.h"
#include "general_libs/math/mathfns.h"
#include "general_libs/io/inout.h"
#include "general_libs/general/machines.h"
#include "general_libs/math/minpack.h"
#include "general_libs/math/linpack.h"
#include "general_libs/math/numrec.h"
#include "general_libs/math/diffeqs.h"
#include "general_libs/math/mathutil.h"
// #include <strings.h>                            // For unix
#include "general_libs/general/strings.h"    // For Visual c
#include "general_libs/math/quads.h"
#else
#include "stdinc.h"
#include "vectmath.h"
#include "getparam.h"
#include "mathfns.h"
#include "inout.h"
#include "machines.h"
#include "minpack.h"
#include "linpack.h"
#include "numrec.h"
#include "diffeqs.h"
#include "mathutil.h"
// #include <strings.h>							// For unix
#include "strings.h"	// For Visual c
#include "quads.h"
#endif


#include "data_struc_defs.h"
#include "constants_defs.h"
#include "models.h"

//========================================================

typedef struct {
    string paramfile;
// Data options
    string inputObsDataFile;
    string obsdatafile;
    string usingcolumns;
    string witherrors;
    int errorstype;
// Model options:
    string model;
    string suffixModel;
    string model_paramfile;
// Output options:
    string xoutmin;
    string xoutmax;
    string Nx;
// Minimization parameters:
    string minMethod;
    int gridN;
    int gridNSimplex;
    real ftolmin;
    bool ucparams;
// Other options:
    real epsq;
	string options;
} cmdline_data, *cmdline_data_ptr;

#define MAXLINES    1000

typedef struct {
	real cpuinit;

	string headline0;
	string headline1;
	string headline2;
	string headline3;

    char model_comment[100];
    real (*modelfun)(real, real *);
    void (*modelfunwd)(double, double [],
                  double *, double [], int);
    real (*modelChi2)(real (*ymodel)(real, real *), real params[]);
    void (*modelend)(void);
    
// Input options
    int nfiles;
    char *filenames[40];
    
    int filecount;
    
    int nwitherrorbars;
    bool errorbars[MAXLINES];

    // ----- Block of definitions for usingcolumns ---------------------------------
    int *vcol1;
    int *vcol2;
    int *vcol3;
    int *vcol4;

//
// Output options:
    real xoutmin;
    real xoutmax;
    int Nx;
//

    char tmpDir[100];
    char outDir[100];
    char logfilePath[100];
    char fpfnamemodeltest[100];
    char fpfnamemodeltable[100];

    int minmethod_int;
    char minmethod_comment[100];

    int col1;
    int col2;
    int col3;
    int col4;

// Structure to add an interpolation mechanism of the model function
    int nInterp;
    real *xInterp;
    real *yInterp;
    real *yInterp2;

	FILE *outlog;

	char mode[2];

} global_data, *global_data_ptr;

global global_data gd;
global cmdline_data cmd;

// Mover a global_data:
global real *yout;

typedef struct _ParamsTable {
    int nparams;
    real *params;
    real *dparams;
    real **rparams;
} ParamsTable, *ParamsTableptr;

global ParamsTable pd;

// Obs Data structure
typedef struct _ObsDataTable {
    int no;
    real *xo;
    real *yo;
    real *sigma;
    real *sigmam;
} ObsDataTable, *ObsDataTableptr;

ObsDataTableptr *pObs;
ObsDataTable pObsT;     // Combined totals

#define nObs(x)    (((ObsDataTableptr) (x))->no)
#define xObs(x)    (((ObsDataTableptr) (x))->xo)
#define yObs(x)    (((ObsDataTableptr) (x))->yo)
#define sigmaObs(x)    (((ObsDataTableptr) (x))->sigma)
#define sigmamObs(x)    (((ObsDataTableptr) (x))->sigmam)


#include "protodefs.h"

#endif // ! _globaldefs_h

