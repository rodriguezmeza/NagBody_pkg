/*==============================================================================
    HEADER: cmdline_defs.h		[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date: January 2018
    Purpose: Definitions for importing arguments from the command line
    Language: C
    Use: '#include "cmdline_defs.h"
    Use in routines and functions: (main)
    External headers: stdinc.h
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
    and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"nchi2 fitting models code"
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
"paramfile=",               ";Parameter input file. Overwrite what follows",
// Data options
"inputObsDataFile=",		";File to process", ":in",
"usingcolumns=1:2:3",       ";Observed data file using columns (x,y,sigma)", ":uc",
"witherrors=true",          ";Plot error bars", ":werr",
"errorstype=1",             ";Error bars type", ":errt",
// Model options:
"modeltype=",               ";Model to study", ":model",
"suffixModel=",             ";Suffix model to add to output filenames", ":suffix",
"modelParamfile=",          ";Parameters of the model can be given through this file", ":mpf",
// Output options:
"xoutmin=",             ";xmin for the output model table",
"xoutmax=",             ";xmax for the output model table",
"Nx=",                  ";Number of x of the output model table", ":nx",
//
// Minimization parameters:
"minMethod=simplex",         ";Minimization method to use", ":mm",
"gridN=10",                  ";Number of elements of the gridding to search an initial point in parameter space", ":gn",
"gridNSimplex=50",            ";Number that control the size of the gridding of the initial point in parameter space (only in simpex method)", ":gnsimplex",
"ftolMinimization=1.0e-6",                  ";Fractional convergence tolerance to achive in the minimization search", ":ftolm",
"useCentralParamValues=true",   ";Use central values of parameteres in the prior region of parameter space", ":ucparams",
// Other options:
"epsQuad=1.0e-6",         ";Parameter error in the numerical quadratures", ":epsq",
//
"options=",                 ";Various control options", ":opt",
"Version=0.3",              ";Mario A. Rodríguez-Meza 2005-2020",
NULL,
};

#endif // ! _cmdline_defs_h
