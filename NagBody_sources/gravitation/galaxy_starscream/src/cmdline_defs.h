/*==============================================================================
	HEADER: cmdline_defs.h			[galaxy_starscream]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: main
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 2007;  November 2008;
	Copyright: (c) 2005-2014 Mar.  All Rights Reserved
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code to generate a galaxy pair collision model"
#define HEAD3	"Starscream (by Jay Billings)"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			";Parameter input file. Overwrite what follows",
//
    "ingal1=galaxy1",       ";Input file with galaxy1 data",
    "ingal2=",              ";Input file with galaxy2 data",
//
    "a=0.0",                ";Euler alpha angle",
    "b=0.5236",             ";Euler beta angle (default pi/6)",
    "g=0.0",                ";Euler gamma angle",
//
    "radius=1.492945497",   ";Radius of the orbit",
    "p=0.1",                ";Impact parameter",
//
    "snapout=galaxies",     ";Output file of N-body frames", ":o",
    "snapoutfmt=",			";Output file format of galaxy model (default is gadget)", ":ofmt",
    "options=",				";Various control options", ":opt",
    "Version=0.2",			";Mar 2005-2014",
    NULL,
};

#endif // ! _cmdline_defs_h
