/*==============================================================================
	HEADER: cmdline_defs.h			[nbody_n2]
	Written by: Mario A. Rodriguez-Meza
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
        http://www.inin.gob.mx/

	Major revisions: July 2007;  November 2008; april 2018
	Copyright: (c) 2005-2018 Mar.  All Rights Reserved
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code for the evolution of an N-Body selfgravitating system"
#define HEAD3	"Direct force calculation"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			";Parameter input file. Overwrite what follows",
    "eps=0.025",			";Force smoothing length",
    "nbody=1024",			";Number of bodies for test run",
    "dtime=1/32",			";Integration time step", ":dt",
    "tstop=2.0",			";Time to stop integration",
    "seed=123",				";Random number seed for test run",
    "icfile=",				";N-Body initial conditions", ":ic",
    "icfilefmt=",			";N-Body initial conditions file format", ":icfmt",
    "snapout=",				";Output file of N-body frames", ":o",
    "snapoutfmt=",			";Output file format of N-body frames", ":ofmt",
    "dtout=4/32",			";Data output time step",
    "dtoutinfo=4/32",		";Info output time step",
    "statefile=",			";Write run state to a file", ":state",
	"stepState=20",			";number of steps to save a state-run file",
    "restorefile=",			";Continue run from state file", ":restore",
    "options=",				";Various control options", ":opt",
    "Version=0.2",			";Mario A. Rodriguez-Meza 2005-2018",
    NULL,
};

#endif // ! _cmdline_defs_h
