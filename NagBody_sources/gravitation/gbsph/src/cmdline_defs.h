/* =============================================================================
	HEADER: cmdline_defs.h				[gbsph]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: main.c
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: February 2005; July 24, 2007; October 04, 2007;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code for evolution of a General Body and SPH (gbsph)"
#define HEAD3	"Hierarchical force calculation"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",				";File with the info of parameters to use",
    "forcecalc_method=barnes",	";Gravitation calculation method to use", ":fmeth",
    "force_models=",			";Forces to include", ":fmod",
    "icfile=",                  ";Input file with N-Body initial conditions", ":ic",
    "icfilefmt=",               ";N-Body initial conditions file format", ":icfmt",
    "snapout=",					";Output file of N-body frames", ":o",
    "snapoutfmt=",				";Snap output file format of N-body system", ":ofmt",
    "dtime=1/32",               ";Integration time step", ":dt",
    "eps=0.025",               	";Density smoothing length",
    "theta=1.0",               	";Force accuracy parameter",
    "usequad=false",            ";If true, use quad moments",
//
    "dm_lambda=1.0",           	";Dark matter lambda, range of interaction",
    "dm_alpha=1.0",             ";Dark matter alpha, intensity of the force",
    "G=1.0",					";Local value of the constant of gravity. Its value sets units",
	"dm_a=1.0",					";Scalar field slope coeficient for perturbation",
	"dm_time=0",				";Scalar field perturbation time",
//
    "eps_pot=1.0",				";External potential energy or mass",
    "sigma_pot=1.0",			";External potential length",
    "x_pot=7.0",				";External potential x-position",
    "y_pot=0.0",				";External potential y-position",
    "z_pot=0.0",				";External potential z-position",
//
	"computeTransport=false",	";Compute transport properties",
//
    "options=",                 ";Various control options", ":opt",
    "tstop=2.0",                ";Time to stop integration",
    "dtout=8/32",				";Data output time step",
    "dtoutinfo=16/32",          ";Info output time step",
    "nbody=4096",               ";Number of bodies for test run",
    "seed=-1",					";Random number seed for test run",
    "statefile=",				";Write state file as code runs", ":state",
	"stepState=20",				";number of steps to save a state-run file",
    "restorefile=",				";Continue run from state file", ":restore",
    "Version=0.2",              ";M.A. Rodriguez-Meza 1999-2014",
    NULL,
};

#endif // ! _cmdline_defs_h
