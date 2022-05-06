/* ==============================================================================
!	HEADER: cmdline_defs.h														!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date: February 2005												!
!	Purpose: Definitions for importing arguments from the command line			!
!	Language: C																	!
!	Use: '#include "...."														!
!	Use in routines and functions: md_lj_tree (main)							!
!	External headers: None														!
!	Comments and notes:															!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
!		e-mail: marioalberto.rodriguez@inin.gob.mx
!		http://www.astro.inin.mx/mar											!
!																				!
!	Major revisions:															!
!	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.							!
!===============================================================================
!	Legal matters:																!
!	The author does not warrant that the program and routines it contains		!
!	listed below are free from error or suitable for particular applications,	!
!	and he disclaims all liability from any consequences arising from their		!
!	use.																		!
!==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code for the evolution of a Lennard-Jones liquid"
#define HEAD3	"Tree force calculation"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",		";Parameter input file. Overwrite what follows",
    "forcecalc_method=normal",	";Force calculation method to use",
	"theta=1.0",				";Force accuracy parameter (without use by now)",
	"usequad=false",			";If true, use quad moments (without use by now)",
	"temperature=0.7143",		";Temperature of the simulation",
	"density=0.844",			";Density of the liquid",
	"stepEquil=100",			";step to begin equilibrium computations",
	"stepVel=5",				";number of step jumps to save a velocity histogram",
	"stepAvgVel=4",				";number of velocity histograms to average",
	"sizeHistVel=50",			";array size for velocity histogram",
	"rangeVel=4.0",				";range of velocities for histogram",
	"stepRdf=5",				";number of step jumps to save a RDF histogram",
	"stepAvgRdf=20",			";number of RDF histograms to average",
	"sizeHistRdf=200",			";array size for RDF histogram",
	"rangeRdf=4.0",				";range of RDF for histogram",
    "nbody=512",	 			";Number of bodies for test run",
    "dtime=1/256",              ";Integration time step",
    "tstop=2.0",                ";Time to stop integration",
    "seed=123",                 ";Random number seed for test run",
    "icfile=",                  ";N-Body initial conditions",
    "icfilefmt=",               ";N-Body initial conditions file format",
    "snapout=",                 ";Output file of N-body frames",
    "snapoutfmt=",              ";Output file format of N-body frames",
    "dtout=5/256",              ";Data output time step",
    "dtoutinfo=5/256",          ";Info output time step",
    "statefile=",               ";Write run state to a file",
    "restorefile=",             ";Continue run from state file",
    "options=",                 ";Various control options",
    "Version=1.2",              ";Mar 2005-2006",
    NULL,
};

#endif /* ! _cmdline_defs_h */
