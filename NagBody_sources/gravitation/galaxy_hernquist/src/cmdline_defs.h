/*==============================================================================
	HEADER: cmdline_defs.h			[galaxy_hernquist]
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
#define HEAD2	"Code to generate a galaxy model"
#define HEAD3	"Hernquist model (by Dubinski)"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			";Parameter input file. Overwrite what follows",
//
    "ngas=0",               ";Number of bodies for gas",
    "nhalo=245760",			";Number of bodies for halo",
    "ndisk=29491",			";Number of bodies for disk",
    "nbulge=1024",			";Number of bodies for bulge",
//
    "mgas=0.",              ";Mass of the gas",
    "mhalo=1.6",            ";Mass of the halo",
    "mdisk=0.1017823",      ";Mass of the disk",
    "mbulge=0.0025",        ";Mass of the bulge",
//
    "masscut=0.95",        ";Mass fraction on which the distribution is cut",
//
    "ag=0.1",               ";Radial length scale of the gas",
    "zg=0.01",              ";z length scale of the gas",
    "rmaxg=0.",             ";Maximu radial length of the gas",
//
    "ad=22.2",               ";Radial length scale of the disk",
    "zd=0.007",              ";z length scale of the disk",
    "rmaxd=0.5",             ";Maximu radial length of the disk",
//
    "ah=0.3",               ";Radial length scale of the halo",
    "gammah=0.0",              ";z length scale of the halo",
//    "rmaxh=0.",             ";Maximu radial length of the halo",
//
    "ab=0.008",             ";Radial length scale of the bulge",
    "gammab=0.0",           ";z length scale of the bulge",
//    "rmaxb=0.",             ";Maximu radial length of the bulge",
//
    "snapout=galaxy.bin",   ";Output file of N-body frames", ":o",
    "snapoutfmt=",			";Output file format of galaxy model (default is gadget)", ":ofmt",
    "seed=-1",				";Random number seed for gsl_init",
    "options=",				";Various control options", ":opt",
    "Version=0.2",			";Mar 2005-2014",
    NULL,
};

#endif // ! _cmdline_defs_h
