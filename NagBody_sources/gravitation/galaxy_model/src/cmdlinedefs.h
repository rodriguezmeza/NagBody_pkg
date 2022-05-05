/* ==============================================================================
	HEADER: cmdlinedefs.h		[galaxy_models]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "cmdline_defs.h"
	Use in routines and functions: (main)
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: June 06, 2008;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved.
=================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
===============================================================================*/

#ifndef _cmdlinedefs_h
#define _cmdlinedefs_h

//#include "../../../General_libs/general/stdinc.h"

#define HEAD1	"NagBody"
#define HEAD2	"galaxy_model"
#define HEAD3	"Initial condition generator of several types of galaxies"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",				";Parameter input file. Overwrite what follows",
//
    "filename=spiral",			";Output file of N-body frames", ":o",
    "filenamefmt=snap-bdhg",	";Output file format of N-body frames", ":ofmt",
    "statfile=stats",           ";Output file for stat of N-body frames",
    "seed=123",                 ";Random number seed for test run",
// Bulge's parameters ----------------------------------------------------------
	"usebulge=true",			";If true, include bulge",
    "nbulge=1364",				";Number of bodies for bulge",
	"bulgmass=0.333",             ";Bulge mass (Disk mass = 1)",
	"abulge=0.1",				";Bulge scale - length (Disk h = 1)",
	"selfgbul=true",			";If true, include bulge self-gravity",
	"rmaxbulg=0.5",				";Maximum radius for bulge particles",
	"epsbulge=0.01",			";Softening length for bulge particles",
	"axibulge=false",			";If true, the bulge is non-spherical",
	"cbulge=0.075",				";c-bulge",
	"zmaxbulg=0.1",				";zmax-bulge (suggest c*rmaxbulg/a)",
    "nsimpson=100",				";nsimpson",
	"bulgerot=false",			";If true, include bulge rotation",
	"brotfrac=0.75",			";Fraction of reversed particles",
// Disk's parameters -----------------------------------------------------------
	"usedisk=true",				";If true, include disk",
    "ndstars=4096",				";Number of bodies for disk",
	"z0=0.2",					";z-disk length scale",
	"rsolar=8.5/3.5",			";Solar radius",
	"qsolar=1.5",				";Q-solar",
	"epsdisk=0.01",				";Softening length for disk particles",
	"zmax=2.0",					";zmax (suggest 10*z0)",
	"rmax=15",					";Maximum radius for disk particles (Suggest 10*h)",
	"usegas=false",				";Gas in the disk?",
	"ndgas=512",				";Number of disk gas-particles",
	"gasmass=1",				";Total mass of the gas (= 1)",
	"gastemp=10000",			";Temperature of gas (Suggest 10^4)",
	"z0gas=0.2",				";z-gas disk length scale",
	"zmaxgas=2.0",				";zmax gas (suggest 10*z0gas)",
	"rmaxgas=15",				";Maximum radius for gas disk particles (Suggest rmax)",
	"rmingas=0.15",				";Minimum radius for gas disk particles",
	"selfggas=false",			";If true, include gas self-gravity",
// Halo's parameters -----------------------------------------------------------
	"usehalo=true",				";If true, include halo",
	"selfghal=true",			";If true, include halo self-gravity",
	"rmaxhalo=100",				";Maximum radius of halo",
	"nhalo=21844",				";Number of halo particles",
	"epshalo=0.01",				";Softening length for halo particles",
	"halotype=LH",				";Type of halo",
	"halomass=5.333",			";Halo mass (Disk mass = 1)",
	"gamhalo=1",				";Halo core radius (gamma) (Disk h = 1)",
	"rthalo=0.5",				";Halo tidal radius (Disk h = 1)",
	"ahalo=0.5",				";Halo scale-length (Disk h = 1)",
// Satellite's parameters ------------------------------------------------------
	"usesat=false",				";If true, include satellite",
	"satmass=1",				";Satellite mass (Disk mass = 1)",
	"asat=1",					";Satellite scale-length (Disk h = 1)",
	"xsat=50",					";Satellite x-pos",
	"ysat=0",					";Satellite y-pos",
	"zsat=0",					";Satellite z-pos",
	"vxsat=-10",				";Satellite x-vel",
	"vysat=0",					";Satellite y-vel",
	"vzsat=0",					";Satellite z-vel",
	"rmaxsat=15",				";Maximum satellite radius (Disk h = 1)",
	"selfgsat=false",			";If true, include satellite self-gravity",
	"nsat=1024",				";Number of satellite particles",
	"epssat=0.01",				";Softening length of satellite particles",
// Add models composition's parameters -----------------------------------------
	"addmods=false",			";If true, add two BDH galaxies",
	"rp=10",					";Perincenter distance for parabolic orbit",
	"rsep=75",					";Initial center of mass separation",
	"thetmod1=0",				";Rotation angle (degrees) theta for disk 1",
	"phimod1=0",				";Rotation angle (degrees) phi for disk 1",
	"thetmod2=45",				";Rotation angle (degrees) theta for disk 2",
	"phimod2=45",				";Rotation angle (degrees) phi for disk 2",
// ---------------------------------------------------------------------
 	"outpteps=false",			";If true, output eps",
    "Version=0.2",              ";M.A. Rodriguez-Meza 1999-2008",
    NULL,
};

#endif /* ! _cmdlinedefs_h */
