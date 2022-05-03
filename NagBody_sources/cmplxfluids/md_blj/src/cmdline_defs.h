/*==============================================================================
	HEADER: cmdline_defs.h				[md_blj]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: md_blj (main)
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:  November 2008;
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code for the evolution of a General binary mixture liquid"
#define HEAD3	"Tree force calculation"
														// CHECK 2D --- OK!!!
string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",				";Parameter input file. Overwrite what follows",
    "forcecalc_method=cells",	";Force calculation method to use", ":fm",
    "potentialType=0",			";Force model type to use (0 - LJ, 1 - SLJ, 2 - table file)", ":pT",
    "fnamePot=",				";Filename with table force and potential (model type 2 only)",
	"theta=1.0",				";Force accuracy parameter (without use by now)",
	"usequad=false",			";If true, use quad moments (without use by now)",
	"density=0.844",			";Density of the liquid", ":d",
	"temperature=0.7143",		";Temperature of the simulation", ":t",
	"adjustTemperature=true",	";Adjust temperature (Microcanonical ensamble if false)",
	"stepAdjustTemperature=1",	";Number of steps to adjust temperature",
	"adjustCenterOfMass=false",	";Adjust center of mass of the system (To compute bulk viscosity)",
	"stepEquil=100",			";step to begin equilibrium computations",
	"stepAvg=10",				";number of steps to average properties",
//
//
	"eps11=1.0",				";Spieces-1 eps",
	"eps12=1.0",				";Cross eps",
	"eps22=1.0",				";Spieces-2 eps",
	"sigma11=1.0",				";Spieces-1 sigma",
	"sigma12=1.0",				";Cross sigma",
	"sigma22=1.0",				";Spieces-2 sigma",
	"Rcut11=1.12246",			";Spieces-1 Cut radius (2^1/6 sigma)",
	"Rcut12=1.12246",			";Cross Cut radius (2^1/6 sigma)",
	"Rcut22=1.12246",			";Spieces-2 Cut radius (2^1/6 sigma)",
    "dtime=1/256",              ";Integration time step", ":dt",
    "tstop=4.0",                ";Time to stop integration",
    "intMethod=2",				";Integration method (0,1 - fixed dt; 2 - variable dt)",
//
// CHATO:
//    "icModel=2",				";IC method (1, 2, 3, 4, 5 and 6). Run with tstop=0 to see details",
// FCC:
    "icModel=6",				";IC method (1, 2, 3, 4, 5 and 6). Run with tstop=0 to see details",
    "seed=123",                 ";Random number seed for test run",
//
#ifdef THREEDIM
	"unitCells=8:8:8",			";Number of unit cells along axes (icModel=3)", ":uC",
#else
	"unitCells=8:8",			";Number of unit cells along axes (icModel=3)", ":uC",
#endif
// CHATO:
//    "nbodyprop=512/512",		";Number of bodies of the mixture for test run (nbody1/nbody2)",
#ifdef THREEDIM
// FCC:
    "nbodyprop=2/2",			";Proportion of number of bodies of the mixture for test run (nbody1/nbody2)",
#else
    "nbodyprop=1/1",			";Proportion of number of bodies of the mixture for test run (nbody1/nbody2)",
#endif
    "massprop=1/1",				";Proportion of masses of the spieces for test run (mass1/mass2)",
#ifdef THREEDIM
	"LxLyprop=8.4653/8.4653",	";Base sides of the parallelepiped",
#else
	"Lx=5",						";x-side of base rectangle",
#endif
//
    "icfile=",                  ";N-Body initial conditions", ":ic",
    "icfilefmt=snap-blj-ascii", ";N-Body initial conditions file format", ":icfmt",
//
    "snapout=",                 ";Output file of N-body frames", ":o",
    "snapoutfmt=snap-blj-ascii",";Output file format of N-body frames", ":ofmt",
    "dtout=5/256",              ";Data output time step",
    "dtoutinfo=5/256",          ";Info output time step",
    "statefile=",               ";Write run state to a file", ":state",
	"stepState=20",				";number of steps to save a state-run file",
    "restorefile=",             ";Continue run from state file", ":restore",
    "options=",                 ";Various control options", ":opt",
//
// OUTPUT:
//
//	"printSnap=false",			";Print snap (column form) every stepSnap=dtout/dtime",	// borrar cuando ya este seguro de que no sirven...
//	"stepSnap=10",				";number of steps to save a snap",	// borrar cuando ya este seguro de que no sirven...
//
	"computeRhoAxes=false",		";Compute density profile",
	"stepRhoAxes=5",			";number of step jumps to save a rho_axes histogram",
	"stepAvgRhoAxes=4",			";number of rho_axes histograms to average",
	"sizeHistRhoAxes=100",		";array size for rho_axes histogram",
//
	"computeNFrecAxes=false",	";Compute frequency distribution of density profile",
	"stepNFrecAxes=5",			";number of step jumps to save a nfrec_axes histogram",
	"stepAvgNFrecAxes=4",		";number of nfrec_axes histograms to average",
	"sizeHistNFrecAxes=100",	";array size for nfrec_axes histogram",
//
	"computeVelDist=false",		";Compute velocity distribution",
	"stepVel=5",				";number of step jumps to save a velocity histogram",
	"stepAvgVel=4",				";number of velocity histograms to average",
	"sizeHistVel=50",			";array size for velocity histogram",
	"rangeVel=4.0",				";range of velocities for histogram",
//
	"computeRdf=false",			";Compute radial distribution function",
	"stepRdf=5",				";number of step jumps to save a RDF histogram",
	"stepAvgRdf=20",			";number of RDF histograms to average",
	"sizeHistRdf=200",			";array size for RDF histogram",
	"rangeRdf=4.0",				";range of RDF for histogram",
//
	"computePressAxes=false",	";Compute pressure profile",
	"stepPressAxes=5",			";number of step jumps to save a pressure profile measurement",
	"stepAvgPressAxes=20",		";number of pressure profile measurements to average",
	"sizeHistPressAxes=10",		";array size for pressure profile histogram",
//
	"computeChemPot=false",		";Compute chemical potential",
	"stepChemPot=5",			";number of step jumps to save a ChemPot measurement",
	"stepAvgChemPot=20",		";number of ChemPot measurements to average",
	"sizeHistChemPot=10",		";array size for chemical potential histogram",
	"numTestBodies=20",			";number of test bodies to make a ChemPot measurement",
//
	"computeDiffusion=false",	";Compute diffusion coeficient",
	"stepDiffuse=4",				";number of step jumps to save a diffusion measurement",
	"stepAvgDiffuse=200",			";number of diffusion measurements to average",
	"nBuffDiffuse=10",				";size of buffer to save a diffusion measurements",
	"nValDiffuse=250",				";number of values to save of diffusion measurement",
//
	"computeVelAcf=false",		";Compute velocity autocorrelation function",
	"computeBulkViscosity=false",		";Compute bulk viscosity",
//
	"computeTransport=false",	";Compute transport properties",
	"stepAcf=3",				";number of step jumps to save an Acf measurement",
	"stepAvgAcf=200",			";number of Acf measurements to average",
	"nBuffAcf=10",				";size of buffer to save an Acf measurements",
	"nValAcf=200",				";number of values to save of Acf measurement",
//
	"computeSTCorr=false",		";Compute space-time correlations",
	"stepCorr=5",				";number of step jumps to save a Corr measurement",
	"stepAvgCorr=500",			";number of Corr measurements to average",
	"nBuffCorr=80",				";size of buffer to save a Corr measurements",
	"nValCorr=1025",				";number of values to save of Corr measurement",
	"nFunCorr=4",				";number of values to save of Corr measurement",
//
	"lattCorr_kx=1",			";x-wavenumber for lattice correlation computation", ":lCkx",
	"lattCorr_ky=-1",			";y-wavenumber for lattice correlation computation", ":lCky",
#ifdef THREEDIM
	"lattCorr_kz=1",			";z-wavenumber for lattice correlation computation", ":lCkz",
#endif
//
// ... END OUTPUT ...
//	
    "Version=0.2",              ";Mar 2005-2011",
    NULL,
};

#endif /* ! _cmdline_defs_h */

