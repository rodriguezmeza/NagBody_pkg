/*==============================================================================
   NAME: main.c					[gbsph]
	Written by: Mario Alberto Rodriguez-Meza
	Starting date: October 11, 1999
	Purpose: Simulation of general body systems and smoothed particle
			hydrodynamics (Hierarquical tree code)
	Language: C
	Info:	M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revision:	May 30, 2000; November 20, 2001; July 24, 2007;
	Copyright: (c) 2005-2008 Mar.  All Rights Reserved
================================================================================

	Use: gbsph -help
	Input: Command line parameters, Parameters file and/or icfile
	Output: energy.dat, gbsph.log, parameters-usedvalues, snapshot data, ...
	Units: G=1, ...
	History: Version 0.1 Final (SPHC) May 30, 2000;
		Update (GBSPHC) November 20, 2001.
		Version 0.2: SPHC, GBSPH, PGBSPH Integration (in process).
			Transparent use of parameter file and command line
			parameters, support of several SNAP formats,
			several types of forces (i.e., Yukawa). Generation of
			several kind of models.
	Comments and notes: 
	References: Hernquist & katz, Barnes & Hut, Springel et al.
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/vectmath.h"
//#include "../../../General_libs/general/getparam.h"
#define global
#include "globaldefs.h"
#include "cmdline_defs.h"
#include "protodefs.h"

int main(int argc, string argv[])
{
	gd.cpuinit = cputime();
    InitParam(argv, defv);
    StartRun(argv[0], HEAD1, HEAD2, HEAD3);
	MainLoop();
	EndRun();
    return 0;
}

