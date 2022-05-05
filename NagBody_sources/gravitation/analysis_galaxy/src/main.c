/*==============================================================================
   NAME: main.c						[analysis_galaxy]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 1, 2006
	Purpose: Analysis of data coming from galactic dynamics codes 
			(galaxy_models, gbsph, pnbody, ...)
	Language: C
	Info:	M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revision:
	Copyright: (c) 2007-2011 Mar.  All Rights Reserved
==============================================================================

	Use: analysis_galaxy -help
	Input:
	Output:
	Units:
	History:
     		Version 0.2 Final


	Comments and notes:
-----------------------------------------------------------------------------*/

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/general/getparam.h"

#define global
#include "globaldefs.h"

#include "cmdline_defs.h"
#include "protodefs.h"

int main(int argc, string argv[])
{
	gd.cpuinit = cputime();

    InitParam(argv, defv);
    gd.headline0 = argv[0]; gd.headline1 = HEAD1;
    gd.headline2 = HEAD2; gd.headline3 = HEAD3;
    startrun(); 
    startoutput();
    data_analysis();
	code_endrun(gd.outlog, gd.cpuinit);
}

