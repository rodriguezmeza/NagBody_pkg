/*==============================================================================
	NAME: main.c					[nplot2d]
	Written by: Mario A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: Main routine - Plotting data in 2D
	Language: C
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revision:	May 2006; November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.
================================================================================

	Use: nplot2d -help
	Input: 	Command line parameters, Parameters file and/or datafile
	Output: ...
	Units:
	History:
	Comments and notes:
	References:
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#define global

#include "globaldefs.h"
#include "cmdlinedefs.h"

int main(int argc, string argv[])
{
    gd.cpuinit = cputime();

    InitParam(argv, defv);
    StartRun(argv[0], HEAD1, HEAD2, HEAD3);
    StartOutput();
	Plot2DDriver();
	EndRun();
	
    return 0;
}


