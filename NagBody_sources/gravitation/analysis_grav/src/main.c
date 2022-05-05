/*==============================================================================
   NAME: main.c						[analysis_grav]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 1, 2006
	Purpose: Analysis of data coming from molecular dynamics codes (md_lj_tree,
		md_blj_tree, md_lj_n2)
	Language: C
	Info:	M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revision:
	Copyright: (c) 2006-2011 Mar.  All Rights Reserved
================================================================================

	Use: analysis_grav -help
	Input:
	Output:
	Units:
	History:
	Version 1.2 Final
	Comments and notes:
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#define global
#include "globaldefs.h"
#include "cmdline_defs.h"

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

