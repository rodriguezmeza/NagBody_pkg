/*==============================================================================
 NAME: main.c						[templateBessel]
 Written by: Mario Alberto Rodriguez-Meza
 Starting date: January, 2005
 Purpose: Simulation (N-body) of a selfgravitating N-body system.
 N*N force computation is done.
 Language: C
 Info: M.A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.astro.inin.mx/mar
 
 Major revision:	July 2007; November 2008; January 2018
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 
 Use: templateBessel -help
 Input: 	Command line parameters, Parameters file ...
 History:
 Comments and notes: ...
 References:
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

#include "../../../General_libs/general/stdinc.h"
#include "../../../General_libs/general/getparam.h"
#include "../../../General_libs/math/mathfns.h"

#include "protodefs.h"

string defv[] = {	";Test bessel functions",
"xvalue=1.65",		";Argument to Bessel functions", ":x",
"Version=0.2",		";Mar 2005-2018",
NULL,
};

int main(int argc, string argv[])
{
    double x;
	
    InitParam(argv, defv);
    x = GetdParam("x");
	fprintf(stdout,"# x\t\t BesselI(0,x)\t BesselI(1,x)\t BesselK(0,x)\t BesselK(1,x)\n");
	fprintf(stdout,"%f\t%f\t%f\t%f\t%f\n",x,bessi0(x),bessi1(x),bessk0(x),bessk1(x));
	return (0);
}


