/*===============================================================================
!   NAME: md_lj_tree															!
!	Written by: Mario Alberto Rodriguez-Meza.									!
!	Starting date: January, 2005.												!
!	Purpose: Simulation (MD) of Lennard-Jones gas dynamics.						!
!		Hierarchical tree force computation is done.							!
!	Language: C																	!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.							!
!		e-mail: marioalberto.rodriguez@inin.gob.mx
!		http://www.astro.inin.mx/mar											!
!																				!
!	Major revision:	March 2006.													!
!	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.							!
!===============================================================================
!																				!
!	Use: ./md_lj_tree -help														!
!	Input: 	Command line parameters, Parameters file and/or icfile				!
!	Output: thermo.dat rdf.dat vel.dat snap.dat ...								!
!	Units:	eps=sigma=mass(p)=kB=1												!
!	History:																	!
!	Comments and notes: ...														!
!	References:	Barnes Treecode, NEMO project, Gadget, Rapaport's book			!
!===============================================================================
!	Legal matters:																!
!	The author does not warrant that the program and routines it contains		!
!	listed below are free from error or suitable for particular applications,	!
!	and he disclaims all liability from any consequences arising from their		!
!	use.																		!
!==============================================================================*/

#include "stdinc.h"
#include "vectmath.h"
#include "getparam.h"
#define global
#include "global_defs.h"
#include "cmdline_defs.h"
#include "proto_defs.h"

int main(int argc, string argv[])
{
    cpuinit = cputime();					// Set starting cpu time

    InitParam(argv, defv);					// Initialize parameter table
    headline0 = argv[0]; headline1 = HEAD1; // Initialize headers
    headline2 = HEAD2; headline3 = HEAD3;
    startrun();								// Read parameters and start
    startoutput();							// run and output
    if (nstep == 0) {                           
        tree_ljforce();	 					// Compute force
        output();							// Write output
    }
    if (dtime != 0.0)                            
        while (tstop - tnow > 0.01*dtime) {	// Main evolution loop
            stepsystem();                  	// Single advance step
            output();						// Write output
			checkstop();					// Check the run stop status
			if (stopflag) break;
        }
	code_endrun();							// Save restart file, close files
    return 0;
}

