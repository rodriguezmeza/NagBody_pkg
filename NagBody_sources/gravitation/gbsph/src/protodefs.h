/* =============================================================================
	HEADER: protodefs.h				[gbsph]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: forcecalc.c, gbsphio.c, main.c, startrun.c,
									timestep.c
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 24, 2007; October 04, 2007;
	Copyright: (c) 1999-2008 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void forcecalc(bodyptr, int);  
void forcecalc_collective_motion_without_leader(bodyptr, int);  
void external_forcecalc(bodyptr, int);

void inputdata(void);                   
void forcereport(void);                 
void output(void);                      
void restorestate(string);

void force_setkernel(void);

void maketree(bodyptr, int);            

void gravcalc(void);                    

void MainLoop(void);
void StartRun(string, string, string, string);
void StartOutput(void);
void checkstop(void);
void EndRun(void);

#endif	// ! _protodefs_h

