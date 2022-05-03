/*==============================================================================
	HEADER: protodefs.h			[md_blj]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of prototypes
	Language: C
	Use: '#include "...."
	Use in routines and functions: md_lj_tree (main), startrun, timestep,
					forcecalc, md_lj_tree_io
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h
														// CHECK 2D --- OK!!!
void startoutput(void);
void output(void);
void restorestate(string);

void MainLoop(void);
void tree_ljforce(void);
void stepsystem(void);
void AdjustTemp(real);
void Diagnose();

void StartRun(string, string, string, string);
void StartOutput(void);
void checkstop(void);
void EndRun(void);
void testdata(void);

#endif /* ! _protodefs_h */
