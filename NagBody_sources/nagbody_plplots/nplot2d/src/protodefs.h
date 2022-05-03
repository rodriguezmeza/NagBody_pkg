/*==============================================================================
	HEADER: proto_defs.h		[nplot2d]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: Definitions of global prototypes
	Language: C
	Use: '#include "proto_defs.h"
	Use in routines and functions:
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _proto_defs_h
#define _proto_defs_h

void StartOutput(void);
void StartRun(string, string, string, string);
void EndRun(void);
void Plot2DDriver(void);

void InputData(string, int, int, int *);
void InputData_3c(string, int, int, int, int *);
void InputData_4c(string, int, int, int, int, int *);
void TestData(int *);

#endif /* ! _proto_defs_h */
