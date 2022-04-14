/*==============================================================================
	HEADER: proto_defs.h			[model]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global routines
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main), start_run, time_step,
					forcecalc, md_lj_tree_io
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 23, 2007
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _proto_defs_h
#define _proto_defs_h

void InputData(void);
// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
void InputData_long(void);
//
void StartOutput(void);
void output(void);

void outputdata_body(string, int, int);

void StartRun(string, string, string, string);
void EndRun(void);

void testdata(void);


#endif	/* ! _prot_defs_h	*/