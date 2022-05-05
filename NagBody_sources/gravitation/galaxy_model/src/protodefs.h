/* ==============================================================================
	HEADER: protodefs.h			[galaxy_models]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: Definitions of global prototypes
	Language: C
	Use: '#include "proto_defs.h"
	Use in routines and functions:
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: June 11, 2008;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved.
=================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
===============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void StartOutput(void);
void StartRun(string, string, string, string);
void EndRun(void);

void FreeArrays(void);

void InitDisk(void);
void InitBulge(void);
void InitHalo(void);
void DiskVel(void);
void DiskStat(void);
void StackMod(void);
void cmtv(void);
void AddSat(void);

void HaloVel(void);
void BulgeVel(void);

// Bulge protodefs from galmod:
//void inbmass(void);
//void setbulge(void);
//void cmbulge(void);


void obsigma(real, real, realptr, realptr, realptr, real,
					real, real, int, real, real);

real cfunc(real, real, real);

#endif /* ! _protodefs_h */
