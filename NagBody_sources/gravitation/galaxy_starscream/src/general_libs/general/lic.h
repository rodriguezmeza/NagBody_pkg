/*==============================================================================
	HEADER: lic.h					[General_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: definitions and protodefinitios to license the code
	Language: C
	Use: '#include "lic.h"'
	Use in routines and functions:
	External headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: January 2007;
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#ifndef _lic_h
#define _lic_h
#include "stdinc.h"

typedef struct {
    int id;
    char date[12];
    char licid[50];
    char machineid[50];
    char passwd[50];
    char pathusername[200];
    char firstname1[15];
    char firstname2[15];
    char name1[15];
    char name2[15];
	char institution[200];
	char email[100];
} license;

typedef struct {
	int year;
	int month;
	int day;
} date;

void licDriver(char *);

#endif	/* ! _lic_h */
