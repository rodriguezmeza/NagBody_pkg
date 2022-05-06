/* ==============================================================================
!	HEADER: proto_defs.h														!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date: February 2005												!
!	Purpose: Definitions of prototypes											!
!	Language: C																	!
!	Use: '#include "...."														!
!	Use in routines and functions: md_lj_tree (main), startrun, timestep,		!
!					forcecalc, md_lj_tree_io									!
!	External headers: None														!
!	Comments and notes:															!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.							!
!		e-mail: marioalberto.rodriguez@inin.gob.mx
!		http://www.astro.inin.mx/mar											!
!																				!
!	Major revisions:															!
!	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.							!
!===============================================================================
!	Legal matters:																!
!	The author does not warrant that the program and routines it contains		!
!	listed below are free from error or suitable for particular applications,	!
!	and he disclaims all liability from any consequences arising from their		!
!	use.																		!
!==============================================================================*/

#ifndef _proto_defs_h
#define _proto_defs_h

void inputdata(void);       
void startoutput(void);                 
void output(void);
//void savestate(string);                 
void restorestate(string);              

void tree_ljforce(void);               
void stepsystem(void);

void startrun(void);                    
void checkstop(void);
void code_endrun(void);

void maketree(bodyptr, int);            

// (MD LJ Liquids) Several methods and type of forces
//void ljforcecalc_barnes(void);				// Opcion para solo correr barnes
void ljforcecalc_barnes(bodyptr, int);
void ljforcecalc_barnes2(bodyptr, int);
void ljforcecalc_normal(bodyptr, int);
void ljforcecalc_normal2(bodyptr, int);
void ljforcecalc_nblist(bodyptr, int);
void ljforcecalc_direct(bodyptr, int);
void ljforcecalc_direct2(bodyptr, int);
void ljforcecalc_cellsmethod(bodyptr, int);


#endif /* ! _proto_defs_h */
