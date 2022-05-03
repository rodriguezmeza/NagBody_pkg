/* ==============================================================================
	HEADER: nagbody_proto.h			[NagBody]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
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

	Major revisions:
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _nagbody_proto_h
#define _nagbody_proto_h

#include "globaldefs.h"

void maketree(bodyptr, int, char *, global_data_tree *, 
				global_data_tree_bljforcecalc *);
void maketree_t(bodyptr, int, char *, global_data_tree *, 
				global_data_tree_tljforcecalc *);

// I/O prototype definitions: --------------------------------------------------
void inputdata(char *, char *, char *, int, int *, int *, realptr, bool *, 
			void *, char *, char *);
void outputdata(char *, char *, char *, int, int, real, void *, char *);
// -----------------------------------------------------------------------------

// (MD LJ Liquids) Several methods and type of forces --------------------------
void  ljforcecalc_barnes(bodyptr, int, global_data_tree_bljforcecalc *, 
							global_data_tree *);
void ljforcecalc_barnes2(bodyptr, int, global_data_tree_bljforcecalc *, 
							global_data_tree *);
void ljforcecalc_normal(bodyptr, int, global_data_tree_bljforcecalc *, 
							global_data_tree *);
void ind_ljforcecalc_normal(bodyptr, int, global_data_tree_bljforcecalc *, 
	global_data_tree *, bodyptr p, real *, real *);
void ljforcecalc_normal2(bodyptr, int, global_data_tree_bljforcecalc *, 
							global_data_tree *);
void ljforcecalc_normal3(bodyptr, int, global_data_tree_bljforcecalc *, 
							global_data_tree *);
void ljforcecalc_nblist(bodyptr, int, global_data_tree_bljforcecalc *, 
							global_data_tree *);
void ljforcecalc_direct(bodyptr, int, global_data_tree_bljforcecalc *);
void ind_ljforcecalc_direct(bodyptr, int, global_data_tree_bljforcecalc *,
	bodyptr p, real *, real *, short *);
void ljforcecalc_direct2(bodyptr, int, global_data_tree_bljforcecalc *);
void ljforcecalc_bazant(bodyptr, int, global_data_tree_bljforcecalc *);
//void ljforcecalc_sw(bodyptr, int, global_data_tree_bljforcecalc *);
void ljforcecalc_cellsmethod(bodyptr, int, global_data_tree_bljforcecalc *);
void ljforcecalc_cellsmethod3(bodyptr, int, global_data_tree_bljforcecalc *);
void ljpotcalc_cellsmethod(bodyptr, int, global_data_tree_bljforcecalc *);
// -----------------------------------------------------------------------------

// (MD TLJ Liquids) Several methods and type of forces --------------------------
void tljforcecalc_barnes(bodyptr, int, global_data_tree_tljforcecalc *, 
							global_data_tree *);
void tljforcecalc_barnes2(bodyptr, int, global_data_tree_tljforcecalc *, 
							global_data_tree *);
void tljforcecalc_normal(bodyptr, int, global_data_tree_tljforcecalc *, 
							global_data_tree *);
void ind_tljforcecalc_normal(bodyptr, int, global_data_tree_tljforcecalc *, 
	global_data_tree *, bodyptr p, real *, real *);
void tljforcecalc_normal2(bodyptr, int, global_data_tree_tljforcecalc *, 
							global_data_tree *);
void tljforcecalc_normal3(bodyptr, int, global_data_tree_tljforcecalc *, 
							global_data_tree *);
void tljforcecalc_nblist(bodyptr, int, global_data_tree_tljforcecalc *, 
							global_data_tree *);
void tljforcecalc_direct(bodyptr, int, global_data_tree_tljforcecalc *);
void ind_tljforcecalc_direct(bodyptr, int, global_data_tree_tljforcecalc *,
	bodyptr p, real *, real *, short *);
void tljforcecalc_direct2(bodyptr, int, global_data_tree_tljforcecalc *);
void tljforcecalc_cellsmethod(bodyptr, int, global_data_tree_tljforcecalc *);
void tljforcecalc_cellsmethod3(bodyptr, int, global_data_tree_tljforcecalc *);
void tljpotcalc_cellsmethod(bodyptr, int, global_data_tree_tljforcecalc *);
// -----------------------------------------------------------------------------

void PotentialParameters(int, real, real, real,
	real, real, real, real, real, real, global_data_tree_bljforcecalc *);
void PotentialParameters_t(int, real, real, real, real, real, real,
	real, real, real, real, real, real,
	real, real, real, real, real, real, global_data_tree_tljforcecalc *);
// -----------------------------------------------------------------------------


#endif	// ! _nagbody_proto_h
