/*==============================================================================
	HEADER: protodefs.h			[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main), start_run, time_step,
					forcecalc, md_lj_tree_io
	External headers:
	Comments and notes:
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _proto_defs_h
#define _proto_defs_h
														// CHECK 2D --- OK!!!
void data_analysis(void);

#ifdef THREEDIM
void snap_anim_plplot_3d(string, string);
void snap_anim_plplot_trajectory_3d(string, string);
#endif
void snap_anim_plplot(string, string);
void snap_anim_plplot_trajectory(string, string);
void rdf_anim_plplot(string, string);
void vel_anim_plplot(string, string);
void general_fx_anim_plplot(string, string);
void rhoaxes_anim_plplot(string, string);
void nfrecaxes_anim_plplot(string, string);
void bodies_anim_plplot(string, string);
void bodies_anim_plplot_trajectory(string, string);

void EvalVelDist(string, string);
void EvalRdf(string, string);
void EvalRhoAxes(string, string);
void EvalNFrecAxes(string, string);

void thermo_avg(void);

void units_conversion(void);

/* ------------[	Input/Output Prototype definitions	 ]----------- */

void startoutput(void);
void startrun(string, string, string, string);

void snap_conversion(void);
void snap_less_nbody(void);

#ifdef THREEDIM
//void PrintSnap_Slab(string, string, int, char *);
void PrintSnap_Slab(char *, char *, int, char *, char mode[2]);
#endif

void Header_to_Global(void);
void Global_to_Header(void);

#endif	/* ! _proto_defs_h	*/
