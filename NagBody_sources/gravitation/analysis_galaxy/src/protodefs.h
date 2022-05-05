/*==============================================================================
	HEADER: protodefs.h			[analysis_galaxy]
	Written by: M.A. Rodriguez-Meza
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
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void data_analysis(void);

void snap_anim_plplot(string, string);
void snap_anim_plplot_trajectory(string, string);
void general_fx_anim_plplot(string, string);
void rhotheta_anim_plplot(string, string);
void vcr_anim_plplot(string, string);
void acf(string, string);
void bodies_anim_plplot(string, string);
void bodies_anim_plplot_trajectory(string, string);
void locate_bodiesid(string, string);

void EvalRhoTheta(string, string);
void EvalVcR(string, string);
void EvalAcf(string, string);
void LocateBodiesID(string, string);

/* ------------[	Input/Output Prototype definitions	 ]----------- */

void startoutput(void);
void startrun(void);

void snap_conversion(void);
void snap_less_nbody(void);
void two_bdh_galaxies(string, string);

#endif	/* ! _protodefs_h	*/
