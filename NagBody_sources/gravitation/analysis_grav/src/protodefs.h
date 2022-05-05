/*==============================================================================
	HEADER: protodefs.h			[analysis_grav]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: analysis_grav (main), startrun,
					analysis_grav_io
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
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

#ifndef _proto_defs_h
#define _proto_defs_h

void data_analysis(void);

//void readin_pl_snap(string, int, int, int *);
//void readin_pl_snap_3d(string, int, int, int, int *);

//void snap_anim_plplot_3d(string, string);
void snap_anim_plplot(string, string);
void snap_anim_plplot_trajectory(string, string);
void slab_anim_plplot(string, string);
//void snap_anim_plplot_trajectory_3d(string, string);
//void rdf_anim_plplot(string, string);
//void vel_anim_plplot(string, string);
void general_fx_anim_plplot(string, string);
//void rhoaxes_anim_plplot(string, string);
void rhotheta_anim_plplot(string, string);
void vcr_anim_plplot(string, string);
void acf(string, string);
//void nfrecaxes_anim_plplot(string, string);
void bodies_anim_plplot(string, string);
void bodies_anim_plplot_trajectory(string, string);
void locate_bodiesid(string, string);

void PrintSnap_Spherical(string, string, int, char *);
void PrintBodiesSets_Spherical(string, string, int, int, int *, int *, real, char *);

void PrintSnap_Slab(string, string, int, char *);

//void snap_anim_gnuplot(string, string);
//void vel_anim_gnuplot(string, string);
//void rdf_anim_gnuplot(string, string);

//void EvalVelDist(string, string);
//void EvalRdf(void);
//void EvalRdf(string, string);
//void EvalRhoAxes(string, string);
void EvalRhoTheta(string, string);
void EvalVcR(string, string);
void EvalAcf(string, string);
void LocateBodiesID(string, string);
//void EvalNFrecAxes(string);

void EvalRhoDist(string, string);
void rho_anim_plplot(string, string);
//void vc_anim_plplot(string, string);


//void PrintSnap(FILE *);
//void readin_snap(int, bool *);

//void thermo_avg(void);


/* ------------[	Input/Output Prototype definitions	 ]----------- */

void startoutput(void);
void startrun(void);

//void inputdata(int, bool *);
//void outputdata(int);
void snap_conversion(void);
void snap_less_nbody(void);

void groups_catalog(void);
void GroupsCatalog(int);

//void Header_to_Global(void);
//void Global_to_Header(void);

//void maketree_grav(bodyptr, int);
//void normal_gravcalc(bodyptr, int);
//void ind_normal_gravcalc(bodyptr, int, bodyptr);

#endif	/* ! _proto_defs_h	*/
