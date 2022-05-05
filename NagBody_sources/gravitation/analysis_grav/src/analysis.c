/*==============================================================================
	MODULE: analysis.c				[analysis_grav]
	Written by: M.A. Rodriguez-Meza
	Starting date:	May 2006
	Purpose:
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
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
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

#define SNAP_ANIM				1
#define SNAP_ANIMTRAJECTORY		3
#define SNAP_ANIMGENERALFX		7
#define SLAB_ANIM				8
#define RHOTHETA_ANIM			20
#define VCR_ANIM				21
#define ACF_ANIM				22
#define BODIES_ANIM				10
#define BODIES_ANIMTRAJECTORY	11
#define LOCATEBODIESID			12

#define RHO_ANIM				50
//#define VC_ANIM					51

#define SNAP_CONVERSION			200
#define SNAP_LESSNBODY			201
#define GROUPS_CATALOG			202

local void datanaly_type_string_to_int(string, int *);

void data_analysis(void)
{
	string snapfilename, snapfilenametmp;
	int datanaly_type_int;

	snapfilename="snap-analysis.dat";
	snapfilenametmp="snap-analysis.tmp";

	datanaly_type_string_to_int(cmd.data_analysis_type, &datanaly_type_int);
	switch (datanaly_type_int) {
	case SNAP_ANIM: snap_anim_plplot(snapfilename, snapfilenametmp);	break;
	case SNAP_ANIMTRAJECTORY: snap_anim_plplot_trajectory(snapfilename, snapfilenametmp); break;
	case SNAP_ANIMGENERALFX: general_fx_anim_plplot(snapfilename, snapfilenametmp); break;
	case SLAB_ANIM: slab_anim_plplot(snapfilename, snapfilenametmp); break;
	case RHOTHETA_ANIM: rhotheta_anim_plplot(snapfilename, snapfilenametmp); break;
	case VCR_ANIM: vcr_anim_plplot(snapfilename, snapfilenametmp); break;
	case ACF_ANIM: acf(snapfilename, snapfilenametmp); break;
	case BODIES_ANIM: bodies_anim_plplot(snapfilename, snapfilenametmp); break;
	case BODIES_ANIMTRAJECTORY: bodies_anim_plplot_trajectory(snapfilename, snapfilenametmp); break;
	case LOCATEBODIESID: locate_bodiesid(snapfilename, snapfilenametmp); break;

	case RHO_ANIM: rho_anim_plplot(snapfilename, snapfilenametmp); break;
//	case VC_ANIM: vc_anim_plplot(snapfilename, snapfilenametmp); break;

	case SNAP_CONVERSION: snap_conversion(); break;
	case SNAP_LESSNBODY: snap_less_nbody(); break;
	case GROUPS_CATALOG: groups_catalog(); break;

	default: printf("No data analysis type chosen : %d\n",datanaly_type_int); break;
	}
}

local void datanaly_type_string_to_int(string datanaly_type_str,int *datanaly_type_int)
{
	if (strcmp(datanaly_type_str,"exit")		== 0)		*datanaly_type_int = 0;

	if (strcmp(datanaly_type_str,"snap-anim")	== 0)		*datanaly_type_int = 1;
	if (strcmp(datanaly_type_str,"snap-animtrajectory")	== 0) *datanaly_type_int = 3;
	if (strcmp(datanaly_type_str,"slab-anim")	== 0) *datanaly_type_int = SLAB_ANIM;
	if (strcmp(datanaly_type_str,"general-fx-anim")	== 0)	*datanaly_type_int = 7;
	if (strcmp(datanaly_type_str,"rhotheta-anim")	== 0)	*datanaly_type_int = 20;
	if (strcmp(datanaly_type_str,"vcr-anim")	== 0)		*datanaly_type_int = 21;
	if (strcmp(datanaly_type_str,"acf-anim")	== 0)		*datanaly_type_int = 22;
	if (strcmp(datanaly_type_str,"bodies-anim")	== 0)			*datanaly_type_int = 10;
	if (strcmp(datanaly_type_str,"bodies-animtrajectory") == 0) *datanaly_type_int = 11;
	if (strcmp(datanaly_type_str,"locate-bodiesid")	== 0)	*datanaly_type_int = 12;

	if (strcmp(datanaly_type_str,"rho-anim")	== 0) *datanaly_type_int = RHO_ANIM;
//	if (strcmp(datanaly_type_str,"vc-anim")	== 0) *datanaly_type_int = VC_ANIM;

	if (strcmp(datanaly_type_str,"snap-conversion")	== 0)	*datanaly_type_int = 200;
	if (strcmp(datanaly_type_str,"snap-lessnbody")	== 0) *datanaly_type_int = SNAP_LESSNBODY;
	if (strcmp(datanaly_type_str,"groups-catalog")	== 0) *datanaly_type_int = GROUPS_CATALOG;
}

#undef SNAP_ANIM
#undef SNAP_ANIMTRAJECTORY
#undef SNAP_ANIMGENERALFX
#undef SLAB_ANIM
#undef RHOTHETA_ANIM
#undef VCR_ANIM
#undef ACF_ANIM
#undef BODIES_ANIM
#undef BODIES_ANIMTRAJECTORY
#undef LOCATEBODIESID

#undef RHO_ANIM
//#undef VC_ANIM

#undef SNAP_CONVERSION
#undef SNAP_LESSNBODY
#undef GROUPS_CATALOG

