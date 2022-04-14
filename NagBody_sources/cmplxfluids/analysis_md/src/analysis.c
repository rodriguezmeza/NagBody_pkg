/*==============================================================================
	MODULE: analysis.c				[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date:	May 2006
	Purpose:
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes:
	Info: Mario A. Rodriguez-Meza
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

#include "globaldefs.h"

#define SNAP_ANIM				1
#define SNAP_ANIMTRAJECTORY		3
#if (NDIM==3)
#define SNAP_ANIMTRAJECTORY3D	4
#define SNAP_ANIM3D				2
#define SLAB_ANIM				18
#define SLAB_ANIMTRAJECTORY		19
#endif
#define RDF_ANIM				5
#define VEL_ANIM				6
#define SNAP_ANIMGENERALFX		7
#define RHOAXES_ANIM			8
#define NFRECAXES_ANIM			9
#define BODIES_ANIM				10
#define BODIES_ANIMTRAJECTORY	11

#define SNAP_CONVERSION			200
#define SNAP_LESSNBODY			201

#define THERMO_AVG				300

#define UNITS_CONVERSION		500

local void datanaly_type_string_to_int(string, int *);

void data_analysis(void)								// CHECK 2D --- OK!!!
{
	string snapfilename, snapfilenametmp;
	int datanaly_type_int;

	snapfilename="snap-analysis.dat";
	snapfilenametmp="snap-analysis.tmp";

	datanaly_type_string_to_int(cmd.data_analysis_type, &datanaly_type_int);
	switch (datanaly_type_int) {
	case SNAP_ANIM: snap_anim_plplot(snapfilename, snapfilenametmp);	break;
	case SNAP_ANIMTRAJECTORY: snap_anim_plplot_trajectory(snapfilename, snapfilenametmp); break;
#if (NDIM==3)
	case SNAP_ANIMTRAJECTORY3D: snap_anim_plplot_trajectory_3d(snapfilename, snapfilenametmp); break;
	case SNAP_ANIM3D: snap_anim_plplot_3d(snapfilename, snapfilenametmp); break;
    case SLAB_ANIM: slab_anim_plplot(snapfilename, snapfilenametmp); break;
    case SLAB_ANIMTRAJECTORY: slab_anim_plplot_trajectory(snapfilename, snapfilenametmp); break;
#endif

    case RDF_ANIM: rdf_anim_plplot(snapfilename, snapfilenametmp); break;
	case VEL_ANIM: vel_anim_plplot(snapfilename, snapfilenametmp); break;
	case SNAP_ANIMGENERALFX: general_fx_anim_plplot(snapfilename, snapfilenametmp); break;
	case RHOAXES_ANIM: rhoaxes_anim_plplot(snapfilename, snapfilenametmp); break;
	case NFRECAXES_ANIM: nfrecaxes_anim_plplot(snapfilename, snapfilenametmp); break;
	case BODIES_ANIM: bodies_anim_plplot(snapfilename, snapfilenametmp); break;
	case BODIES_ANIMTRAJECTORY: bodies_anim_plplot_trajectory(snapfilename, snapfilenametmp); break;

	case SNAP_CONVERSION: snap_conversion(); break;
	case SNAP_LESSNBODY: snap_less_nbody(); break;

	case THERMO_AVG: thermo_avg(); break;

	case UNITS_CONVERSION: units_conversion(); break;

	default: printf("No data analysis type chosen.\n"); break;
	}
}

														// CHECK 2D --- OK!!!
local void datanaly_type_string_to_int(string datanaly_type_str,int *datanaly_type_int)
{
	if (strcmp(datanaly_type_str,"exit") == 0)
		*datanaly_type_int = 0;

	if (strcmp(datanaly_type_str,"snap-anim") == 0)	
		*datanaly_type_int = SNAP_ANIM;
	if (strcmp(datanaly_type_str,"snap-animtrajectory")	== 0)
		*datanaly_type_int = SNAP_ANIMTRAJECTORY;
#if (NDIM==3)
	if (strcmp(datanaly_type_str,"snap-anim3d")	== 0)
		*datanaly_type_int = SNAP_ANIM3D;
	if (strcmp(datanaly_type_str,"snap-animtrajectory3d") == 0)
		*datanaly_type_int = SNAP_ANIMTRAJECTORY3D;
	if (strcmp(datanaly_type_str,"slab-anim")	== 0) *datanaly_type_int = SLAB_ANIM;
	if (strcmp(datanaly_type_str,"slab-animtrajectory")	== 0) 
        *datanaly_type_int = SLAB_ANIMTRAJECTORY;
#endif

	if (strcmp(datanaly_type_str,"rdf-anim") == 0)
		*datanaly_type_int = RDF_ANIM;
	if (strcmp(datanaly_type_str,"vel-anim") == 0)
		*datanaly_type_int = VEL_ANIM;
	if (strcmp(datanaly_type_str,"general-fx-anim")	== 0)
		*datanaly_type_int = SNAP_ANIMGENERALFX;
	if (strcmp(datanaly_type_str,"rhoaxes-anim") == 0)
		*datanaly_type_int = RHOAXES_ANIM;
	if (strcmp(datanaly_type_str,"nfrecaxes-anim") == 0)
		*datanaly_type_int = NFRECAXES_ANIM;
	if (strcmp(datanaly_type_str,"bodies-anim")	== 0)
		*datanaly_type_int = BODIES_ANIM;
	if (strcmp(datanaly_type_str,"bodies-animtrajectory") == 0) 
		*datanaly_type_int = BODIES_ANIMTRAJECTORY;

	if (strcmp(datanaly_type_str,"snap-conversion")	== 0)	
		*datanaly_type_int = SNAP_CONVERSION;
	if (strcmp(datanaly_type_str,"snap-lessnbody")	== 0)
		*datanaly_type_int = SNAP_LESSNBODY;

	if (strcmp(datanaly_type_str,"thermo-avg")	== 0)
		*datanaly_type_int = THERMO_AVG;

	if (strcmp(datanaly_type_str,"units-conversion") == 0)	
		*datanaly_type_int = UNITS_CONVERSION;
}

#undef SNAP_ANIM
#undef SNAP_ANIMTRAJECTORY

#if (NDIM==3)
#undef SNAP_ANIM3D
#undef SNAP_ANIMTRAJECTORY3D
#undef SLAB_ANIM
#undef SLAB_ANIMTRAJECTORY
#endif

#undef RDF_ANIM
#undef VEL_ANIM
#undef SNAP_ANIMGENERALFX
#undef RHOAXES_ANIM
#undef NFRECAXES_ANIM
#undef BODIES_ANIM
#undef BODIES_ANIMTRAJECTORY

#undef SNAP_CONVERSION
#undef SNAP_LESSNBODY

#undef THERMO_AVG

#undef UNITS_CONVERSION
