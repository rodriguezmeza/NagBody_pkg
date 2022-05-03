/*==============================================================================
	HEADER: global_defs.h		[nplot2d]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "global_defs.h"
	Use in routines and functions:
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: marioalberto.rodriguez@inin.gob.mx
        http://www.inin.gob.mx/

	Major revisions: November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#ifndef _global_defs_h
#define _global_defs_h

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include "plevent.h"
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#ifndef NOGNU
#include "./general_libs/general/stdinc.h"
#include "./general_libs/math/vectdefs.h"
#include "./general_libs/NagBody/nagbody.h"
#include "./general_libs/io/inout.h"
#include "./general_libs/math/mathfns.h"
#include "./general_libs/general/getparam.h"
#include "./general_libs/general/lic.h"
#include "./general_libs/visual/pldefs.h"
#else
#include "stdinc.h"
#include "vectdefs.h"
#include "nagbody.h"
#include "inout.h"
#include "mathfns.h"
#include "getparam.h"
#include "lic.h"
#include "pldefs.h"
#endif

#include "protodefs.h"

// -----Block of command line definitions---------------------------------------

typedef struct {										// CHECK 2D --- OK!!!

	string paramfile;
	string inputfile;

	string witherrorbars;
	int errorbarstype;

	string graphicsarray;
	int plottype;
	string xrange;
	string yrange;
	string usingcolumns;
	string xlabel;
	string ylabel;
	string plotlabel;
	int labelcolor;

	string legends;
	int legendspos;
	real legendsfontsize;
	int legendsfontweight;
	int fontset;
	int fontchr;
	real labelfontsize;
	int labelfontweight;

	int nlxdigmax;
	int nlydigmax;
	real nlsize;
	int nlcolor;

	string axesorigin;

	bool xaxis;
	bool yaxis;
	int axestype;
	int axeswidth;
	int axescolor;
	string plotjoined;
	string linetype;
	int linewidth;
	int linecolor;
	string withdots;
	string withsymbols;
	string symboltype;
	int symbolcolor;
	int symbolweight;
	real symbolsize;

	string text1;
	string text1side;
	real text1disp;
	real text1pos;
	real text1just;
	real text1size;
	int text1weight;
	int text1color;

	string text2;
	string text2side;
	real text2disp;
	real text2pos;
	real text2just;
	real text2size;
	int text2weight;
	int text2color;

	string text3;
	string text3side;
	real text3disp;
	real text3pos;
	real text3just;
	real text3size;
	int text3weight;
	int text3color;

	string text4;
	real text4x;
	real text4y;
	real text4dx;
	real text4dy;
	real text4just;
	real text4size;
	int text4weight;
	int text4color;

	bool locatemode;

	bool frame;
	string framestyle;	// More detailed control of axis, frame, ticks, ...
	bool gridlines;
	string epilog;


//---------MENU COMMAND LINE PLPLOT OPTIONS------------------------------------
// Option set to false == use default, they are bool

	bool pl_showall;
	bool pl_h;
	bool pl_v;
	bool pl_verbose;
	bool pl_debug;
	bool pl_hack;
	string pl_dev;
	string pl_o;
	string pl_display;
	string pl_px;
	string pl_py;
	string pl_geometry;
	string pl_wplt;
	string pl_mar;
	string pl_a;
	string pl_jx;
	string pl_jy;
	string pl_ori;
	bool pl_freeaspect;
	bool pl_portrait;
	string pl_width;
	string pl_bg;
	string pl_ncol0;
	string pl_ncol1;
	bool pl_fam;
	string pl_fsiz;
	string pl_fbeg;
	string pl_finc;
	string pl_fflen;
	bool pl_nopixmap;
	bool pl_db;
	bool pl_np;
	string pl_bufmax;
//	string pl_server_name;
//	string pl_plserver;
//	string pl_plwindow;
//	string pl_tcl_cmd;
//	string pl_auto_path;
//	string pl_tk_file;
	string pl_dpi;
	string pl_compression;
	string pl_drvopt;

} cmdline_data, *cmdline_data_ptr;


//------------------------------------------------------------------------------

#define MAXLINES	1000

typedef struct {										// CHECK 2D --- OK!!!

	int nfiles;
	char *filenames[40];

	int nwitherrorbars;
	bool errorbars[MAXLINES];

	int nplotjoined;
	bool plotjoined[MAXLINES];

	int nwithsymbols;
	bool withsymbols[MAXLINES];

	int nwithdots;
	bool withdots[MAXLINES];

	int linetype[MAXLINES];
	int symboltype[MAXLINES];

	int nlegends;
	char *legendnames[40];

	int graphicsarray_nrows;
	int graphicsarray_ncols;

// ----- Block of definitions for xrange and yrange ----------------------------
	bool x_autoscale;
	bool y_autoscale;
	real xmin;
	real xmax;
	real ymin;
	real ymax;

// ----- Block of definitions for usingcolumns ---------------------------------
	int *vcol1;
	int *vcol2;
	int *vcol3;
	int *vcol4;

// ----- Block of definitions for framestyle -----------------------------------
	string frame_xopt;
	string frame_yopt;
	real frame_xtick;
	int frame_nxsub;
	real frame_ytick;
	int frame_nysub;

// ----- Block of definitions for headlines -----------------------------------
	string headline0;
	string headline1;
	string headline2;
	string headline3;
	string comment;			
	FILE *outlog;

// ----- Other definitions -----------------------------------------------------
	real xorigin;
	real yorigin;

	int npoint;

	long seed;										// SE ESTA USANDO? Si... testdata ...

	real *xval;										// ESTA DEFINIDO EN nagbody_struct
	real *yval;										// ESTA DEFINIDO EN nagbody_struct
	real *yminval;
	real *ymaxval;

	real cpuinit;

} global_data, *global_data_ptr;

// DEFINICION GLOBAL DE ESTRUCTURAS DE DATOS...
global global_data gd;
global cmdline_data cmd;

#endif // ! _global_defs_h

