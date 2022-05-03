/*==============================================================================
	HEADER: cmdlinedefs.h			[nplot2d]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "cmdlinedefs.h"
	Use in routines and functions: (main)
	External headers: stdinc.h
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

//#ifndef NOGNU
//#include "../../../General_libs/general/stdinc.h"
//#else
//#include "stdinc.h"
//#endif

#define HEAD1	"NagBody"
#define HEAD2	"plotting 2d graphics code..."
#define HEAD3	"Based on plplot library"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			";Parameter input file. Overwrite what follows",
	"graphicsarray=1x1",	";Graphics array style (nrows x ncols)", ":ga",
    "inputfile=",			";File to process", ":in",
    "outputfile=",			";Output file to save the plot", ":o",
    "epilog=",				";More options to plot",
    "labelcolor=2",			";Label color",					// Color for xlabel, ylabel and plotlabel
	"plottype=0",			";Plot type - normal (0), normal-log (1), log-normal (2), log-log (3)",
	"usingcolumns=1:2",		";Columns to plot", ":uc",
    "xlabel=",				";X-axis label", ":xl",
	"xrange=autoscale",		";Range in x-axis", ":xr",
    "ylabel=",				";Y-axis label", ":yl",
	"yrange=autoscale",		";Range in y-axis", ":yr",
//    "xlabelcolor=2",		";X-axis label color",				// Color for xlabel, ylabel and plotlabel
//    "ylabelcolor=2",		";Y-axis label color",				// are controled by the same instruction...
	"legends=",				";Legends for data",
	"legendspos=0",			";Position of legends: top-left (0), top-right (1), bottom-left(2), botton-right (3)",
    "legendsfontsize=1",		";Legends font scale size",
    "legendsfontweight=1",	";Legends font weight",
    "fontset=0",			";Font set switch",
    "fontchr=1",			";Set character font (1-normal, 2-Roman, 3-Italic, 4-Script)",
//    "xlabelfont=0",			";X-axis label font",
//    "ylabelfont=0",			";Y-axis label font",
    "labelfontsize=1",		";Label font scale size",
    "labelfontweight=1",	";Label font weight",
//    "xlabelfontsize=1",		";X-axis label font scale size",
//    "ylabelfontsize=1",		";Y-axis label font scale size",
    "xaxis=false",			";Drawing axes",
    "yaxis=false",			";Drawing axes",
	"background=",			";Background color (000000=black, FFFFFF=white)", ":bg",
	"aspectratio=",			";Page aspect ratio (default: same as output device)", ":a",
	"frame=true",			";Display frame",
	"framestyle=",			";Frame detailed control style",			// Is controled by axescolor, axeswidth, axestype
//
//
	"axesorigin=0,0",		";Specifies where any axes drawn should cross",  // ACTIVAR ESTA OPCION....
//
// Usar en combinacion con la opcion framestyle ...
//
//
	"nlxdigmax=0",			";Maximum field width for numeric labels on x-axis",
	"nlydigmax=0",			";Maximum field width for numeric labels on y-axis",
	"nlsize=1",				";Numeric label font size",
	"nlcolor=1",				";Numeric label font color",
//
	"axescolor=1",			";Axes color",
	"axeswidth=1",			";Axes width",
	"axestype=1",			";Axes type",
	"gridlines=false",		";Grid lines at the major tick interval (both axes)",
/*
    "axesorigin=automatic",	";Origin's position",		// Axes origin is (0,0) always ... can be changed?
    "prolog=",				";More options to render before the main plot",		// Not active by now...
*/
	"plotlabel=",			";Plot label",
//	"plotlabelcolor=2",		";Plot label color",				// Check it with xlabel and ylabel
//	"plotlabelfont=0",		";Plot label font",
//	"plotlabelfontsize=1",	";Plot label font scale size",
	"plotjoined=true",		";Line or symbol", ":pj",
	"linetype=1",			";Line type",
	"linewidth=1",			";Line width",
	"linecolor=3",			";Line color",
	"withdots=false",		";Plot with dots (set plotjoined to true)", ":wd",
	"withsymbols=false",	";Plot with symbols", ":ws",
	"symboltype=1",			";Symbol type",
	"symbolcolor=4",		";Symbol color",
	"symbolweight=1",		";Symbol weight",
	"symbolsize=1",			";Symbol size",
//
	"witherrorbars=false",		";Plot error bars", ":web",
	"errorbarstype=0",		";Error bars type", ":errt",
//
	"text1=",				";More text to draw",
	"text1side=b",			";More text to draw - specify edge and direction (bottom and vertical, ...)",
	"text1disp=-5",			";More text to draw - position measured outwards from the specifed edge",
	"text1pos=0.3",			";More text to draw - position along the specified edge",
	"text1just=0",			";More text to draw - justification",
	"text1size=1",			";More text to draw - size",
	"text1weight=1",		";More text to draw - weight",
	"text1color=2",			";More text to draw - color",

	"text2=",				";More text to draw",
	"text2side=b",			";More text to draw - specify edge and direction (bottom and vertical, ...)",
	"text2disp=-3",			";More text to draw - position measured outwards from the specifed edge",
	"text2pos=0.3",			";More text to draw - position along the specified edge",
	"text2just=0",			";More text to draw - justification",
	"text2size=1",			";More text to draw - size",
	"text2weight=1",		";More text to draw - weight",
	"text2color=2",			";More text to draw - color",

	"text3=",				";More text to draw",
	"text3side=b",			";More text to draw - specify edge and direction (bottom and vertical, ...)",
	"text3disp=-1",			";More text to draw - position measured outwards from the specifed edge",
	"text3pos=0.3",			";More text to draw - position along the specified edge",
	"text3just=0",			";More text to draw - justification",
	"text3size=1",			";More text to draw - size",
	"text3weight=1",		";More text to draw - weight",
	"text3color=2",			";More text to draw - color",

	"text4=",				";More text to draw",
	"text4x=0",				";More text to draw - specify x position (world coordinates...)",
	"text4y=0",				";More text to draw - specify y position (world coordinates...)",
	"text4dx=1",			";More text to draw - specify dx direction (world coordinates...)",
	"text4dy=0",			";More text to draw - specify dy direction (world coordinates...)",
	"text4just=0",			";More text to draw - justification",
	"text4size=1",			";More text to draw - size",
	"text4weight=1",		";More text to draw - weight",
	"text4color=2",			";More text to draw - color",
//
	"locatemode=false",		";Locate mode",
//
//---------MENU COMMAND LINE PLPLOT OPTIONS------------------------------------
// Option set to false == use default
	"showall=false",		";Turns on invisible options",
	"h=false",				";Print out this message",
	"v=false",				";Print out the PLplot library version number",
	"verbose=false",		";Be more verbose than usual",
	"debug=false",			";Print debugging info (implies -verbose)",
	"hack=false",			";Enable driver-specific hack(s)",
	"dev=xwin",				";Output device name (try NULL string to see a list of devices)",
	"display=",				";X server to contact",
//	"px=1",					";Plots per page in x",
//	"py=1",					";Plots per page in y",
	"geometry=640x480",			";Window size, in pixels", ":geo",
	"wplt=0,0,1,1",			";Relative coordinates [0-1] of window into plot [xl,yl,xr,yr]",
	"mar=0",				";Margin space in relative coordinates (0 to 0.5, def 0)",
	"jx=0",					";Page justification in x (-0.5 to 0.5, def 0)",
	"jy=0",					";Page justification in y (-0.5 to 0.5, def 0)",
	"ori=0",				";Plot orientation (0,1,2,3=landscape,portrait,seascape,upside-down)",
	"freeaspect=false",		";Allow aspect ratio to adjust to orientation swaps",
	"portrait=false",		";Sets portrait mode (both orientation and aspect ratio)",
	"width=",				";Sets pen width (0 <= width)",
	"ncol0=",				";Number of colors to allocate in cmap 0 (upper bound)",
	"ncol1=",				";Number of colors to allocate in cmap 1 (upper bound)",
	"fam=false",			";Create a family of output files",
	"fsiz=",				";Output family file size (e.g. -fsiz 0.5G, def MB)",
	"fbeg=",				";First family member number on output",
	"finc=",				";Increment between family members",
	"fflen=",				";Family member number minimum field width",
	"nopixmap=false",		";Don't use pixmaps in X-based drivers",
	"db=false",				";Double buffer X window output",
	"np=false",				";No pause between pages",
	"bufmax=",				";bytes sent before flushing output",
//	"server_name=",			";Main window name of PLplot server (tk driver)",
//	"plserver=",			";Invoked name of PLplot server (tk driver)",
//	"plwindow=",			";Name of PLplot container window (tk driver)",
//	"tcl_cmd=",				";Depreciated - use -drvopt tcl_cmd= instead",
//	"auto_path=",			";Additional directory(s) to autoload (tk driver)",
//	"tk_file=",				";file for plserver (tk driver)",
	"dpi=",					";Resolution, in dots per inch (e.g. -dpi 360x360)",
	"compression=",			";Sets compression level in supporting devices",
	"drvopt=",				";Driver specific options",
//-----------------------------------------------------------------------------
    "Version=0.2",			";Mar 2005-2011",
    NULL,
};

#endif // ! _cmdline_defs_h
