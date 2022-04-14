/* ==============================================================================
	HEADER: cmdline_defs.h			[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: analysis_md (main)
	External headers: None
	Comments and notes:
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
===============================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code for data analysis of a binary N-body liquid simulation"
#define HEAD3	"Several schemes of force calculation"

														// CHECK 2D --- OK!!!
string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			";Parameter input file. Overwrite what follows",
//
//	"density=",				";Density of the liquid (as was used in snap computation)", ":d",
//	"temperature=",			";Temperature of the liquid (as was used in snap computation)", ":t",
//    "nbodyprop=1/1",		";Number of bodies (as was used in snap computation) (nbody1/nbody2)",
//    "massprop=1/1",			";Masses of the spieces (as was used in snap computation) (mass1/mass2)",
//#ifdef THREEDIM
//	"LxLyprop=1/1",			";Base sides of the parallelepiped (as was used in snap computation) (Lx/Ly)",
//#else
//	"Lx=1",					";x-side of base rectangle (as was used in snap computation)",
//#endif
//	"eps11=",				";Spieces-1 eps (as was used in snap computation)",
//	"eps12=",				";Cross eps (as was used in snap computation)",
//	"eps22=",				";Spieces-2 eps (as was used in snap computation)",
//	"sigma11=",				";Spieces-1 sigma (as was used in snap computation)",
//	"sigma12=",				";Cross sigma (as was used in snap computation)",
//	"sigma22=",				";Spieces-2 sigma (as was used in snap computation)",
//	"Rcut11=",				";Spieces-1 Cut radius (as was used in snap computation)",
//	"Rcut12=",				";Cross Cut radius (as was used in snap computation)",
//	"Rcut22=",				";Spieces-2 Cut radius (as was used in snap computation)",
//	"stepEquil=100",			";step to begin equilibrium computations",
//
	"stepAvg=1",				";number of histograms to average (RhoAxes, NFrecAxes, Vel, RDF)",
	"sizeHist=50",				";array size for histogram (RhoAxes, NFrecAxes, Vel, RDF)",
	"rangeVal=4.0",				";range of values for histogram (RhoAxes, NFrecAxes, Vel, RDF)",
//
//	"stepAvgRhoAxes=1",			";number of rho histograms to average",
//	"sizeHistRhoAxes=50",		";array size for rho histogram",
//
//	"stepAvgNFrecAxes=1",		";number of nfrec_axes histograms to average",
//	"sizeHistNFrecAxes=50",		";array size for nfrec_axes histogram",
//
//	"stepAvgVel=1",				";number of block of steps to average velocity",
//	"sizeHistVel=50",			";array size for velocity histogram",
//	"rangeVel=4.0",				";range of velocities for histogram",
//
//	"stepAvgRdf=1",			";number of block of steps to average Rdf",
//	"sizeHistRdf=200",			";array size for Rdf histogram",
//	"rangeRdf=4.0",				";range of Rdf for histogram",
//	
//    "nbody=512",				";Number of bodies for test run (must be the same as in snaps)",
    "in=",						";Input N-Body snap base name (snaps/snap%04d)",
    "infmt=snap-blj-ascii",		";N-Body snap file format",
    "out=",						";Output base name file of N-body frames or output file name", ":o",
    "outfmt=",					";Output file format of N-body frames (can be the extension)", ":ofmt",
	"basedir=.",				";Base directory for output files", ":bdir",
	"isnap=0",					";Initial snap number",
	"fsnap=10000",				";Final snap number",
//
	"unitsset=3.405e-10:6.6327e-26:1.64339e-21",	";Unit set: unit-length:unit-mass:unit-energy (MKS)",
//
    "options=",                 ";Various control options (see man pages for options)", ":opt",
    "analysis_type=snap-anim",	";Data analysis type to be done (see man pages for options)", ":at",
//
	"reductionFac=1",			";Reduction factor for analysis_type=snap-lessnbody",
	"bodiesID=",				";Body IDs to use with bodyies-anim (b1,b2,...)", ":bid",
	"xrange=autoscale",			";Range in x-axis (xmin:xmax)", ":xr",
	"yrange=autoscale",			";Range in y-axis (ymin:ymax)", ":yr",
#ifdef THREEDIM
	"zrange=autoscale",			";Range in z-axis (zmin:zmax)", ":zr",
//
	"xmin=0",			";Xmin for slab-anim (Choose this value acording snap)",
	"xmax=1",			";Xmax for slab-anim (Choose this value acording snap)",
	"ymin=0",			";Ymin for slab-anim (Choose this value acording snap)",
	"ymax=1",			";Ymax for slab-anim (Choose this value acording snap)",
	"zmin=0",			";Zmin for slab-anim (Choose this value acording snap)",
	"zmax=1",			";Zmax for slab-anim (Choose this value acording snap)",
//
#endif
	"usingcolumns=1:2",			";Columns to plot or to average", ":uc",
	"usingrows=1:10",			";Rows to plot or to average (interval of lines)", ":ur",
    "xlabel=",					";X-axis label", ":xl",
    "ylabel=",					";Y-axis label", ":yl",
	"plotlabel=",				";Plot label", ":pl",
    "labelfontsize=1",			";Label font scale size", ":lfs",
	"plotjoined=true",			";Plot with points joined with a line", ":pj",
	"withdots=false",			";Plot with dots (set plotjoined to false)", ":wd",
	"withsymbols=true",			";Plot with symbols", ":ws",
	"symbolsize=1",				";Symbol size", ":ss",
//
    "labelfontweight=1",	";Label font weight",
	"nlsize=1",				";Numeric label font size",
	"linewidth=1",			";Line width",
	"axeswidth=1",			";Axes width",
	"symbolweight=1",		";Symbol weight",
	"symbolcolor=1",		";Symbol color",
//
	"aspectratio=",				";Page aspect ratio [1.5] (def: same as output device)", ":a",
	"dev=xwin",					";Output device name (xwin, aqt, ...)",
	"geometry=640x480",			";Window size, in pixels", ":geo",
	"ori=0",					";Plot orientation (0,1,2,3=landscape,portrait,seascape,upside-down)",
	"background=",				";Background color (000000=black, FFFFFF=white)", ":bg",
	"ncol0=",					";Number of colors to allocate in cmap 0 (upper bound)",
	"ncol1=",					";Number of colors to allocate in cmap 1 (upper bound)",
    "Version=0.2",              ";Mar 2005-2009",
    NULL,
};


#endif /* ! _cmdline_defs_h */
