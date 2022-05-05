/* ==============================================================================
	HEADER: cmdline_defs.h			[analysis_grav]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main)
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
===============================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code for data analysis of a autogravitating system simulation"
#define HEAD3	"-----"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			";Parameter input file. Overwrite what follows",
//
    "in=",						";Input N-Body snap base name (snaps/snap%04d)",
    "infmt=",					";N-Body snap file format (not active by now)",
    "out=",						";Output base name file of N-body frames", ":o",
    "outfmt=",					";Output file format of N-body frames (is the extension)", ":ofmt",
	"basedir=.",				";Base directory for output files", ":bdir",
	"isnap=0",					";Initial snap number",
	"fsnap=100",				";Final snap number",
    "options=",                 ";Various control options (save, snap_cgm, rdf_cgm, vel_cgm)", ":opt",
    "analysis_type=",			";Data analysis type to be done (snap-anim, snap-animtrajectory... see the manual)", ":at",
//
// Seran activados cuando se necesite calcular la fuerza entre particulas.
//	"usequad=false",		";use quadrupoles",
//	"theta=1",		";force resolution parameter",
//	"eps=0.025",		";smoothing length",
//
	"stepAvgRhoTheta=1",			";number of rho_theta histograms to average",
	"sizeHistRhoTheta=50",			";array size for rho_theta histogram",
	"RhoDeltaZ=0.2",			";thickness of the slab",
	"RhoR=0.2",			";Radius of the slab",
	"RhoDeltaR=0.2",			";Radial size of the slab",
	"RMax=0.5",			";Maximum radius of the slab",
//
	"xmin=0",			";Xmin for slab-anim (Choose this value acording snap)",
	"xmax=1",			";Xmax for slab-anim (Choose this value acording snap)",
	"ymin=0",			";Ymin for slab-anim (Choose this value acording snap)",
	"ymax=1",			";Ymax for slab-anim (Choose this value acording snap)",
	"zmin=0",			";Zmin for slab-anim (Choose this value acording snap)",
	"zmax=1",			";Zmax for slab-anim (Choose this value acording snap)",
//
	"stepAvgRho=1",			";number of rho histograms to average",
	"sizeHistRho=50",			";array size for rho histogram",
//
	"stepAvgVcR=1",				";number of block of steps to average rot velocity",
	"sizeHistVcR=50",			";array size for rot velocity histogram",
	"rangeR=0.4",				";range of rot velocities for histogram",
	"ThetaMin=0.",				";Azimuthal min angle",
	"ThetaMax=3.1415",			";Azimuthal max angle",
//
//	"computeTransport=false",		";Compute transport properties",
//	"stepAcf=3",			";number of step jumps to save an Acf measurement",
	"limitAcfAv=20",			";number of Acf measurements to average",
	"nBuffAcf=10",			";size of buffer to save an Acf measurements",
	"nValAcf=20",			";number of values to save of Acf measurement",
//
//    "nbody=512",				";Number of bodies for test run (must be the same as in snaps)",
//    "dtime=1/256",				";Integration time step",
//    "dtout=5/256",              ";Data output time step",
//
	"bodiesID=",				";Body IDs to use with bodyies-anim (b1,b2,...)", ":bid",
	"bodiesSets=",			";Set of bodies Ids to plot",
	"reductionFac=1.0",			";Reduction Factor for snap_lessnbody option",
	"foffile=",					";File with friend-of-friends info",
	"gravityConstant=1.0",		";Gravitational constant value", ":G",
//
	"xrange=autoscale",			";Range in x-axis (xmin:xmax)", ":xr",
	"yrange=autoscale",			";Range in y-axis (ymin:ymax)", ":yr",
	"zrange=autoscale",			";Range in z-axis (zmin:zmax)", ":zr",
	"usingcolumns=1:2",			";Columns to plot or to average", ":uc",
	"usingrows=1:10",			";Rows to plot or to average (interval of lines)", ":ur",
    "xlabel=",					";X-axis label", ":xl",
    "ylabel=",					";Y-axis label", ":yl",
	"plotlabel=",				";Plot label", ":pl",
    "labelfontsize=1",			";Label font scale size", ":lfs",
	"withdots=false",			";Plot with dots (set plotjoined to true)", ":wd",
	"withsymbols=true",			";Plot with symbols", ":ws",
	"symboltype=1",			";Symbol type", ":st",
	"symbolsize=1",				";Symbol size", ":ss",
//
    "labelfontweight=1",	";Label font weight",
	"nlsize=1",				";Numeric label font size",
	"linewidth=1",			";Line width",
	"axeswidth=1",			";Axes width",
	"symbolweight=1",		";Symbol weight",
	"symbolcolor=1",		";Symbol color",
//
	"a=",						";Page aspect ratio [1.5] (def: same as output device)",
	"dev=xwin",					";Output device name (xwin, aqt, ...)",
	"geo=640x480",				";Window size, in pixels (e.g. -geo 400x300)",
	"ori=0",					";Plot orientation (0,1,2,3=landscape,portrait,seascape,upside-down)",
	"bg=",						";Background color (000000=black, FFFFFF=white)",
	"ncol0=",					";Number of colors to allocate in cmap 0 (upper bound)",
	"ncol1=",					";Number of colors to allocate in cmap 1 (upper bound)",
    "Version=0.2",              ";Mar 2005-2011",
    NULL,
};


#endif /* ! _cmdline_defs_h */
