/*==============================================================================
	HEADER: cmdline_defs.h			[models]
	Written by: Mario A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: model (main)
	External headers: None
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

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"Code for generating initial conditions of an N-General-Body system"
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",			";Parameter input file. Overwrite what follows",
    "in=",                  ";Input file or files with N-Body data",
    "infmt=",               ";N-Body data file format",
    "out=snapout.dat",		";Output file of N-body frames", ":o",
    "outfmt=",				";Snap output file format of N-body system", ":ofmt",
    "model-type=plummer-finite",	";Model to generate", ":mt",
    "options=",             ";Various control options (not in use by now)", ":opt",
    "nbody=4096",			";Number of bodies of the model",
    "Mtotal=1.0",			";Total mass of the system",
    "Rmax=1.0",				";Maximum radius of the system",
    "vcmx=0.",				";x-component of the velocity of the center of mass of the system",
    "vcmy=0.",				";y-component of the velocity of the center of mass of the system",
    "vcmz=0.",				";z-component of the velocity of the center of mass of the system",
    "cmx=0.",				";x-component of the position of the center of mass of the system",
    "cmy=0.",				";y-component of the position of the center of mass of the system",
    "cmz=0.",				";z-component of the position of the center of mass of the system",
    "absvel=0.",			";Maximum absolute velocity for a random motion of the bodies",
    "omega0=0.",			";Angular velocity of the system around z-axis (omega0=0.698783 for isothermal collapse of a sphere)",
    "a_p=0.1",				";Amplitud of the azimuthal perturbation (Useful isothermal collapse of a sphere)",
    "m_p=2.0",				";Mode of the azimuthal perturbation (Useful isothermal collapse of a sphere)",
    "SoundSpeed=0.394338",	";Sound speed (Useful for isothermal collapse of a sphere)",
    "a_spheroid=2.0",		";Spheroid x-axis (Bonnort-Ebert case default, also give Mtotal=40)",
    "b_spheroid=2.0",		";Spheroid x-axis (Bonnort-Ebert case default)",
    "c_spheroid=1.0",		";Spheroid x-axis (Bonnort-Ebert case default)",
    "factor=1.0",			";Numerical factor to convert one type of snapdata",
    "seed=-1",              ";Random number seed for the model",
    "Version=0.2",          ";Mario A. Rodriguez-Meza 1999-2018",
    NULL,
};

#endif // !_cmdline_defs_h
