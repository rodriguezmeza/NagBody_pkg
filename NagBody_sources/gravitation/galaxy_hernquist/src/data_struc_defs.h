/*==============================================================================
	HEADER: data_struc_defs.h			[galaxy_hernquist]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: main, direct_gravcalc,
					nbody_n2_io, startrun, timestep
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:  November 2008;
	Copyright: (c) 2005-2014 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#if !defined(global)                    // global def question must be here
#  define global extern
#endif


//#include "globaldefs.h"

#ifndef _data_struc_defs_h
#define _data_struc_defs_h


//typedef struct
global struct io_header_1
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
} header1;
//io_header_1; //

//global io_header_1 header1;

//typedef struct
global struct particle_data
{
    float  Pos[3];
    float  Vel[3];
    float  Mass;
    int    Type;
    
    float  Rho, U, Temp, Ne;
} *P;
//particle_data; //
//particle_data *P;




#define IPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&(param);												\
  id[nt++]=INT;}

#define RPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=DOUBLE;}

#define BPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)									\
  {strcpy(tag[nt],paramtext);										\
  param=(string) malloc(n);											\
  addr[nt]=param;													\
  id[nt++]=STRING;}

#endif // ! _data_struc_defs_h

