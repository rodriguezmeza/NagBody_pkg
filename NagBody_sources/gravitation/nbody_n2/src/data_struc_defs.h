/*==============================================================================
	HEADER: data_struc_defs.h			[nbody_n2]
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
        http://www.inin.gob.mx/

	Major revisions:  November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.																		!
==============================================================================*/
 
#ifndef _data_struc_defs_h
#define _data_struc_defs_h

typedef struct _node {
    short type;                 
    bool update;                
    real mass;                  
    vector pos;                 
    struct _node *next;            
} node, *nodeptr;
 
#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Mass(x)   (((nodeptr) (x))->mass)
#define Pos(x)    (((nodeptr) (x))->pos)
#define Next(x)   (((nodeptr) (x))->next)

#define BODY 01

typedef struct {
    node bodynode;              
    vector vel;                 
    vector acc;                 
    real phi;                   
} body, *bodyptr;

#define Vel(x)    (((bodyptr) (x))->vel)
#define Acc(x)    (((bodyptr) (x))->acc)
#define Phi(x)    (((bodyptr) (x))->phi)

#if !defined(global)					// global def question must be here
#  define global extern
#endif

#define DO_BODY(p,start,finish)  for (p = start; p < finish; p++)
#define DO_COORD(k)				for (k=0; k<NDIM; k++)

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

