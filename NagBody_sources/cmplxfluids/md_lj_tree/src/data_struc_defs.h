/* ==============================================================================
!	HEADER: data_struc_defs.h													!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date: February 2005												!
!	Purpose: Definition of N-Body data structure								!
!	Language: C																	!
!	Use: '#include "...."														!
!	Use in routines and functions: md_lj_tree (main), tree_ljforcecalc,			!
!					md_lj_tree_io, start_run, time_step							!
!	External headers: None														!
!	Comments and notes:															!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.							!
!		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar											!
!																				!
!	Major revisions:															!
!	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.							!
! =============================================================================*/
 
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
#define CELL 02                 

typedef struct {
    node bodynode;              
    vector vel;                 
    vector acc;                 
    real phi;
	int Id;
} body, *bodyptr;

#define Vel(x)    (((bodyptr) (x))->vel)
#define Acc(x)    (((bodyptr) (x))->acc)
#define Phi(x)    (((bodyptr) (x))->phi)
#define Id(x)    (((bodyptr) (x))->Id)

#define NSUB (1 << NDIM)        
 
typedef struct {
    node cellnode;              
    real rcrit2;                
    nodeptr more;                  
    union {
        nodeptr subp[NSUB];     
        matrix quad;            
    } sorq;
} cell, *cellptr;
 
#define Rcrit2(x) (((cellptr) (x))->rcrit2)

#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->sorq.subp)
#define Quad(x)   (((cellptr) (x))->sorq.quad)

#if !defined(global)			// global def question must be here
#  define global extern
#endif


#define DO_BODY(p,start,finish)  for (p = start; p < finish; p++)
#define DO_DESCENDENTS(q,p)		for (q = More(p); q != Next(p); q = Next(q))
#define DO_COORD(k)				for (k=0; k<NDIM; k++)

#endif /* ! _data_struc_defs_h */

