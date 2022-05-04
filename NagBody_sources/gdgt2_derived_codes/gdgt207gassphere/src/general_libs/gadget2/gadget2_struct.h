/* ==============================================================================
	HEADER: gadget2_struct.h			[NagBody]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main), start_run, time_step,
					forcecalc, md_lj_tree_io
	External headers:
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

#ifndef _gadget2_struct_h
#define _gadget2_struct_h


// ------------------------START STRUCTURE DEFINITIONS--------------------------

typedef struct {
	FLOAT s[3];                     
	FLOAT vs[3];                    
	FLOAT mass;                     
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	int   bitflags;                 
#else
	FLOAT maxsoft;                  
	
#endif
#endif
} gadget2_DomainNODE, *gadget2_DomainNODE_ptr;

global gadget2_DomainNODE_ptr DomainMoment;


typedef struct {
	int Daughter;                   
	int Pstart;                     
	int Blocks;                     
	int Leaf;                       
	peanokey Size;                  
	peanokey StartKey;              
	long long Count;                
} gadget2_topnode_data

global gadget2_topnode_data *TopNodes;                       


//------------------------------------------------------------------------------

#endif	// ! _gadget2_struct_h
