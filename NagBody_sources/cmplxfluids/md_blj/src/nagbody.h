/* ==============================================================================
	HEADER: nagbody.h				[NagBody]
	Written by: Mario A. Rodriguez-Meza
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
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _nagbody_h
#define _nagbody_h

#include "nagbody_struct.h"
#include "nagbody_proto.h"

// MACROS - Definitions to add a parameter in the scheme of parameterfile ------

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

// -----------------------------------------------------------------------------

// MACROS - Definitions to handle periodic boundaries and cells method...-------

#define DO_CELL(j, m)  for (j = gdforce.cellList[m]; j >= 0;		\
								j = gdforce.cellList[j])

#define DO_CELL_ptr(j, m)  for (j = gdforce->cellList[m];			\
							j >= 0; j = gdforce->cellList[j])

/*
#define VWrap(v, t)													\
	v[t] = v[t] - ( (real) (nint(v[t]/Box[t])) )*Box[t]
*/
// Alternative definition for VWrap

#define VWrap(v, t)													\
   if (v[t] >= 0.5 * gdforce.Box[t])  v[t] -= gdforce.Box[t];		\
   else if (v[t] < -0.5 * gdforce.Box[t]) v[t] += gdforce.Box[t]

#define VWrap_ptr(v, t)												\
   if (v[t] >= 0.5 * gdforce->Box[t])  v[t] -= gdforce->Box[t];		\
   else if (v[t] < -0.5 * gdforce->Box[t]) v[t] += gdforce->Box[t]

// Esta definicion esta casada con la definicion de m2v y shift 
// dentro de la rutina cellsmethod ... en el codigo.
// Cambios en esas definiciones alteraran el funcionamiento de
// VCellWrap ...
#define VCellWrap(t)												\
   if (m2v[t] >= gdforce.cells[t]) {								\
     m2v[t] = 0;													\
     shift[t] = gdforce.Box[t];										\
   } else if (m2v[t] < 0) {											\
     m2v[t] = gdforce.cells[t] - 1;									\
     shift[t] = - gdforce.Box[t];									\
   }

#define VCellWrap_ptr(t)											\
   if (m2v[t] >= gdforce->cells[t]) {								\
     m2v[t] = 0;													\
     shift[t] = gdforce->Box[t];									\
   } else if (m2v[t] < 0) {											\
     m2v[t] = gdforce->cells[t] - 1;								\
     shift[t] = - gdforce->Box[t];									\
   }

#if NDIM == 2
#define VWrapAll(v)													\
   {VWrap (v, 0);													\
   VWrap (v, 1);}
#define VWrapAll_ptr(v)												\
   {VWrap_ptr (v, 0);												\
   VWrap_ptr (v, 1);}
#define VCellWrapAll()												\
   {VCellWrap (0);													\
   VCellWrap (1);}
#define VCellWrapAll_ptr()											\
   {VCellWrap_ptr (0);												\
   VCellWrap_ptr (1);}
#define OFFSET_VALS													\
   {{0,0}, {1,0}, {1,1}, {0,1}, {-1,1}}
#define N_OFFSET  5
#endif

#if NDIM == 3

#define VWrapAll1D(v)												\
   {VWrap (v, 0);}

#define VWrapAll2D(v)												\
   {VWrap (v, 0);													\
   VWrap (v, 1);}

#define VWrapAll(v)													\
   {VWrap (v, 0);													\
   VWrap (v, 1);													\
   VWrap (v, 2);}
#define VWrapAll_ptr(v)												\
   {VWrap_ptr (v, 0);												\
   VWrap_ptr (v, 1);												\
   VWrap_ptr (v, 2);}
#define VCellWrapAll()												\
   {VCellWrap (0);													\
   VCellWrap (1);													\
   VCellWrap (2);}
#define VCellWrapAll_ptr()											\
   {VCellWrap_ptr (0);												\
   VCellWrap_ptr (1);												\
   VCellWrap_ptr (2);}
#define OFFSET_VALS													\
   {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0},					\
	{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1},			\
	{-1,-1,1}, {0,-1,1}, {1,-1,1}									\
   }

#define N_OFFSET  14
#endif

// -----------------------------------------------------------------------------

// TAGS FOR I/O CONTROL --------------------------------------------------------

#define IO_NULL_FMT							1

#define IO_SNAP_BLJ_FMT						11
#define IO_SNAP_BLJ_PV_FMT					12
#define IO_SNAP_BLJ_BIN_FMT					13
#define IO_SNAP_TLJ_FMT						14
// -----------------------------------------------------------------------------

#endif	// ! _nagbody_h
