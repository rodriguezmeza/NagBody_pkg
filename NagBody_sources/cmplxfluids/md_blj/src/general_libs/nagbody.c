/*==============================================================================
	MODULE: nagbody.c		[NagBody]
	Written by: M.A. Rodriguez-Meza
	Starting date:	January, 2005
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:	stdinc.h, mathfns.h, vectmath
						vectmath.h, getparam.h
						types.h, stat.h, inout.h
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "stdinc.h"
#include "mathfns.h"
#include "inout.h"
#include "vectdefs.h"
#include "vectmath.h"
#include "nagbody.h"
#include "physconstants.h"

#include <string.h>
// #include <strings.h>							// For unix
//#include "../../../../General_libs/strings.h"	// For Visual C

#include <sys/stat.h>


//------------- COMIENZA BLOQUE DE RUTINAS DE RESERVA DE MEMORIA ---------------
void AllocateMemory(short allocate_mode)					// CHECK 2D --- OK!!
{
    static bool firstcall = TRUE;

	switch(allocate_mode) {
		case 1:
			allocate_memory(); break;
		case 2:
			allocate_memory_double(); break;
		case 3:											// Output single
			if (!firstcall) return;
			firstcall = FALSE;                     
			allocate_memory(); break;
		case 4:											// Output double
			if (!firstcall) return;
			firstcall = FALSE;                     
			allocate_memory_double(); break;
		case 5:											// io special particle structure
			allocate_memory_long(); break;				// N > 10^6 ...
	}
}

void allocate_memory(void)								// CHECK 2D --- OK!!
{
//    static bool firstcall = TRUE;
	int bytes,bytes_tot=0;
//	long int bytes,bytes_tot=0;		// Para poder trabajar con mas de 16x10^6 particulas


// SE COMENTAN PARA QUE FUNCIONE analysis_grav y analysis_galaxy en cuanto
// a intercambio de formatos o lectura de formatos tipo gadget11-ascii por ejemplo...
// Pero este par de lineas fue puesto con un proposito que mas tarde vere de que
// se trato...
//
//	if (!firstcall) return;
//	firstcall = FALSE;                     

	if(All.MaxPart>0) {
		if(!(P_data=malloc(bytes=All.MaxPart*sizeof(particle_data)))) {
			printf("failed to allocate memory for `P_data' (%ld bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		if(!(PTimeTree=malloc(bytes=All.MaxPart*sizeof(timetree_data)))) {
			printf("failed to allocate memory for `PTimeTree' (%d bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		P= P_data-1;
		PTimeTree--;

		printf("\nAllocated %g MByte for particle storage.\n\n",
			bytes_tot/(1024.0*1024.0));
    }
// IBERO COMIENZO //////////////////////////////////////////////////////////
/*
// SE PERDERA MEMORIA ... COMENTAR CUANDO NO SE OCUPE ...    
	if(All.MaxPart>0) {
		if(!(P_IBERO_data=malloc(bytes+=All.MaxPart*sizeof(particle_data)))) {
			printf("failed to allocate memory for `P_data' (%ld bytes).\n",
                   bytes);
			endrun(1);
		}
		bytes_tot+=bytes;
        
//		if(!(PTimeTree=malloc(bytes=All.MaxPart*sizeof(timetree_data)))) {
//			printf("failed to allocate memory for `PTimeTree' (%d bytes).\n",
//                   bytes);
//			endrun(1);
//		}
//		bytes_tot+=bytes;
        
		P_IBERO= P_IBERO_data-1;
//		PTimeTree--;
        
		printf("\nAllocated %g MByte for particle storage.\n\n",
               bytes_tot/(1024.0*1024.0));
    }
*/
// IBERO FIN ///////////////////////////////////////////////////////////

	if(All.MaxPartSph>0) {
		bytes_tot=0;

		if(!(SphP_data=malloc(bytes=All.MaxPartSph*sizeof(sph_particle_data)))){
			printf("failed to allocate memory for `SphP_data' (%d bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		SphP= SphP_data-1; 

		printf("Allocated %g MByte for storage of SPH data.\n\n",
			bytes_tot/(1024.0*1024.0));
	}
}

// Particle data structure to manipulate I/O
// N > 10^6 purpose...
void allocate_memory_long(void)
{
	int bytes,bytes_tot=0;
	if(All.MaxPart>0) {
		if(!(P_data_long=malloc(bytes=All.MaxPart*sizeof(particle_data_long)))) {
			printf("failed to allocate memory for `P_data_long' (%ld bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		if(!(PTimeTree=malloc(bytes=All.MaxPart*sizeof(timetree_data)))) {
			printf("failed to allocate memory for `PTimeTree' (%d bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		P_long= P_data_long-1;
		PTimeTree--;

		printf("\nAllocated %g MByte for particle storage.\n\n",
			bytes_tot/(1024.0*1024.0));
    }

	if(All.MaxPartSph>0) {
		bytes_tot=0;

		if(!(SphP_data=malloc(bytes=All.MaxPartSph*sizeof(sph_particle_data)))){
			printf("failed to allocate memory for `SphP_data' (%d bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		SphP= SphP_data-1; 

		printf("Allocated %g MByte for storage of SPH data.\n\n",
			bytes_tot/(1024.0*1024.0));
	}
}
//


void allocate_memory_double(void)						// CHECK 2D --- OK!!
{
	int bytes,bytes_tot=0;
//	long int bytes,bytes_tot=0;		// Para poder trabajar con mas de 16x10^6 particulas

	if(All.MaxPart>0) {
		if(!(P_data_double=malloc(bytes=All.MaxPart
										*sizeof(particle_data_double)))) {
			printf("failed to allocate memory for `P_data' (%d bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		if(!(PTimeTree=malloc(bytes=All.MaxPart*sizeof(timetree_data)))) {
			printf("failed to allocate memory for `PTimeTree' (%d bytes).\n",
				bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		P_double= P_data_double-1;
		PTimeTree--;

		printf("\nAllocated %g MByte for particle storage.\n\n",
			bytes_tot/(1024.0*1024.0));
	}

	if(All.MaxPartSph>0) {
		bytes_tot=0;

		if(!(SphP_data_double=malloc(bytes=All.MaxPartSph
										*sizeof(sph_particle_data_double)))) {
			printf("failed to allocate memory for `SphP_data' (%d bytes).\n",
					bytes);
			endrun(1);
		}
		bytes_tot+=bytes;

		SphP_double= SphP_data_double-1; 

		printf("Allocated %g MByte for storage of SPH data.\n\n",
			bytes_tot/(1024.0*1024.0));
	}
}

void free_memory(void)									// CHECK 2D --- OK!!
{
	if(All.MaxPart>0) {
		printf("\nFreeing memory\n\n");
		free(P_data);
//    free(P);
//	free(PTimeTree);				// ADDED ... CHECK!!!!
	}

	if(All.MaxPartSph>0)
		free(SphP_data);
}

//-------------- TERMINA BLOQUE DE RUTINAS DE RESERVA DE MEMORIA ---------------


//------------------- COMIENZA BLOQUE DE RUTINAS GENERALES ---------------------

void code_endrun(stream outlog, real cpuinit)			// CHECK 2D --- OK!!
{
	fprintf(stdout,"\n\nTotal running cpu time: %gm\n\n",cputime()-cpuinit);
	fclose(outlog);
}

//-------------------- TERMINA BLOQUE DE RUTINAS GENERALES ---------------------

