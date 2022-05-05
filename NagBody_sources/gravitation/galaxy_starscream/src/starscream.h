//
// ESTE ES "data_struc_defs.h"
//
/*-----------------------------------------------------------------------------
/
/ Filename: starscream.h
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: This is the header file for Starscream. Starscream creates
/              galaxies.
/
/	       Starscream uses the GNU Scientific Library (GSL). You can
/	       download the GSL source code from:
/
/		http://www.gnu.org/software/gsl
/
/	       or replace it with another math library.
/
/ Copyright Information:
/
/ Copyright (c) 2008, 2009       Jay Jay Billings
/
/ This program is free software; you can redistribute it and/or modify
/ it under the terms of the GNU General Public License as published by
/ the Free Software Foundation; either version 2 of the License, or
/ (at your option) any later version.
/
/ This program is distributed in the hope that it will be useful,
/ but WITHOUT ANY WARRANTY; without even the implied warranty of
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/ GNU General Public License for more details.
/
/ You should have received a copy of the GNU General Public License
/ along with this program; if not, write to the Free Software
/ Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
/
/ The license is also available at:
/
/		http://www.gnu.org/copyleft/gpl.html .
/
/ Date: 2009/06/07
/
*///---------------------------------------------------------------------------

//
// MAR :: NagBody
#ifndef _starscream_defs_h
#define _starscream_defs_h
// END :: NagBody

/* Header files to include							      */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
//#include <fftw3.h>

/* This is a type definition of a galaxy. The thought here is to create galaxies
as objects and, hopefully, make the code extremely clean. It is essentially a
a collection of arrays.
        lambda: The spin paramater.
        m_d: The mass fraction of the disk relative to the total mass.
        j_d: The angular momentum fraction of the disk relative to the total.
        total_mass: The total mass of the galaxy.
        disk_mass: The mass of the disk.
        disk_scale_length: The scale length of the disk.
        halo_mass: The mass of the halo.
        halo_scale_length: The scale length of the halo.
        halo_concentration: The concentration of the halo. 
        halo_a_value: The value that relates a Hernquist profile to a NFW profile
                      through the halo scale length and halo concentration.
        v200: The halo's virial velocity.
        r200: The halo's virial radius.
	x: The x positions of the particle.
        y: The y positions of the particle.
        z: The z positions of the particle.
        vel_x, vel_y, vel_z: The velocity components of the particles.
        mass: The masses of the particle.					
        potential: The potential of the disk evaluated at the grid points. The
		   size of this array should be gx*gy*gz.	      
        storage1: A special variable for storing important info. Could be anything!   
	id: The particle id number. 
	num_part: An integer array containing the number of particles in the disk, [0],
		  and in the halo, [1], and the total number of particles (optional). */
typedef struct {
    double lambda;
    double m_d;
    double j_d;
    double total_mass;
    double disk_mass;
    double disk_scale_length;
    double halo_mass;
    double bulge_mass;
    double halo_scale_length;
    double halo_concentration;
    double halo_a_value;
    double v200;
    double r200;
    double *x;
    double *y;
    double *z;
    double *vel_x;
    double *vel_y;
    double *vel_z;
    double *mass;
    double ***potential;
    double space[3];
    double *storage1;
    int *id;
    int num_part[3];
    int potential_defined;
} galaxy;

//Gadget2-style header for Gadget2 snapshots.
//
// MAR :: NagBody
//struct io_header_1 {
typedef struct io_header_1 {
   int npart[6];
   double mass[6];
   double time;
   double redshift;
   int flag_sfr;
   int flag_feedback;
   int npartTotal[6];
   int flag_cooling;
   int num_files;
   double BoxSize;
   double Omega0;
   double OmegaLambda;
   double HubbleParam;
   char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
//} header1;
} io_header_1;

//Gadget2-style particle data structure.
//struct particle_data {
typedef struct {
   float Pos[3];
   float Vel[3];
   float Mass;
   int Type;

   float Rho, U, Temp, Ne;
//} *P;
} particle_data, *particle_data_ptr;

#if !defined(global)					// global def question must be here
#  define global extern
#endif

global particle_data_ptr P;
global io_header_1 header1;
// END :: NagBody


//
// MAR :: NagBody :: "global" added
//
//Particle IDs.
global int *Id;

//Time and Redshift variables. Both will probably be 0.
global double Time, Redshift;

//Good numbers to have on file: Total number of particles and total number of gas
//particles.
global int NumPart, Ngas;

// Global variables for the random number environment.
global const gsl_rng_type *T;
global gsl_rng *r;
global int random_number_set;


// A variable to check if the galaxy is being copied.
global int a_copy;

// A list of universal constants. They are set in create_galaxy(). 
global double pi, G;

// The particle-mesh grid size and the Green's function and potential 
// storage buffers. Global to keep from hitting the stack limit for large grid 
// size.
global int Ng;

//
// END :: NagBody


/*----- These are the function prototypes for Starscream. 		-----*/

// Initialization and destruction functions
int create_galaxy(galaxy *, int [], double, double, double, double, 
                 double, int, double, int);
void destroy_galaxy(galaxy *, int);
void destroy_galaxy_system(galaxy *, int);
int create_galaxy_system(galaxy *, galaxy *, galaxy *); 

// Input, output, and manipulation functions
int write_gadget_ics(galaxy *, char *);
void copy_galaxy(galaxy *, galaxy *, int);
void set_orbit_parabolic(galaxy *, galaxy *, double, double);
int load_snapshot(char *, int);
int allocate_memory();
int reordering();
int unload_snapshot();


//
//
#define IPName(param,paramtext)									\
{strcpy(tag[nt],paramtext);										\
    addr[nt]=&(param);											\
    id[nt++]=INT;}

#define RPName(param,paramtext)									\
{strcpy(tag[nt],paramtext);										\
    addr[nt]=&param;											\
    id[nt++]=DOUBLE;}

#define BPName(param,paramtext)									\
{strcpy(tag[nt],paramtext);										\
    addr[nt]=&param;											\
    id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)								\
{strcpy(tag[nt],paramtext);										\
    param=(string) malloc(n);									\
    addr[nt]=param;                                             \
    id[nt++]=STRING;}

//
// MAR :: NagBody
#endif	/* ! _starscream_defs_h	*/
//