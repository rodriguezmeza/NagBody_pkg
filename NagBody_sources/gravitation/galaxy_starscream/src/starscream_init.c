/*-----------------------------------------------------------------------------
/
/ Filename: starscream_init.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: These are the initialization routines for Starscream. 
/              They include routines to set up the memory structure and 
/              allocate particle positions and velocities. Starscream 
/              creates galaxies. 
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

// Get the configuration header
//
// MAR :: NagBody
//#include "../config.h"
// END :: NagBody

// Include the starscream header
#include "starscream.h"

// Create a galaxy
int create_galaxy(galaxy *gal, int parts[3], double m_d, double j_d, double lambda, 
                 double c, double v200, int Ngrid, double space, int info) {

    int i,j,num_part;
    long seed;

//    These values are global, but assigned values here.
    pi = 3.14159274101257;     // pi = 4.0*atan(1.0)
    //G = 6.67428E-8;            // G = 6.67428E-8 cm^3 g^-1 s^-2
    G = 1;
    Ng = 2*Ngrid;

    // Create the random number generator environment.
    if (random_number_set != 1) {
       seed = time(NULL);
       gsl_rng_env_setup();
       T = gsl_rng_default;
       r = gsl_rng_alloc(T);
       gsl_rng_set(r,seed);
       random_number_set = 1;
    }

    // Set the parameters.
    gal->num_part[0] = parts[0];
    gal->num_part[1] = parts[1];
    gal->num_part[2] = parts[2];
    gal->m_d = m_d;
    gal->j_d = j_d;
    gal->lambda = lambda;
    gal->halo_concentration = c;
    gal->v200 = v200;
    gal->space[0] = space;
    
    gal->num_part[3] = gal->num_part[0] + gal->num_part[1]+gal->num_part[2];
    num_part = gal->num_part[0]+gal->num_part[1]+gal->num_part[2];


    // Allocate particle id numbers array.
    if (!(gal->id=calloc(num_part,sizeof(int)))) {
       fprintf(stderr,"Unable to allocate particle ID numbers.\n");
       return -1;
    }
    // Allocate x coordinates for all the particles.
    if (!(gal->x=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle x coordinates.\n");
       return -1;
    }
    // Allocate y coordinates for all the particles.
    if (!(gal->y=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle y coordinates.\n");
       return -1;
    }
    // Allocate z coordinates for all the particles.
    if (!(gal->z=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle z coordinates.\n");
       return -1;
    }
    // Allocate x velocities for all the particles.
    if (!(gal->vel_x=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle x coordinates.\n");
       return -1;
    }
    // Allocate y velocities for all the particles.
    if (!(gal->vel_y=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle y coordinates.\n");
       return -1;
    }
    // Allocate z velocities for all the particles.
    if (!(gal->vel_z=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle z coordinates.\n");
       return -1;
    }
    // Allocate masses for all the particles.
    if (!(gal->mass=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle masses.\n");
       return -1;
    }
    // Turn on and allocate the potential grid, start with x-axis
    gal->potential_defined = 1;
    if (!(gal->potential=calloc(Ng,sizeof(double *)))) {
       fprintf(stderr,"Unable to create potential x axis.\n");
       return -1;
    }
    for (i = 0; i < Ng; ++i) {
        // y-axis
        if (!(gal->potential[i] = calloc(Ng,sizeof(double *)))) {
           fprintf(stderr,"Unable to create potential y axis.\n");
           return -1;
        }
        // z-axis
        for (j = 0; j < Ng; ++j) {
            if (!(gal->potential[i][j] = calloc(Ng,sizeof(double)))) {
               fprintf(stderr,"Unable to create potential z axis.\n");
               return -1;
            }
        }
    }

    return 0;
}


// Destroy a galaxy
void destroy_galaxy(galaxy *gal, int info) {

     int i,j;

     if (info != 0) {
        fprintf(stderr,"Destroying galaxy...");
     }

     // Deallocate the potential grid to be really nice to the memory.
     for (i = 0; i < Ng; ++i) {
         for (j = 0; j < Ng; ++j) {
             free(gal->potential[i][j]);
         }
         free(gal->potential[i]);
     }
     free(gal->potential);

     // Deallocate all the small parts of the galaxy to be nice to the memory.
     free(gal->id); free(gal->x); free(gal->y); free(gal->z); free(gal->mass);
     free(gal->vel_x); free(gal->vel_y); free(gal->vel_z);

     if (info != 0) {
        fprintf(stderr," All memory deallocated.\n");
     }
     
     return;
}

// Destroy a system of galaxies.
void destroy_galaxy_system(galaxy *gal, int info) {

     if (info != 0) {
        fprintf(stderr,"Destroying galaxy system...");
     }

     free(gal->id); free(gal->x); free(gal->y); free(gal->z); free(gal->mass);
     free(gal->vel_x); free(gal->vel_y); free(gal->vel_z);

     if (info != 0) {
        fprintf(stderr," All memory deallocated.\n");
     }
     
     return;
}

// In order to perform a galaxy collision, it is necessary to combine 
// two galaxies into one set of orbiting galaxies. The orbits should 
// be defined by a user-defined orbit or the set_orbit_parabolic()
// function and the following function should be called to combine 
// the galaxies into one object.
// 
// It is necessary to allocate the different parts of the memory here
// because this galactic system does not have quantities like halo
// or disk scale lengths. These are quantities in the parent galaxy,
// not the galaxy-orbit-galaxy combination.
int create_galaxy_system(galaxy *gal_1, galaxy *gal_2, galaxy *gal_3) {

     int i, a, num_part;

     // Create the galaxy system first!
     gal_3->num_part[0] = gal_1->num_part[0] + gal_2->num_part[0];
     gal_3->num_part[1] = gal_1->num_part[1] + gal_2->num_part[1];
     gal_3->num_part[2] = gal_1->num_part[2] + gal_2->num_part[2];
     num_part = gal_3->num_part[0]+gal_3->num_part[1]+gal_3->num_part[2];   
     // Allocate particle id numbers array.
     if (!(gal_3->id=calloc(num_part,sizeof(int)))) {
        printf("Unable to allocate particle ID numbers.\n");
        return -1;
     }
     // Allocate x coordinates for all the particles.
     if (!(gal_3->x=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle x coordinates.\n");
        return -1;
     }
     // Allocate y coordinates for all the particles.
     if (!(gal_3->y=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle y coordinates.\n");
        return -1;
     }
     // Allocate z coordinates for all the particles.
     if (!(gal_3->z=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle z coordinates.\n");
        return -1;
     }
     // Allocate x velocities for all the particles.
     if (!(gal_3->vel_x=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle x coordinates.\n");
        return -1;
     }
     // Allocate y velocities for all the particles.
     if (!(gal_3->vel_y=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle y coordinates.\n");
        return -1;
     }
     // Allocate z velocities for all the particles.
     if (!(gal_3->vel_z=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle z coordinates.\n");
        return -1;
     }
     // Allocate masses for all the particles.
     if (!(gal_3->mass=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle masses.\n");
        return -1;
     }
     // Turn off the galaxy potential.
     gal_3->potential_defined = 1;

     // Copy all of the galaxy information.
     // Disk 1
     a = 0;
     for (i = 0; i < gal_1->num_part[0]; ++i) {
         gal_3->mass[a] = gal_1->mass[i];
         gal_3->id[a] = i;
         gal_3->x[a] = gal_1->x[i];
         gal_3->y[a] = gal_1->y[i];
         gal_3->z[a] = gal_1->z[i];
         gal_3->vel_x[a] = gal_1->vel_x[i];
         gal_3->vel_y[a] = gal_1->vel_y[i];
         gal_3->vel_z[a] = gal_1->vel_z[i];
         ++a;
     }
     // Disk 2
     for (i = 0; i < gal_2->num_part[0]; ++i) {
         gal_3->mass[a] = gal_2->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_2->x[i];
         gal_3->y[a] = gal_2->y[i];
         gal_3->z[a] = gal_2->z[i];
         gal_3->vel_x[a] = gal_2->vel_x[i];
         gal_3->vel_y[a] = gal_2->vel_y[i];
         gal_3->vel_z[a] = gal_2->vel_z[i];
         ++a;
     }
     // Halo 1
     for (i = gal_1->num_part[0]; i < gal_1->num_part[0]+gal_1->num_part[1]; ++i) {
         gal_3->mass[a] = gal_1->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_1->x[i];
         gal_3->y[a] = gal_1->y[i];
         gal_3->z[a] = gal_1->z[i];
         gal_3->vel_x[a] = gal_1->vel_x[i];
         gal_3->vel_y[a] = gal_1->vel_y[i];
         gal_3->vel_z[a] = gal_1->vel_z[i];
         ++a;
     }
     // Halo 2
     for (i = gal_2->num_part[0]; i < gal_2->num_part[0]+gal_2->num_part[1]; ++i) {
         gal_3->mass[a] = gal_2->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_2->x[i];
         gal_3->y[a] = gal_2->y[i];
         gal_3->z[a] = gal_2->z[i];
         gal_3->vel_x[a] = gal_2->vel_x[i];
         gal_3->vel_y[a] = gal_2->vel_y[i];
         gal_3->vel_z[a] = gal_2->vel_z[i];
         ++a;
     }
     //Bulge 1
 for (i = gal_1->num_part[0]+gal_1->num_part[1]; i < gal_1->num_part[0]+gal_1->num_part[1]+gal_1->num_part[2]; ++i) {
         gal_3->mass[a] = gal_1->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_1->x[i];
         gal_3->y[a] = gal_1->y[i];
         gal_3->z[a] = gal_1->z[i];
         gal_3->vel_x[a] = gal_1->vel_x[i];
         gal_3->vel_y[a] = gal_1->vel_y[i];
         gal_3->vel_z[a] = gal_1->vel_z[i];
         ++a;
     }
     // Bulge 2
  for (i = gal_2->num_part[0]+gal_2->num_part[1]; i < gal_2->num_part[0]+gal_2->num_part[1]+gal_2->num_part[2]; ++i) {
         gal_3->mass[a] = gal_2->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_2->x[i];
         gal_3->y[a] = gal_2->y[i];
         gal_3->z[a] = gal_2->z[i];
         gal_3->vel_x[a] = gal_2->vel_x[i];
         gal_3->vel_y[a] = gal_2->vel_y[i];
         gal_3->vel_z[a] = gal_2->vel_z[i];
         ++a;
     }

    return 0;
}

