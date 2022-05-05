/*-----------------------------------------------------------------------------
/
/ Filename: starscream.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: Starscream creates galaxies. 
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
/ Date: 2008/04/15
//Adaptado Juan Carlos Luna //Marzo 2013
*///---------------------------------------------------------------------------

/* The following program is an example of how to use Starscream to create initial
   conditions for a galaxy.                                                           */

// Be sure to include the starscream header file, wherever it is located.
//
// MAR :: NagBody
//#define global
//

//#include "../../../General_libs/general/stdinc.h"
#include "globaldefs.h"
#include "protodefs.h"


void MainLoop (void) {

    int j;

     // Galaxy pointers
     galaxy *galaxy_1,*galaxy_2,*galaxy_3;

     if (!(galaxy_1=calloc(1,sizeof(galaxy)))) {
        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
         exit(0);
     }
     if (!(galaxy_2=calloc(1,sizeof(galaxy)))) {
        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
         exit(0);
     }
     if (!(galaxy_3=calloc(1,sizeof(galaxy)))) {
        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
         exit(0);
     }

     //Carga snapshot y transfiere la información en una galaxia de formato de Starscream en Disco, Halo y Bulbo.
    load_gadget2_galaxy(cmd.ingal1,galaxy_1);
     galaxy_1->total_mass=0;
 for(j = 0; j < galaxy_1->num_part[0]+galaxy_1->num_part[1]+galaxy_1->num_part[2]; ++j)
{
    galaxy_1->total_mass+=galaxy_1->mass[j];
}

    // Create a second galaxy by copying the first.
     copy_galaxy(galaxy_1,galaxy_2,0);

printf("\n Afuera de copy galaxy ... [1]\n");

    galaxy_2->total_mass=0;
 for(j = 0; j < galaxy_2->num_part[0]+galaxy_2->num_part[1]+galaxy_2->num_part[2]; ++j)
{
    galaxy_2->total_mass+=galaxy_2->mass[j];
}

printf("\n Afuera de copy galaxy ... [2]\n");

    // This function sets two galaxies on a parabolic orbit.
//    set_orbit_parabolic(galaxy_1,galaxy_2,1.492945497,0.1);
    set_orbit_parabolic(galaxy_1,galaxy_2,cmd.radius,cmd.p);

printf("\n Afuera de set orbit parabolic ...\n");
    
    // Now that the galaxies have been put on an orbit with each other,
     // store them in the same galaxy object.
     create_galaxy_system(galaxy_1,galaxy_2,galaxy_3);

printf("\n Afuera de create galaxy system ...\n");
    destroy_galaxy(galaxy_1,0);
printf("\n Afuera de destroy galaxy 1...\n");
    destroy_galaxy(galaxy_2,0);
printf("\n Afuera de destroy galaxy 2...\n");

    // Write initial conditions for the massively parallel N-body code Gadget2.
    write_gadget_ics(galaxy_3,cmd.snapoutfile);

printf("\n Afuera de write gadget ics ...\n");

    // Destroy a galaxy. If the galaxy can not be destroyed, return an error. This
     // function will SEGFAULT if the arrays in the galaxy can not be freed.
//     destroy_galaxy(galaxy_1,0);
//     destroy_galaxy(galaxy_2,0);
     destroy_galaxy_system(galaxy_3,0);

     // Destroy the random number environment.
     gsl_rng_free(r);

     // Print an exit statement and call it quits.
     fprintf(stderr,"\nDone.\n");
//     return 0;
}

void EndRun(void)
{
    char   buf[200];
    FILE *fd;
    
    fclose(gd.outlog);
    printf("\nFinal CPU time : %lf\n\n", cputime() - gd.cpuinit);
}

