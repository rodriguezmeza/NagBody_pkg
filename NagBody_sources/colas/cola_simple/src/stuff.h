/*
    Copyright (c) 2011-2013       Svetlin Tassev
                           Harvard University, Princeton University
 
    This file is part of COLAcode.

    COLAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    COLAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with COLAcode.  If not, see <http://www.gnu.org/licenses/>.
*/

// Mar
#ifndef _stuff_h
#define _stuff_h

#if !defined(global)                    // global def question must be here
#  define global extern
#endif

//#define local     static
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

#include <fftw3.h>

// Mar
global char input_fname[200];
global float Scale,BoxSize;
global int NROW,NGRID,NumPart;
global int GridScale;
global float Om;
global float nLPT;
global float subtractLPT;
global int fullT;
global int StdDA;
//


struct quadratic_params
       {
         double y;
       };


struct io_header_1
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
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} ;

// Mar
global struct particle_data
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
}*P,*PR,*PZ,*PI,*P2,*POrig,*PInit,*PS;



// Mar
global float *N11;
global float *N12;
global float *N13;

global fftwf_complex *FN11;
global fftwf_complex *FN12;
global fftwf_complex *FN13;
//

// Mar
global fftwf_plan p11,p12,p13;

global fftwf_complex *P3D;
global float * density;

global float * dX2;
global float * dY2;
global float * dZ2;

global float * dX;
global float * dY;
global float * dZ;

global float * dXz;
global float * dYz;
global float * dZz;
//



void PtoMesh(float *densityz, struct particle_data *Pz);
void MtoParticles(struct particle_data *P, float *particles,float*mesh);
float putInBox(float din);
int load_snapshot(char *shortfile,int empty, struct particle_data **Ppointer);
int allocate_memory(int**IdPointer, struct particle_data **Ppointer);
int reordering(int*Id, struct particle_data *P);
float fix(float x,float q);
void obtainDisplacements(struct particle_data *P,float *dX,float *dY,float *dZ);
float meanM(float *arr,int N,float mean);

void forces(fftwf_complex *P3D,int filter,int NGRID);


void fft(float *arr,fftwf_complex *fft,int NROW);



void readSnapshots(void);
void Displacements(void);






int putInGrid(int din);

void WRtNROW(float *d,int i,int j,int k,float f);
float REdN(float *d,int i,int j,int k,int NGRID);




void rearrange(float *dX,float*dY,float*dZ);


float growthD(float a);
float growthDtemp(float a);
float Qfactor(float a);


void ReconstructNormalOrder(float *dXr,float *dYr,float *dZr,struct particle_data * PR);
void slice(struct particle_data *P);
float growthD2temp(float a);
float growthD2(float a);
void velRSD(struct particle_data *P,float A);



double Sq(double ai,double af,double aRef);
double fun (double x, void * params);
double gpQ(double a);
double DERgpQ(double a);
double Sphi(double ai,double af,double aRef);


float decayD(float a);
double DprimeQ(double a,float nGrowth);
double DprimeQprimeQ(double a,int nGrowth);
float growthD1Sq(float a);
float growthD2v(float a);



double   SphiStd(double ai,double af);
double   SqStd(double ai,double af);
float AofTime(float y) ;
double   CosmoTime(double af);

#endif // ! _stuff_h
