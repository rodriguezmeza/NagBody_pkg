/* =============================================================================
	MODULE: galmod.c		[galaxy_models]
	Written by: M.A. Rodriguez-Meza
	Starting date:	May 2006
	Purpose:
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: June 11, 2008;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved.
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/


#include "globaldefs.h"
#include "protodefs.h"

local void indmass(void);
local void setdisk(void);
local void cmdisk(void);
local void surfden(void);


// ---------------------------------------------


void InitDisk(void)
{
/*
!     Subroutine to initialize density structure of the disk.  Disks are
!     modeled as isothermal sheets (Spitzer 1942, Ap.J. 95, 329):
!
!        rho[D](R,z) = rho[D](0,0) exp(-R/h) [sech(z/z0)]**2 ,
!
!     where rho[D](R,z) is the volume mass density of the disk at
!     cylinderical coordinates (R,z), h is the  exponential scale
!     length of the disk and z0 is the z scale height.
*/
	printf("\nInitDisk...\n");
	indmass();
	setdisk();
	cmdisk();
	surfden();
}


local void indmass(void)
{
// Subroutine to initialize arrays, etc. having to do with disk
// masses.

	int i;
	real dgasmass;

	printf("\nindmass...\n");

	gd.xgasmass=cmd.gasmass-cmd.gasmass*((1.-(1.+cmd.rmingas/gd.h)*
			rexp(-cmd.rmingas/gd.h))/(1.-(1.+cmd.rmaxgas/gd.h)*
			rexp(-cmd.rmaxgas/gd.h)));
	dgasmass=cmd.gasmass-gd.xgasmass;
	if (dgasmass < 0.0) dgasmass=0.0;
	cmd.gasmass=cmd.gasmass-dgasmass;

	gd.xgasmass=cmd.gasmass;

	if (!cmd.selfggas) cmd.gasmass=0.0;

	for (i=cmd.ndgas+1; i<=gd.nbodies; i++)
           pmass[i]=(gd.diskmass-cmd.gasmass)/cmd.ndstars;

	for (i=1; i<=cmd.ndgas; i++)
           pmass[i]=cmd.gasmass/cmd.ndgas;
}


local void setdisk(void)
{
// Establish an exponental surface density distribution with scale 
// length h and constant scale height z0.

	int n,m;
	real phi,p,r,r2,anum,zp,zps,bnum,pz;

	printf("\nsetdisk...");

	n=0;
 
NDIEZ:
 
	r2=xrandom(0.0,1.0)*cmd.rmaxgas*cmd.rmaxgas;

	if (r2 < rsqr(cmd.rmingas))
		goto NDIEZ;
 
	phi=xrandom(0.0,1.0)*2.*PI;
	r=rsqrt(r2);
	p=rexp(-r/gd.h);
	anum=xrandom(0.0,1.0);

	if (anum > p)
		goto NDIEZ;

	n=n+1;
 
	x[n]=r*rcos(phi);
	y[n]=r*rsin(phi);
 
	if (n < cmd.ndgas) 
		goto NDIEZ;

NQUINCE:
 
	r2=xrandom(0.0,1.0)*cmd.rmax*cmd.rmax;
	phi=xrandom(0.0,1.0)*2.*PI;
	r=rsqrt(r2);
	p=rexp(-r/gd.h);
	anum=xrandom(0.0,1.0);

	if (anum > p)
		goto NQUINCE;
 
	n=n+1;
	x[n]=r*rcos(phi);
	y[n]=r*rsin(phi);
 
	if (n < gd.ndisk)
		goto NQUINCE;
 
	m=0;

NVEINTE:
 
	if (m < cmd.ndgas) {
		zp=2.*(xrandom(0.0,1.0)-0.5)*cmd.zmaxgas;
		zps=zp/cmd.z0gas;
     } else {
		zp=2.*(xrandom(0.0,1.0)-0.5)*cmd.zmax;
		zps=zp/cmd.z0;
	}

	bnum=xrandom(0.0,1.0);

	pz=4./((rexp(zps)+rexp(-zps))*(rexp(zps)+rexp(-zps)));

	if (bnum > pz)
		goto NVEINTE;

	m=m+1;
	z[m]=zp;

	if (m < gd.ndisk) 
		goto NVEINTE;

	printf("\nsetdisk:  Disk established <<\n");
}

local void cmdisk(void)
{
//  Subroutine to transform disk coordinates to center of mass.

	int i,j;
	real xsum,ysum,zsum,psum,xcm,ycm,zcm;

	printf("\ncmdisk...\n");

	psum=0.0;
	xsum=0.0;
	ysum=0.0;
	zsum=0.0;

	for (i=1; i<=gd.ndisk; i++) {
		psum=psum+pmass[i];
		xsum=xsum+pmass[i]*x[i];
		ysum=ysum+pmass[i]*y[i];
		zsum=zsum+pmass[i]*z[i];
	}

	xcm=xsum/psum;
	ycm=ysum/psum;
	zcm=zsum/psum;

	for (j=1; j<=gd.ndisk; j++) {
		x[j]=x[j]-xcm;
		y[j]=y[j]-ycm;
		z[j]=z[j]-zcm;
		radcyl[j]=rsqrt(x[j]*x[j]+y[j]*y[j]);
		radsph[j]=rsqrt(radcyl[j]*radcyl[j]+z[j]*z[j]);
	}

	printf("\ncmdisk  CMDISK completed <<\n");
	printf("\ncmdisk  Corrections : x y z : %g %g %g\n",xcm,ycm,zcm);
}

local void surfden(void)
{
    // Subroutine to compute the surface density corresponding to each
    // particle.
    
    int i;
    real corr;
    
    printf("\nsurfden...\n");
    
    corr=1.-rexp(-cmd.rmax/gd.h)*(1.+cmd.rmax/gd.h);
    
    gd.surfd0=gd.diskmass/(2.*PI*gd.h*gd.h*corr);
    
    for (i=1; i<=gd.ndisk; i++)
        surfd[i]=gd.surfd0*rexp(-radcyl[i]/gd.h);
    
    printf("\nsurfden  Disk surface density computed <<\n");
}
