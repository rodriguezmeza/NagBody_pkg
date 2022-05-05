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

#include "../../../General_libs/general/stdinc.h"
#include "../../../General_libs/general/constant.h"
#include "../../../General_libs/math/mathfns.h"
#include "../../../General_libs/math/numrec.h"

#include "../../../General_libs/math/vectmath.h"
#include "../../../General_libs/NagBody/nagbody.h"

#include "globaldefs.h"
#include "protodefs.h"

local void inbmass(void);
local void setbulge(void);
local void cmbulge(void);
local void indmass(void);
local void setdisk(void);
local void cmdisk(void);
local void surfden(void);
local void inhmass(void);
local void sethalo(void);
local void cmhalo(void);

local void bulgepot(void);

local void sigmap(void);
local void sigmar(void);
local void sigmaz(void);

local void circv(void);

local void radacc(void);
local void insmass(void);
local void satvel(void);

local void cmsat(void);
local void cmsatmod(void);

local void setsat(void);

local void forceb(void);
local void forcebs(void);
local void forced(void);
local void forced_direct_original(void);
local void forced_direct_new(void);
local void forced_tree(void);
local void forceds(void);
local void forceh(void);
local void forcehs(void);

local void setmesh(int nphi, realptr, realptr, realptr);
local void setrot(void);

local void setsigma(void);
local void sigalar(void);
local void sigcheck(void);

local void meanrot(void);
local void meshacc(int nphi, realptr, realptr, realptr, realptr);

local void dadr(void);
local void getkappa(void);
local void halopot(void);

// Structure data and routines to compute forces and potential

#define NINTERP 30000

static acsmooth_st[NINTERP+1];
static phsmooth_st[NINTERP+1];

local void set_interpolation_force_table(void);

// ---------------------------------------------

void InitBulge(void)
{
/*
     Subroutine to initialize density structure of the bulge.  Bulges
     are modeled by the density profile analyzed by Hernquist (1990;
     Ap. J. 356, 359).  Namely,

        rho[B](r) = M_bulge * a_bulge /(r * (a_bulge + r)**3),

     where rho[B](r) is the volume mass density of the bulge at
     spherical coordinate (r), a_bulge is the scale-length of the 
     bulge and M_bulge is the bulge mass.
*/
	printf("\nInit_Bulge...\n");
	if (cmd.usebulge && cmd.selfgbul) {
		inbmass();
		setbulge();
		cmbulge();
	}
}

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

void InitHalo(void)
{
/*
     Subroutine to initialize density structure of the halo.  Halos are
     modeled either according to the density profile

        rho[H](r)   EXP[-r*r/(r_c*r_c)]
       --------- = -------------------
        rho[H](0)   (1 +  [r/gamma]**2)
 
     where rho[H](r) is the volume mass density of the halo at
     spherical radius r, gamma is the halo mass scale length,
     and r_c is a cut-off radius; or according to the model analyzed
     by Hernquist (1990; Ap. J. 356; 359).  Namely,

        rho[H](r) = M_halo * a_halo /(r * (a_halo + r)**3),

     where rho[H](r) is the volume mass density of the halo at
     spherical coordinate (r), a_halo is the scale-length of the
     halo and M_halo is the halo mass.
*/

	int i, j;
	real drhalo,qhalo,alphhalo,xj,deltax;				// iso halo

	printf("\nInitHalo...\n");

	if (!cmd.usehalo) {
	} else {
		if (scanopt(cmd.halotype, "IS")) {
			printf("\nHalo type is IS\n");
			gd.ntabhalo=maxtabh;
			drhalo=cmd.rmaxhalo/((real)(gd.ntabhalo-1));

			qhalo=cmd.gamhalo/cmd.rthalo;
			alphhalo=1./(1.-rsqrt(PI)*qhalo*rexp(qhalo*qhalo)*erfcc(qhalo));
 
			rhalo[1]=1.e-10;

            for (i=2; i<=gd.ntabhalo; i++)
                 rhalo[i]=((real)(i-1))*drhalo;
 
			xmhalo[1]=0.0;

			for (i=2; i<=gd.ntabhalo; i++) {
                 xmhalo[i]=xmhalo[i-1];

				for (j=1; j<=100; j++) {
					deltax=(rhalo[i]-rhalo[i-1])/(100.0*cmd.rthalo);
                    xj=rhalo[i-1]/cmd.rthalo+j*deltax-0.5*deltax;
                    xmhalo[i]=xmhalo[i]+rsqr(xj)*rexp(-rsqr(xj))*deltax/
							(rsqr(xj)+rsqr(qhalo));
				}
			}

			for (i=1; i<=gd.ntabhalo; i++) 
                 xmhalo[i]=xmhalo[i]*2.0*cmd.halomass*alphhalo/SQRTPI;

			for (i=1; i<=gd.ntabhalo; i++)
                 uhalo[i]=-xmhalo[i]/rhalo[i]-cmd.halomass*alphhalo*
						rexp(rsqr(qhalo))*expint(rsqr(rhalo[i])/rsqr(cmd.rthalo)
						+ rsqr(qhalo))/(SQRTPI*cmd.rthalo);

			if (rabs(xmhalo[gd.ntabhalo]-cmd.halomass)/cmd.halomass > 1.e-3)
				error("\nHalo mass error in readhalo\n");

		} else {
			printf("\nHalo type must be LH\n");
			if (!scanopt(cmd.halotype, "LH")) {
				error("halo: Halo type error\n");
			} else { 
			}
		}
	}

	if (cmd.usehalo && cmd.selfghal) {
		inhmass();
		sethalo();
		cmhalo();
	}
}

void StackMod(void)
{
//	Init parameters to optionally stack together two DBH models on a
//  parabolic orbit.  By convention, the orbital plane is the z = 0 plane.

	int i;
	real rmod1, rmod2, mtotal, vsep;

	printf("\nStackMod...\n");

	if (cmd.usesat) {
		cmd.addmods=FALSE;
		printf("\nIf usesat is TRUE addmods is FALSE\n");
		return;
	}

	if (cmd.usehalo && !cmd.selfghal) {
		printf("\nSorry, self-gravitating halo required\n");
		cmd.addmods=FALSE;
		return;
	}

	if (cmd.usebulge && !cmd.selfgbul) {
		printf("\nSorry, self-gravitating bulge required\n");
		cmd.addmods=FALSE;
		return;
	}

	rmod1=cmd.rsep/2.0;
	rmod2=cmd.rsep/2.0;

	mtotal=0.0;

	for (i=1; i<=gd.nbodies; i++)
           mtotal=mtotal+pmass[i];

	gd.xmod1=rmod1-cmd.rp;
	gd.ymod1=rsqrt(2.0*rmod1*cmd.rp-cmd.rp*cmd.rp);
	gd.zmod1=0.0;

	gd.vxmod1= -rsqrt(mtotal*(2.0*rmod1-cmd.rp))/(2.0*rmod1);
	gd.vymod1= -rsqrt(mtotal*cmd.rp)/(2.0*rmod1);
	gd.vzmod1=0.0;

	gd.xmod2=-gd.xmod1;
	gd.ymod2=-gd.ymod1;
	gd.zmod2=0.0;

	gd.vxmod2= -gd.vxmod1;
	gd.vymod2= -gd.vymod1;
	gd.vzmod2=0.0;

	vsep=rsqrt(rsqr(gd.vxmod1-gd.vxmod2)+rsqr(gd.vymod1-gd.vymod2)
			+rsqr(gd.vzmod1-gd.vzmod2));

	printf("\nvsep = %g\n",vsep);
}


void AddSat(void)
{
	real menclose,fx,fpx,xold,xnew,rhobar,acoef,rtidal;			// add sat
	int i;

	printf("\nAddSat...\n");

	if (cmd.usesat) {
		gd.radsat=rsqrt(cmd.xsat*cmd.xsat+cmd.ysat*cmd.ysat+cmd.zsat*cmd.zsat);

		menclose=0.0;

		for (i=1; i<=gd.nbodies-cmd.nsat; i++)
			if (radsph[i] <= gd.radsat) 
				menclose=menclose+pmass[i];

		if (cmd.usebulge && (!cmd.selfgbul)) 
			menclose=menclose+cmd.bulgmass*rsqr(gd.radsat)
					/rsqr(gd.radsat+cmd.abulge);

		if (cmd.usehalo && (!cmd.selfghal)) {

			if (scanopt(cmd.halotype, "LH"))
				menclose=menclose+cmd.halomass*rsqr(gd.radsat)
						/rsqr(gd.radsat+cmd.ahalo);

			if (scanopt(cmd.halotype, "IS")) {
 
				for (i=2; i<=gd.ntabhalo; i++)
					if (gd.radsat >= rhalo[i-1] && gd.radsat < rhalo[i])
						menclose=menclose+xmhalo[i-1];

				if (gd.radsat >= rhalo[gd.ntabhalo]) 
					menclose=menclose+xmhalo[gd.ntabhalo];
 			}
		}

		rhobar=3.0*menclose/(4.0*PI*rqbe(gd.radsat));
		acoef=9.0*cmd.satmass/(4.0*PI*rqbe(cmd.asat)*rhobar);
		printf("\n%s : %g %g %g %g %g %g\n\n",
				"rhobar, menclose, radsat, satmass, asat, acoef",
				rhobar, menclose, gd.radsat, cmd.satmass, cmd.asat, acoef);

		xnew=1.0;

SETENTA:
		xold=xnew;
		fx=rqbe(xold)+2.0*rsqr(xold)+xold-acoef;
		fpx=3.0*rsqr(xold)+4.0*xold+1.0;

		xnew=xold-fx/fpx;

		printf("\nfx, fpx, xnew : %g %g %g %g",fx, fpx, xnew, acoef);
 
		if (rabs((xold-xnew)/xold) > 1.e-6) 
			goto SETENTA;
 
		rtidal=xnew*cmd.asat;

		printf("\nrtidal = %g",rtidal);

		insmass();
		setsat();
		satvel();
		cmsat();
		cmsatmod();
	}
}

local void inbmass(void)
{
/*
     Subroutine to initialize arrays, etc. having to do with bulge
     masses.
*/
	int i;
	real mtrunc;

	printf("\ninbmass...");
	for (i=gd.nbodies+1-cmd.nbulge-cmd.nhalo-cmd.nsat; 
					i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++)
		pmass[i] = cmd.bulgmass/((real) cmd.nbulge);
	if (!cmd.axibulge)
		cmd.bulgmass = cmd.bulgmass
					/(cmd.rmaxbulg*cmd.rmaxbulg/rsqr(cmd.rmaxbulg+cmd.abulge));
	else {
		mtrunc = rsqrt(rsqr(cmd.rmaxbulg)/rsqr(cmd.abulge)
				+rsqr(cmd.zmaxbulg)/rsqr(cmd.cbulge));
		cmd.bulgmass = cmd.bulgmass/(rsqr(mtrunc)/rsqr(1.0+mtrunc));
	}
	printf("\ninbmass corrected bulge mass = %g\n",cmd.bulgmass);
}

local void setbulge(void)
{
/*
   Establish a bulge corresponding to model analyzed by Hernquist
   (1991).   This subroutine initializes only the spatial coordinates;
   velocities are set in bulgevel.f.  The variable const represents 
   the maximum value of r**2 * v**2 * f(q), which occurs at 
   r/r_0 = 0.63798179 and v/v_g = v_circ.
*/
	int ntot;
	real constante,xr,xv,e,p,fq,pn,q,cth,sth,signs,phi,
		cthv,sthv,phiv,Rtest,ztest,mtest,ftest;

	printf("\nsetbulge...");

	ntot=0;
 
	constante=1.20223157581242;
 
	if (!cmd.axibulge) {
		while (ntot<cmd.nbulge) {
			do {
				do {
					xr=xrandom(0.0,1.0)*cmd.rmaxbulg;
					xv=xrandom(0.0,1.0)*rsqrt(2.*cmd.bulgmass/cmd.abulge);
					e=0.5*xv*xv-cmd.bulgmass/(xr+cmd.abulge);
				} while (e<=0 && e>= -cmd.bulgmass/cmd.abulge);

				q=rsqrt(-cmd.abulge*e/cmd.bulgmass);
				p=xrandom(0.0,1.0);
 
				fq=(3.*asin(q)+q*rsqrt(1.-q*q)*(1.-2.*q*q)*(8.*rpow(q,4)-
					8.*rsqr(q)-3.))/rpow(1.-rsqr(q),2.5);
				pn=(xr*xr/(cmd.abulge*cmd.abulge))
					*(xv*xv/(cmd.bulgmass/cmd.abulge))*fq/constante;
 
				if (pn>1.0) {
					fprintf(gd.outlog,
							"\nsetbulge: pn=%g, error in setbulge\n",pn);
					pn=1.0;
				}
			} while (p>pn);

			cth=2.*(xrandom(0.0,1.0)-.5);
			sth=rsqrt(1.-cth*cth);
			signs=2.*(xrandom(0.0,1.0)-.5);
			cth=signs*cth/rabs(signs);
			phi=TWOPI*xrandom(0.0,1.0);
 
			cthv=2.*(xrandom(0.0,1.0)-.5);
			sthv=rsqrt(1.-cthv*cthv);
			signs=2.*(xrandom(0.0,1.0)-.5);
			cthv=signs*cthv/rabs(signs);
			phiv=TWOPI*xrandom(0.0,1.0);
 
			ntot=ntot+1;

			x[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=xr*sth*rcos(phi);
			y[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=xr*sth*rsin(phi);
			z[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=xr*cth;
			vx[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=xv*sthv*rcos(phiv);
			vy[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=xv*sthv*rsin(phiv);
			vz[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=xv*cthv;

		}
	} else {
		while (ntot<cmd.nbulge) {
			do {
				Rtest=cmd.rmaxbulg*xrandom(0.0,1.0);
				ztest=cmd.zmaxbulg*xrandom(0.0,1.0);
				p=xrandom(0.0,1.0);
				if(p<=0.5) ztest = -ztest;
				mtest=rsqrt( rsqr(Rtest)/rsqr(cmd.abulge)
						+rsqr(ztest)/rsqr(cmd.cbulge) );
				ftest=Rtest/(cmd.abulge*cmd.abulge*cmd.cbulge
							*mtest*rpow(1.0+mtest,3.0));
				ftest=ftest*cmd.abulge*cmd.cbulge;
				p=xrandom(0.0,1.0);
			} while (p>ftest);

			ntot=ntot+1;

			z[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=ztest; 
			p=xrandom(0.0,1.0);

			x[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=Rtest*rcos(TWOPI*p);
			y[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=Rtest*rsin(TWOPI*p);

			vx[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=0.0;
			vy[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=0.0;
			vz[ntot+gd.nbodies-cmd.nbulge-cmd.nhalo-cmd.nsat]=0.0;
		}
	}
	printf("\nsetbulge  Bulge established <<\n");
}

local void cmbulge(void)
{
/*
   Subroutine to transform bulge coordinates to center of mass.
*/
	int i,j;
	real xsum,ysum,zsum,psum,xcm,ycm,zcm;

	printf("\ncmbulge...");

	psum=0.0;
	xsum=0.0;
	ysum=0.0;
	zsum=0.0;
        
	for (i=gd.nbodies-cmd.nbulge+1-cmd.nhalo-cmd.nsat; 
					i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++) {
		psum=psum+pmass[i];
		xsum=xsum+pmass[i]*x[i];
		ysum=ysum+pmass[i]*y[i];
		zsum=zsum+pmass[i]*z[i];
	}

	xcm=xsum/psum;
	ycm=ysum/psum;
	zcm=zsum/psum;

	for (j=gd.nbodies-cmd.nbulge+1-cmd.nhalo-cmd.nsat;
					j<=gd.nbodies-cmd.nhalo-cmd.nsat; j++) {
		x[j]=x[j]-xcm;
		y[j]=y[j]-ycm;
		z[j]=z[j]-zcm;
		radcyl[j]=rsqrt(x[j]*x[j]+y[j]*y[j]);
		radsph[j]=rsqrt(radcyl[j]*radcyl[j]+z[j]*z[j]);
	}

	printf("\ncmbulge:  CMBULGE completed <<");
	printf("\ncmbulge:  Corrections : x y z : %g %g %g\n",xcm,ycm,zcm);
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


local void inhmass(void)
{
/*
C     Subroutine to initialize arrays, etc. having to do with halo
C     masses.
*/

	int i;

	printf("\ninhmass...\n");

	for (i=gd.nbodies+1-cmd.nhalo-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++) 
           pmass[i]=cmd.halomass/cmd.nhalo;

	if (scanopt(cmd.halotype, "LH"))
		cmd.halomass=cmd.halomass
					/(cmd.rmaxhalo*cmd.rmaxhalo/rsqr(cmd.rmaxhalo+cmd.ahalo));
	else
		if (cmd.rmaxhalo >= rhalo[gd.ntabhalo])
			cmd.halomass=cmd.halomass*cmd.halomass/xmhalo[gd.ntabhalo];
		else
			for (i=2; i<=gd.ntabhalo; i++)
				if (cmd.rmaxhalo > rhalo[i-1] && cmd.rmaxhalo <= rhalo[i])
					cmd.halomass=cmd.halomass*cmd.halomass/xmhalo[i];

	printf("\ninhmass:  corrected halo mass = %g\n",cmd.halomass);
}


local void sethalo(void)
{
/*
Establish a halo corresponding either to the model analyzed by 
Hernquist (1991) or a truncated, non-singular isothermal sphere.
This subroutine initializes only the spatial coordinates; velocities
are set in halovel.f.
*/

	int ntot;
	real constante,xr,xv,e,p,fq,pn,q,cth,sth,signs,phi,cthv,sthv,phiv;

	printf("\nsethalo...\n");

	if (scanopt(cmd.halotype, "LH")) {

	printf("\nhalo type is LH...\n");

		ntot=0;

		constante=1.20223157581242;

NDIEZ: 
		xr=xrandom(0.0,1.0)*cmd.rmaxhalo;
		xv=xrandom(0.0,1.0)*rsqrt(2.*cmd.halomass/cmd.ahalo);
		e=0.5*xv*xv-cmd.halomass/(xr+cmd.ahalo);

		if (e > 0.0 || e < -cmd.halomass/cmd.ahalo) {
			goto NDIEZ;
		}

		q=rsqrt(-cmd.ahalo*e/cmd.halomass);

		p=xrandom(0.0,1.0);
 
		fq=(3.*asin(q)+q*rsqrt(1.-q*q)*(1.-2.*q*q)*(8.*rpow(q,4)-
			8.*q*q-3.))/rpow(1.-q*q,2.5);
		pn=(xr*xr/(cmd.ahalo*cmd.ahalo))*(xv*xv/(cmd.halomass/cmd.ahalo))*fq/constante;

		if (pn > 1.0) {
			printf("\npn error in sethalo %g",pn);
			pn=1.0;
		}

		if (p <= pn) {
 
			cth=2.*(xrandom(0.0,1.0)-.5);
			sth=rsqrt(1.-cth*cth);
			signs=2.*(xrandom(0.0,1.0)-.5);
			cth=signs*cth/rabs(signs);
			phi=TWOPI*xrandom(0.0,1.0);

			cthv=2.*(xrandom(0.0,1.0)-.5);
			sthv=rsqrt(1.-cthv*cthv);
			signs=2.*(xrandom(0.0,1.0)-.5);
			cthv=signs*cthv/rabs(signs);
			phiv=TWOPI*xrandom(0.0,1.0);

			ntot=ntot+1;

			x[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xr*sth*rcos(phi);
			y[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xr*sth*rsin(phi);
			z[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xr*cth;
			vx[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xv*sthv*rcos(phiv);
			vy[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xv*sthv*rsin(phiv);
			vz[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xv*cthv;
 
		} else {
			goto NDIEZ;
		}

		if (ntot < cmd.nhalo) goto NDIEZ;

	printf("\nend halo type LH...\n");
	} else {
  
		ntot=0;

		constante=rsqrt(0.5*(-rsqr(cmd.gamhalo)+rsqrt(rpow(cmd.gamhalo,4)+4.0*rsqr(cmd.gamhalo)*
                rsqr(cmd.rthalo))));
		constante=rsqr(constante)*rexp(-rsqr(constante)/rsqr(cmd.rthalo))/(1.+rsqr(constante)/
                rsqr(cmd.gamhalo));
 
NTREINTA:
		xr=xrandom(0.0,1.0)*cmd.rmaxhalo;

		p=xrandom(0.0,1.0);

		fq=xr*xr*rexp(-rsqr(xr)/rsqr(cmd.rthalo))/(1.+rsqr(xr)/rsqr(cmd.gamhalo));

		pn=fq/constante;

		if (pn > 1.0) {
			printf("\npn error in sethalo %g\n",pn);
			pn=1.0;
		}

		if (p <= pn) {

			cth=2.*(xrandom(0.0,1.0)-.5);
			sth=rsqrt(1.-cth*cth);
			signs=2.*(xrandom(0.0,1.0)-.5);
			cth=signs*cth/rabs(signs);
			phi=TWOPI*rand();

			ntot=ntot+1;
 
			x[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xr*sth*rcos(phi);
			y[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xr*sth*rsin(phi);
			z[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=xr*cth;
			vx[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=0.0;
			vy[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=0.0;
			vz[ntot+gd.nbodies-cmd.nhalo-cmd.nsat]=0.0;
 
		} else
			goto NTREINTA;
 
		if (ntot < cmd.nhalo) goto NTREINTA;

	}
 
	printf("\nsethalo:  Halo established <<\n");
}

local void cmhalo(void)
{
//   Subroutine to transform halo coordinates to center of mass.

	int i,j;
	real xsum,ysum,zsum,psum,xcm,ycm,zcm;

	printf("\ncmhalo...");

	psum=0.0;
	xsum=0.0;
	ysum=0.0;
	zsum=0.0;
        
	for (i=gd.nbodies-cmd.nhalo+1-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++) {
		psum=psum+pmass[i];
		xsum=xsum+pmass[i]*x[i];
		ysum=ysum+pmass[i]*y[i];
		zsum=zsum+pmass[i]*z[i];
	}

	xcm=xsum/psum;
	ycm=ysum/psum;
	zcm=zsum/psum;

	for (j=gd.nbodies-cmd.nhalo+1-cmd.nsat; j<=gd.nbodies-cmd.nsat; j++) {
		x[j]=x[j]-xcm;
		y[j]=y[j]-ycm;
		z[j]=z[j]-zcm;
		radcyl[j]=rsqrt(x[j]*x[j]+y[j]*y[j]);
		radsph[j]=rsqrt(radcyl[j]*radcyl[j]+z[j]*z[j]);
	}

	printf("\ncmhalo:  CMHALO completed <<");
	printf("\ncmhalo:  Corrections : x y z : %g %g %g\n",xcm,ycm,zcm);
}

#define maxnbin		10000

void BulgeVel(void)
{
/*
Subroutine to initialize velocities of bulge particles, using the
spherical Jeans equations.  That is, the radial velocity dispersion
at radius r is:

                             infinity
							 /
					   1     |  G M(r)
	<v_r ^2 (r)> =  -------  |  ------ rho_h(r) dr
				    rho_h(r) /   r^2
					 		 r


where rho_h(r) is the bulge density at r.
*/
 	int nbin,i,j,ir,ntest,irlower,irupper,irindex,nturn,ntot,nxsimp;
	int irtmp;

	real radmax,ria,ri,rj,r,dr,rhoi,rhoj,rja,xv,gspeed,
		p,cth,sth,signs,phi,epsilon,
		rlower,rupper,tmass,fracmass,vescape,sigfix,
		vescsig,vphib,radcylb,zb,sigzb,sigpb,sigrcylb,
		vzb,vrcylb,vbtot,axbulge,cxbulge,bxmass;

	realptr xintbin, xmbin, gridmass;

    double cpustart;
    cpustart = cputime();

	xintbin=nr_dvector(0,maxnbin+1000);
	xmbin=nr_dvector(0,maxnbin+1000);
	gridmass=nr_dvector(0,maxnbin+1000);

	printf("\nBulgeVel...\n");
 
	nbin=1000;

	if (nbin>=maxnbin) error("galmod : 001 : nbin error in bulgevel");

	if (!cmd.axibulge)
		radmax=1.02*cmd.rmaxbulg;
	else
		radmax=1.02*rsqrt(rsqr(cmd.rmaxbulg)+rsqr(cmd.zmaxbulg));

	dr=radmax/nbin;
 
	for (i=0; i<=nbin; i++) {
		xintbin[i]=0.0;
		xmbin[i]=0.0;
		gridmass[i]=0.0;
	}

	for (i=1; i<=gd.nbodies-cmd.nsat; i++) {		// SATELLITE NOT INCLUDED...
		r=rsqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
		epsilon=cmd.epsdisk;
		if (i>gd.ndisk && i<gd.nbodies-cmd.nhalo+1-cmd.nsat) 
			epsilon=cmd.epsbulge;

		if (i>gd.nbodies-cmd.nhalo-cmd.nsat) epsilon=cmd.epshalo;
 
		if (cmd.axibulge && i>gd.ndisk && i<gd.nbodies-cmd.nhalo+1) continue;

		rlower=r-epsilon/2.0;
		rupper=r+epsilon/2.0;
		ir=(int) (r/dr);
		ir=ir+1;
		irlower=(int) (rlower/dr);
		irlower=irlower+1;
 
		if (rlower<0.0) {
			irlower=(int) (rabs(rlower)/dr);
			irlower= -irlower;
		}

		irupper=(int) (rupper/dr);
		irupper=irupper+1;
		tmass=0.0;

		for (j=irlower; j<=irupper; j++) {
			irindex=j;
			if (j<ir) irindex=ir;
			if (j>maxnbin+1000) irindex=maxnbin+1000;

			if (irlower==irupper) {
				tmass=tmass+pmass[i];
				gridmass[irindex]=gridmass[irindex]+pmass[i];
			} else {
				if (j!=irlower && j!=irupper) {
					fracmass=(dr/epsilon)*pmass[i];
                    tmass=tmass+fracmass;
                    gridmass[irindex]=gridmass[irindex]+fracmass;
				} else {
					if (j==irupper) {
						fracmass=((rupper-(irupper-1)*dr)/epsilon)*pmass[i];
						tmass=tmass+fracmass;
						gridmass[irindex]=gridmass[irindex]+fracmass;
					} else {
						fracmass=((irlower*dr-rlower)/epsilon)*pmass[i];
						tmass=tmass+fracmass;
						gridmass[irindex]=gridmass[irindex]+fracmass;
                    }
				}
			}
 
		}

		if (rabs(tmass-pmass[i])/pmass[i] > 1.e-4)
			error("\ngalmod : 002 : mass assignment error in bulgevel");
	}

	for (i=1; i<=nbin; i++) {
           xmbin[i]=xmbin[i-1]+gridmass[i];
	}
 
	xmbin[nbin+1]=xmbin[nbin];

	for (i=1; i<=nbin; i++) { 
		ri=i*dr;
		ria=ri/cmd.abulge;
		rhoi=cmd.bulgmass/(2.*PI*rpow(cmd.abulge,3.0)*ria*rpow(1.+ria,3.0));

		for (j=i; j<=nbin; j++) { 
			rj=j*dr+0.5*dr;
			rja=rj/cmd.abulge;
			rhoj=cmd.bulgmass/(2.*PI*rpow(cmd.abulge,3.0)*rja*rpow(1.+rja,3.0));
			xintbin[i]=xintbin[i]+rhoj*0.5*(xmbin[j]+xmbin[j+1])*dr/(rj*rj);
		}

		xintbin[i]=xintbin[i]/rhoi;
	}

	xintbin[nbin+1]=0.0;

	for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++) {
		r=rsqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
		ir=(int) (r/dr);
		irtmp=MIN(ir,nbin); ir=irtmp;
		sigr[i]=((xintbin[ir+1]-xintbin[ir])*r+xintbin[ir]*(ir+1)*dr-xintbin[ir+1]*ir*dr)/dr;
		sigr[i]=rsqrt(sigr[i]);
	}

//   Initialize velocities isotropically so that the distribution of
//   speeds is proportional to v^2 EXP[-v^2/(2*sigma_r^2)].  Limit
//   speed to the local escape speed.
 
	bulgepot();
 
	ntest=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat;
 
	if (!cmd.axibulge) {
		printf("\nnot axibulge ...\n");

CIENTOSESENTA:

		vescape=rsqrt(2.0*rabs(pot[ntest]));
		xv=xrandom(0.0,1.0)*vescape;
 
		vescsig=vescape/(rsqrt(2.0)*sigr[ntest]);
		sigfix=1.0-erfcc(vescsig)
			-8.0*vescsig*(0.75+0.5*rsqr(vescsig))*rexp(-rsqr(vescsig))/(3.0*rsqrt(PI));
		sigfix=sigr[ntest]/rsqrt(sigfix);

		gspeed=0.5*xv*xv*rexp(1.0-xv*xv/(2.0*rsqr(sigfix)))/rsqr(sigfix);

		p=xrandom(0.0,1.0);

		if (p<=gspeed) {
              cth=2.0*(xrandom(0.0,1.0)-0.5);
              sth=rsqrt(1.0-cth*cth);
              signs=2.0*(xrandom(0.0,1.0)-0.5);
              cth=signs*cth/rabs(signs);
              phi=TWOPI*xrandom(0.0,1.0);
              vx[ntest]=xv*sth*rcos(phi);
              vy[ntest]=xv*sth*rsin(phi);
              vz[ntest]=xv*cth;
              ntest=ntest+1;
		} else {
		goto CIENTOSESENTA;
		}
 
		if (ntest<=gd.nbodies-cmd.nhalo-cmd.nsat) goto CIENTOSESENTA;

	} else {
		printf("\naxibulge ...\n");

		axbulge=cmd.abulge;
		cxbulge=cmd.cbulge;
		bxmass=cmd.bulgmass;
		nxsimp=cmd.nsimpson;

		for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++) {
			radcylb=rsqrt(x[i]*x[i]+y[i]*y[i]);
			zb=z[i];

			obsigma(radcylb,zb, &sigrcylb, &sigzb, &sigpb, axbulge,
					cxbulge,bxmass,nxsimp,cmd.rmaxbulg,cmd.zmaxbulg);

			sigrcylb=rsqrt(rsqr(sigrcylb)+rsqr(sigr[i]));
			sigpb=rsqrt(rsqr(sigpb)+rsqr(sigr[i]));
			sigzb=rsqrt(rsqr(sigzb)+rsqr(sigr[i]));

CIENTOSETENTAYCINCO:        

			vrcylb=grandom(ZERO,sigrcylb);
			vphib=grandom(ZERO,sigpb);
			vzb=grandom(ZERO,sigzb);

			vx[i]=(vrcylb*x[i]-vphib*y[i])/radcylb;
			vy[i]=(vrcylb*y[i]+vphib*x[i])/radcylb;
			vz[i]=vzb;
			vbtot=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
			if (vbtot>=2.0*rabs(pot[i])) goto CIENTOSETENTAYCINCO;
		}
 
	}
 
	if (cmd.bulgerot) {
		ntot=0;
		nturn=0.5*cmd.nbulge*(2.0*cmd.brotfrac-1.0);
		if (cmd.brotfrac==0.5) nturn=0;
		if (cmd.brotfrac==1.0) nturn=cmd.nbulge;

		for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++) {
			radcylb=rsqrt(x[i]*x[i]+y[i]*y[i]);
			vphib=(x[i]*vy[i]-y[i]*vx[i])/radcylb;
			if (vphib<0.0 && ntot<nturn) {
				ntot=ntot+1;
				vx[i] = -vx[i];
				vy[i] = -vy[i];
			}
		}
 
	}

	free_dvector(xintbin,0,maxnbin+1000);
	free_dvector(xmbin,0,maxnbin+1000);
	free_dvector(gridmass,0,maxnbin+1000);

	fprintf(stdout,"\n\nBulgeVel CPU time: %g\n",cputime()-cpustart);
	fflush(stdout);
}

#undef maxnbin

local void bulgepot(void)
{
//  Subroutine to compute potentials of bulge particles.

	int ninterp,j,smindex,i,ihalo;
	real deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,drdeldrg,
		drsm,drdotdr,phsm,drhalo;
	realptr phsmooth;

    double cpustart;
    cpustart = cputime();

	phsmooth=nr_dvector(0,30001);

	ninterp=30000;

	deldrg=2./((real)ninterp);
 
	for (i=0; i<=1+ninterp; i++) {
		xw=i*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0)
			phsmooth[i]=-2.*xw3*(ONE/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.;
		else
			phsmooth[i]=-ONE/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.);

		if (xw >= 2.0)
			phsmooth[i]=ONE;
	}
 
	for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++)
           pot[i]=0.0;
 

//   Compute bulge-disk interaction.
//   ------------------------------

	for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++) {	// DO 50 
		for (j=1; j<=gd.ndisk; j++) {
			aa=x[i]-x[j];
			bb=y[i]-y[j];
			cc=z[i]-z[j];
			drdotdr=aa*aa+bb*bb+cc*cc;
			sdrdotdr=rsqrt(drdotdr);
			rinveff=1./(sdrdotdr+1.e-10);
			drdeldrg=sdrdotdr*((real) ninterp)/(cmd.epsdisk+cmd.epsdisk);
			smindex=((int) drdeldrg);
			if (ninterp < smindex) smindex=ninterp;
			if (1.0 < drdeldrg-smindex)
				drsm=1.0;
			else
				drsm=drdeldrg-smindex;
			phsm=(1.-drsm)*phsmooth[smindex]+drsm*phsmooth[1+smindex];
			rinveff=phsm*rinveff;
			pot[i]=pot[i]-pmass[j]*rinveff;
		}
 
		if (cmd.usebulge && cmd.selfgbul && cmd.axibulge) {
			for (j=gd.ndisk+1; j<=i-1; j++) { // DO 45 
				aa=x[i]-x[j];
				bb=y[i]-y[j];
				cc=z[i]-z[j];
				drdotdr=aa*aa+bb*bb+cc*cc;
				sdrdotdr=rsqrt(drdotdr);
				rinveff=1./(sdrdotdr+1.e-10);
				drdeldrg=sdrdotdr*((real) ninterp)/(cmd.epsbulge+cmd.epsbulge);
				smindex=((int) drdeldrg);
				if (ninterp < smindex) smindex=ninterp;
				if (1.0 < drdeldrg-smindex)
					drsm=1.0;
				else
					drsm=drdeldrg-((real) smindex);
				phsm=(1.-drsm)*phsmooth[smindex]+drsm*phsmooth[1+smindex];
				rinveff=phsm*rinveff;
				pot[i]=pot[i]-pmass[j]*rinveff;
			}
 
			for (j=i+1; j<=gd.ndisk+cmd.nbulge; j++) { //DO 47 J=I+1,ndisk+nbulge
				aa=x[i]-x[j];
				bb=y[i]-y[j];
				cc=z[i]-z[j];
				drdotdr=aa*aa+bb*bb+cc*cc;
				sdrdotdr=rsqrt(drdotdr);
				rinveff=1./(sdrdotdr+1.e-10);
				drdeldrg=sdrdotdr*((real) ninterp)/(cmd.epsbulge+cmd.epsbulge);
				smindex=drdeldrg;
				if (ninterp < smindex) smindex=ninterp;
				if (1.0 < drdeldrg-smindex)
					drsm=1.0;
				else
                    drsm=drdeldrg-smindex;
				phsm=(1.-drsm)*phsmooth[smindex]+drsm*phsmooth[1+smindex];
				rinveff=phsm*rinveff;
				pot[i]=pot[i]-pmass[j]*rinveff;
			}
 
		}
 
	}
 

//   Compute bulge-halo interaction.
//   -------------------------------

	if (scanopt(cmd.halotype, "IS")) {
 
           drhalo=rhalo[5]-rhalo[4];
 
		for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++) {
			ihalo=radsph[i]/drhalo;
			ihalo=ihalo+1;
 
			if (radsph[i] > rhalo[gd.ntabhalo])
                 pot[i]=pot[i]+uhalo[gd.ntabhalo]*rhalo[gd.ntabhalo]/radsph[i];
			else
				if (radsph[i] <= rhalo[2])
                    pot[i]=pot[i]+uhalo[2];
				else
					pot[i]=pot[i]+((radsph[i]-rhalo[ihalo])*
							uhalo[ihalo+1]/drhalo-(radsph[i]-
							rhalo[ihalo+1])*uhalo[ihalo]/drhalo); 
		}
 
	} else
		for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++)
              pot[i]= pot[i]-cmd.halomass/(radsph[i]+cmd.ahalo);


//   Compute bulge self-interaction.
//   -------------------------------

	if (cmd.usebulge && cmd.selfgbul && (!cmd.axibulge))
		for (i=gd.nbodies-cmd.nhalo-cmd.nbulge+1-cmd.nsat; i<=gd.nbodies-cmd.nhalo-cmd.nsat; i++)
              pot[i]= pot[i]-cmd.bulgmass/(radsph[i]+cmd.abulge);

	printf("\nbulgepot:  Bulge potentials computed <<\n");

	free_dvector(phsmooth,0,30001);

	fprintf(stdout,"\n\nbulgepot CPU time: %g\n",cputime()-cpustart);
	fflush(stdout);
}

local void sigmap(void)
{
/*
Subroutine to compute the azimuthal dispersion for each particle
given some functional form for the way the shape of the velocity 
ellipsoid (<theta^2 /<pi^2 ) varies with radius (C(x), x=R/h).
The function C is constrained to be one at the origin and approach 
0.5 at large R where we approach the epicyclic limit of small 
perturbations on circular orbits for a flat rotation curve 
potential.

At present the C function is that given by epicyclic theory:

           C = <theta^2 /<pi^2 = 0.25 KAPPA^2/OMEGA^2

where OMEGA is the angular rotational frequency in the disk and
halo field.
*/

	int i;
	real sigratio;

	for (i=1; i<=gd.ndisk; i++) {
		sigratio=rsqrt(cfunc(rotcirc[i],radcyl[i],kappa[i]));
		if (sigratio > 1.0)			// This condition is reported in meanrot
			sigratio=1.;

		sigphi[i]=sigratio*sigr[i];
	}
 
	printf("\nsigmap:  Phi velocity dispersion computed <<\n");
 }

local void sigmar(void)
{
/*
Subroutine to compute the radial velocity dispersion as a function
of radius. The radial dispersion drops like exp(-R/(2*h)) in 
accordance with the van der Kruit and Searle conjecture and direct 
observations of the disk of the Milky Way by Freeman and Lewis.  
The dispersion is normalized to that required by the Toomre formula
for axisymmetric stability at the solar radius for a given value 
of Q.    

Note : This dispersion profile will be modified near the center of 
		the disk to insure no net radial motion and a real mean 
		motion.
*/

	int i,nsum;
	real rtoll,rminsol,rmaxsol,rsum,sigmean,sigzero;

	rtoll=0.25;

	rminsol=gd.rsolar-rtoll;
	rmaxsol=gd.rsolar+rtoll;

	rsum=0.;
	nsum=0;

	for (i=1; i<=gd.ndisk; i++)
		if (radcyl[i] > rminsol && radcyl[i] <= rmaxsol) {
			nsum=nsum+1;
			rsum=rsum+sigt[i];
		}
  
	sigmean=rsum/((real)(nsum));
	sigzero=sigmean/rexp(-gd.rsolar/(2.0*gd.h));
	gd.sigr0=sigzero;

	printf("\nsigmar:  Reference dispersion : %g\n",sigmean);
	printf("sigmar:  Central Toomre dispersion : %g",sigzero);
 
	for (i=1; i<=gd.ndisk; i++)
		if (i <= cmd.ndgas)
			sigr[i]=0.0;
		else
			sigr[i]=sigzero*rexp(-radcyl[i]/(2.0*gd.h));
 
	printf("\nsigmar:  Radial velocity dispersion computed <<\n");
}

local void sigmaz(void)
{
/*
Subroutine to compute the z velocity dispersion.  For the 
isothermal sheet, the z dispersion is given by :

           sigma[z](z=0) = sqrt [ pi g SIGMA(R) z0 ].

NOTE: The Z dispersion is modified by the same function as the 
		radial dispersion profile so the ratio of Z to R dispersion 
		is everywhere equal.
*/

	int i;
	real term;

	term=PI*cmd.z0;

	for (i=1; i<=gd.ndisk; i++)
		if (i <= cmd.ndgas)
			sigz[i]=0.0;
		else
			sigz[i]=rsqrt(term*surfd[i]);
 
	printf("\nsigmaz:  Z velocity dispersion computed <<\n");
}

local void circv(void)
{
// Subroutine to compute centripital rotation velocities for each
// particle.

	int i;
	real smallnum;

	smallnum=1.e-07;

	for (i=1; i<=gd.ndisk; i++)
		if(radcyl[i] > smallnum)
			if(aradcyl[i] < 0.0)
				rotcirc[i]=rsqrt(rabs(aradcyl[i])*radcyl[i]);
			else {
				printf("\ncircv:  particle: %d had arad= %g\n",i,aradcyl[i]);
				rotcirc[i]=0.;
			}
		else {
			printf("\ncircv:  particle : %d had rad = %g",i,radcyl[i]);
			rotcirc[i]=0.;
		}
 
	printf("\ncircv:  Centripital velocities computed <<\n");
}


local void radacc(void)
{
// Subroutine to compute radial acceleration on each disk particle.

	int i;
	real smallnum;

	smallnum=1.e-07;

	for (i=1; i<=gd.ndisk; i++)
		if (radcyl[i] > smallnum)
              aradcyl[i]=(x[i]*ax[i]+y[i]*ay[i])/radcyl[i];
		else {
			printf("\nradacc:  particle : %d had rad = %g",i,radcyl[i]);
			aradcyl[i]=0.;
		}
 
	printf("\nradacc:  Radial acceleration components computed <<\n"); 
}


local void insmass(void)
{
// Subroutine to initialize arrays, etc. having to do with satellite
// masses.

	int i;

	for (i=gd.nbodies+1-cmd.nsat; i<=gd.nbodies; i++)
           pmass[i]=cmd.satmass/cmd.nsat;

	cmd.satmass=cmd.satmass/(cmd.rmaxsat*cmd.rmaxsat/rsqr(cmd.rmaxsat+cmd.asat));

	printf("\ninsmass:  Corrected satellte mass = %g\n",cmd.satmass);
}


local void satvel(void)
{
//  Subroutine to initialize satellite bulk velocity.

	real smallnum,aradsat,vcircsat,vescsat;

	forceds();
	if (cmd.usebulge) forcebs();
	if (cmd.usehalo) forcehs();

	smallnum=1.e-07;

	if (gd.radsat > smallnum)
		aradsat=(cmd.xsat*gd.axsat+cmd.ysat*gd.aysat+cmd.zsat*gd.azsat)/gd.radsat;
	else {
		fprintf(gd.outlog,"\nsatvel:  satellite had rad = %g\n",gd.radsat);
		aradsat=0.;
	}
 
	vcircsat=rsqrt(gd.radsat*rabs(aradsat));
	vescsat=rsqrt(2.0*rabs(gd.potsat));

	printf("\nSatellite circular velocity = %g",vcircsat);
	printf("\nSatellite escape velocity = %g\n",vescsat);
}

local void cmsat(void)
{
// Subroutine to transform satellite coordinates to center of mass
// and to place satellite at its starting point.

	int i,j;
	real xsum,ysum,zsum,vxsum,vysum,vzsum,psum,xcm,ycm,zcm,vxcm,vycm,vzcm;

	psum=0.0;
	xsum=0.0;
	ysum=0.0;
	zsum=0.0;
	vxsum=0.0;
	vysum=0.0;
	vzsum=0.0;

	for (i=gd.nbodies-cmd.nsat+1; i<=gd.nbodies; i++) { 
		psum=psum+pmass[i];
		xsum=xsum+pmass[i]*x[i];
		ysum=ysum+pmass[i]*y[i];
		zsum=zsum+pmass[i]*z[i];
		vxsum=vxsum+pmass[i]*vx[i];
		vysum=vysum+pmass[i]*vy[i];
		vzsum=vzsum+pmass[i]*vz[i];
	}

	xcm=xsum/psum;
	ycm=ysum/psum;
	zcm=zsum/psum;
	vxcm=vxsum/psum;
	vycm=vysum/psum;
	vzcm=vzsum/psum;

	for (j=gd.nbodies-cmd.nsat+1; j<=gd.nbodies; j++) {
		x[j]=x[j]-xcm+cmd.xsat;
		y[j]=y[j]-ycm+cmd.ysat;
		z[j]=z[j]-zcm+cmd.zsat;
		vx[j]=vx[j]-vxcm+cmd.vxsat;
		vy[j]=vy[j]-vycm+cmd.vysat;
		vz[j]=vz[j]-vzcm+cmd.vzsat;
	}

	printf("\ncmsat:  CMSAT completed <<");
	printf("\ncmsat:  Corrections : x y z vx vy vz : %g %g %g %g %g %g\n",
				xcm,ycm,zcm,vxcm,vycm,vzcm);
}

local void cmsatmod(void)
{
// Subroutine to transform entire model to center-of-mass coordinates,
// including the halo.

	int i,j;
	real xsum,ysum,zsum,vxsum,vysum,vzsum,psum,vxcm,vycm,vzcm,xcm,ycm,zcm;

	psum=0.0;

	if (cmd.usehalo && (!cmd.selfghal)) psum=psum+cmd.halomass;
	if (cmd.usebulge && (!cmd.selfgbul)) psum=psum+cmd.bulgmass;

	xsum=0.0;
	ysum=0.0;
	zsum=0.0;
	vxsum=0.0;
	vysum=0.0;
	vzsum=0.0;

	for (i=1; i<=gd.nbodies; i++) {
		psum=psum+pmass[i];
		xsum=xsum+pmass[i]*x[i];
		ysum=ysum+pmass[i]*y[i];
		zsum=zsum+pmass[i]*z[i];
		vxsum=vxsum+pmass[i]*vx[i];
		vysum=vysum+pmass[i]*vy[i];
		vzsum=vzsum+pmass[i]*vz[i];
	}
 
	xcm=xsum/psum;
	ycm=ysum/psum;
	zcm=zsum/psum;
	vxcm=vxsum/psum;
	vycm=vysum/psum;
	vzcm=vzsum/psum;
 
	for (j=1; j<=gd.nbodies; j++) {
		x[j]=x[j]-xcm;
		y[j]=y[j]-ycm;
		z[j]=z[j]-zcm;
		vx[j]=vx[j]-vxcm;
		vy[j]=vy[j]-vycm;
		vz[j]=vz[j]-vzcm;
	}
 
	printf("\ncmsatmod:  CMSATMOD completed <<");
	printf("\ncmsatmod:  Corrections : x y z vx vy vz : %g %g %g %g %g %g",
				xcm,ycm,zcm,vxcm,vycm,vzcm);
}

void cmtv(void)
{
// Subroutine to center-of-mass transform velocities.

	int i,j;
	real vxsum,vysum,vzsum,psum,vxcm,vycm,vzcm;

	vxsum=0;
	vysum=0;
	vzsum=0;
	psum=0;

	for (i=1; i<=gd.nbodies; i++) {
		vxsum=vxsum+pmass[i]*vx[i];
		vysum=vysum+pmass[i]*vy[i];
		vzsum=vzsum+pmass[i]*vz[i];
		psum=psum+pmass[i];
	}

	vxcm=vxsum/psum;
	vycm=vysum/psum;
	vzcm=vzsum/psum;

	for (j=1; j<=gd.nbodies; j++) { 
		vx[j]=vx[j]-vxcm;
		vy[j]=vy[j]-vycm;
		vz[j]=vz[j]-vzcm;
	}
 
	printf("\ncmtv:  CMTV completed <<");
	printf("\ncmtv:  Corrections : vx vy vz : %g %g %g\n",vxcm,vycm,vzcm);
}

local void setsat(void)
{
/*
Establish a satellite corresponding to model analyzed by Hernquist
(1991).  The variable const represents the maximum value of 
r**2 * v**2 * f(q), which occurs at r/r_0 = 0.63798179 and 
v/v_g = v_circ.
*/

	int ntot;
	real constante,xr,xv,e,p,fq,pn,q,cth,sth,signs,phi,cthv,sthv,phiv;

	ntot=0;

	constante=1.20223157581242;

NDIEZ:
	xr=xrandom(0.0,1.0)*cmd.rmaxsat;
	xv=xrandom(0.0,1.0)*rsqrt(2.*cmd.satmass/cmd.asat);
	e=0.5*xv*xv-cmd.satmass/(xr+cmd.asat);
	if (e > 0.0 || e < -cmd.satmass/cmd.asat) goto NDIEZ;

	q=rsqrt(-cmd.asat*e/cmd.satmass);

	p=xrandom(0.0,1.0);

	fq=(3.*asin(q)+q*rsqrt(1.-q*q)*(1.-2.*q*q)*(8.*rpow(q,4)-8.*q*q-3.))/
		rpow(1.-q*q,2.5);
	pn=(xr*xr/(cmd.asat*cmd.asat))*(xv*xv/(cmd.satmass/cmd.asat))*fq/constante;

	if (pn > 1.0) {
		printf("\npn error in setsat: pn = %g\n",pn);
		pn=1.0;
	}
 
	if (p <= pn) {
		cth=2.*(xrandom(0.0,1.0)-.5);
		sth=rsqrt(1.-cth*cth);
		signs=2.*(xrandom(0.0,1.0)-.5);
		cth=signs*cth/rabs(signs);
		phi=TWOPI*xrandom(0.0,1.0);
		cthv=2.*(xrandom(0.0,1.0)-.5);
		sthv=rsqrt(1.-cthv*cthv);
		signs=2.*(xrandom(0.0,1.0)-.5);
		cthv=signs*cthv/rabs(signs);
		phiv=TWOPI*xrandom(0.0,1.0);
		ntot=ntot+1;
		x[ntot+gd.nbodies-cmd.nsat]=xr*sth*rcos(phi);
		y[ntot+gd.nbodies-cmd.nsat]=xr*sth*rsin(phi);
		z[ntot+gd.nbodies-cmd.nsat]=xr*cth;
		vx[ntot+gd.nbodies-cmd.nsat]=xv*sthv*rcos(phiv);
		vy[ntot+gd.nbodies-cmd.nsat]=xv*sthv*rsin(phiv);
		vz[ntot+gd.nbodies-cmd.nsat]=xv*cthv;
	} else
		goto NDIEZ;

	if (ntot < cmd.nsat) goto NDIEZ;
 
	printf("\nsetsat  Satellite established <<\n");
}

local void forceb(void)
{
// Subroutine to compute forces on disk particles from the bulge.

//	int ninterp;
	int i,j,smindex;
//	real deldrg,xw,xw2,xw3,xw4;
	real aa,bb,cc,sdrdotdr,rinveff,
		r3inveff,drdeldrg,drsm,accsm,drdotdr,arad;
//	realptr acsmooth;

//	acsmooth=nr_dvector(0,30001);
 
//	ninterp=30000;

//	deldrg=2./ninterp;

/*	for (i=0; i<=1+ninterp; i++) {
		xw=i*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0)
			acsmooth[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
		else
			acsmooth[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;

		if (xw >= 2.0)
			acsmooth[i]=1.0;
	}
*/

	if (cmd.usebulge && cmd.selfgbul && cmd.axibulge)
		for (i=1; i<=gd.ndisk; i++)
			for (j=gd.ndisk+1; j<=gd.ndisk+cmd.nbulge; j++) {
				aa=x[i]-x[j];
				bb=y[i]-y[j];
				cc=z[i]-z[j];
				drdotdr=aa*aa+bb*bb+cc*cc;
				sdrdotdr=rsqrt(drdotdr);
				rinveff=1./(sdrdotdr+1.e-10);
				r3inveff=rinveff/(drdotdr+1.e-10);
				drdeldrg=sdrdotdr*NINTERP/(cmd.epsdisk+cmd.epsbulge);
				smindex=drdeldrg;
				if (NINTERP < smindex) smindex=NINTERP;
				if (1.0 < drdeldrg-smindex)
					drsm=1.0;
				else
					drsm=drdeldrg-smindex;

				accsm=(1.-drsm)*acsmooth_st[smindex]+drsm*acsmooth_st[1+smindex];
				r3inveff=accsm*r3inveff;
				ax[i]=ax[i]-aa*pmass[j]*r3inveff;
				ay[i]=ay[i]-bb*pmass[j]*r3inveff;
				az[i]=az[i]-cc*pmass[j]*r3inveff;
			}
	else
		for (i=1; i<=gd.ndisk; i++) {
              arad= -cmd.bulgmass/rsqr(radsph[i]+cmd.abulge);
              ax[i]=ax[i]+arad*x[i]/(radsph[i]+1.e-10);
              ay[i]=ay[i]+arad*y[i]/(radsph[i]+1.e-10);
              az[i]=az[i]+arad*z[i]/(radsph[i]+1.e-10);
		}

//	free_dvector(acsmooth,0,30001);
	printf("\nforceb:  Bulge accelerations computed <<\n");
}

local void forcebs(void)
{
// Subroutine to compute forces on satellite from the bulge.

//	int ninterp;
	int j,smindex,i;
//	real deldrg,xw,xw2,xw3,xw4;
	real aa,bb,cc,sdrdotdr,rinveff,
		r3inveff,drdeldrg,drsm,accsm,drdotdr,phsm,arad;
//	realptr acsmooth, phsmooth;

//	acsmooth=nr_dvector(0,30001);
//	phsmooth=nr_dvector(0,30001);

//	ninterp=30000;

//	deldrg=2./ninterp;
 
/*	for (i=0; i<=1+ninterp; i++) {
		xw=i*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0) {
			phsmooth[i]=-2.*xw3*(ONE/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.;
			acsmooth[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
		} else {
			phsmooth[i]=-ONE/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.);
			acsmooth[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;
		}
		if (xw >= 2.0) {
			phsmooth[i]=ONE;
			acsmooth[i]=1.0;
		}
	}
*/

	if (cmd.usebulge && cmd.selfgbul && cmd.axibulge) {

		for (j=gd.ndisk+1; j<=gd.ndisk+cmd.nbulge; j++) {
			aa=cmd.xsat-x[j];
			bb=cmd.ysat-y[j];
			cc=cmd.zsat-z[j];
			drdotdr=aa*aa+bb*bb+cc*cc;
			sdrdotdr=rsqrt(drdotdr);
			rinveff=1./(sdrdotdr+1.e-10);
			r3inveff=rinveff/(drdotdr+1.e-10);
			drdeldrg=sdrdotdr*NINTERP/(cmd.epsbulge+cmd.epsbulge);
			smindex=drdeldrg;
			if (NINTERP < smindex) smindex=NINTERP;
			if (1.0 < drdeldrg-smindex)
				drsm=1.0;
			else
				drsm=drdeldrg-smindex;

			phsm=(1.-drsm)*phsmooth_st[smindex]+drsm*phsmooth_st[1+smindex];
			accsm=(1.-drsm)*acsmooth_st[smindex]+drsm*acsmooth_st[1+smindex];
			rinveff=phsm*rinveff;
			r3inveff=accsm*r3inveff;
			gd.potsat=gd.potsat-pmass[j]*rinveff;
			gd.axsat=gd.axsat-aa*pmass[j]*r3inveff;
			gd.aysat=gd.aysat-bb*pmass[j]*r3inveff;
			gd.azsat=gd.azsat-cc*pmass[j]*r3inveff;
		}
 
	} else {

		gd.radsat=rsqrt(cmd.xsat*cmd.xsat+cmd.ysat*cmd.ysat+cmd.zsat*cmd.zsat);
 
		gd.potsat= gd.potsat-cmd.bulgmass/(gd.radsat+cmd.abulge);

		arad= -cmd.bulgmass/rsqr(gd.radsat+cmd.abulge);
		gd.axsat=gd.axsat+arad*cmd.xsat/(gd.radsat+1.e-10);
		gd.aysat=gd.aysat+arad*cmd.ysat/(gd.radsat+1.e-10);
		gd.azsat=gd.azsat+arad*cmd.zsat/(gd.radsat+1.e-10);
	}
 
//	free_dvector(acsmooth,0,30001);
//	free_dvector(phsmooth,0,30001);
	printf("\nforcebs  Bulge-satellite accelerations computed <<\n");
}

local void forced(void)
{
	forced_direct_original();
//	forced_direct_new();
//	forced_tree();
}

local void forced_direct_original(void)
{
//  Subroutine to compute forces on disk particles from the disk's self-
//  gravity. Direct method (N^2).

	int ninterp;
	int i,l,j,smindex;
	real deldrg,xw,xw2,xw3,xw4;
	real aa,bb,cc,sdrdotdr,rinveff,
			r3inveff,drdeldrg,drsm,accsm,drdotdr;
	realptr acsmooth;

	acsmooth=nr_dvector(0,30001);

	ninterp=30000;

	deldrg=2./ninterp;

	for (i=0; i<=1+ninterp; i++) {
		xw=((real)i)*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0)
			acsmooth[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
		else
			acsmooth[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;

		if (xw >= 2.0)
			acsmooth[i]=1.0;
	}

	for (l=1; l<=gd.ndisk; l++) {
           ax[l]=0.;
           ay[l]=0.;
           az[l]=0.;
	}
 
	for (i=1; i<=gd.ndisk-1; i++)
		for (j=i+1; j<=gd.ndisk; j++) {
			aa=x[i]-x[j];
			bb=y[i]-y[j];
			cc=z[i]-z[j];
			drdotdr=aa*aa+bb*bb+cc*cc;
			sdrdotdr=rsqrt(drdotdr);
			rinveff=1./(sdrdotdr+1.e-10);
			r3inveff=rinveff/(drdotdr+1.e-10);
			drdeldrg=sdrdotdr*ninterp/(cmd.epsdisk+cmd.epsdisk);
			smindex=drdeldrg;
			if (ninterp < smindex) smindex=ninterp;
			if (1.0 < drdeldrg-smindex)
				drsm=1.0;
			else
				drsm=drdeldrg-smindex;

			accsm=(1.-drsm)*acsmooth[smindex]+drsm*acsmooth[1+smindex];
			r3inveff=accsm*r3inveff;
			ax[j]=ax[j]+aa*pmass[i]*r3inveff;
			ay[j]=ay[j]+bb*pmass[i]*r3inveff;
			az[j]=az[j]+cc*pmass[i]*r3inveff;
			ax[i]=ax[i]-aa*pmass[j]*r3inveff;
			ay[i]=ay[i]-bb*pmass[j]*r3inveff;
			az[i]=az[i]-cc*pmass[j]*r3inveff;
		}
 
	free_dvector(acsmooth,0,30001);
	printf("\nforced:  Disk accelerations computed <<\n");
}

local void forced_direct_new(void)
{
//  Subroutine to compute forces on disk particles from the disk's self-
//  gravity. Direct method (N^2).

//	int ninterp;
	int i,l,j,smindex;
//	real deldrg,xw,xw2,xw3,xw4;
	real aa,bb,cc,sdrdotdr,rinveff,
			r3inveff,drdeldrg,drsm,accsm,drdotdr;
/*	realptr acsmooth;

	acsmooth=nr_dvector(0,30001);

	ninterp=30000;

	deldrg=2./ninterp;

	for (i=0; i<=1+ninterp; i++) {
		xw=((real)i)*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0)
			acsmooth[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
		else
			acsmooth[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;

		if (xw >= 2.0)
			acsmooth[i]=1.0;
	}
*/

	for (l=1; l<=gd.ndisk; l++) {
           ax[l]=0.;
           ay[l]=0.;
           az[l]=0.;
	}
 
	for (i=1; i<=gd.ndisk-1; i++)
		for (j=i+1; j<=gd.ndisk; j++) {
			aa=x[i]-x[j];
			bb=y[i]-y[j];
			cc=z[i]-z[j];
			drdotdr=aa*aa+bb*bb+cc*cc;
			sdrdotdr=rsqrt(drdotdr);
			rinveff=1./(sdrdotdr+1.e-10);
			r3inveff=rinveff/(drdotdr+1.e-10);
			drdeldrg=sdrdotdr*NINTERP/(cmd.epsdisk+cmd.epsdisk);
			smindex=drdeldrg;
			if (NINTERP < smindex) smindex=NINTERP;
			if (1.0 < drdeldrg-smindex)
				drsm=1.0;
			else
				drsm=drdeldrg-smindex;

			accsm=(1.-drsm)*acsmooth_st[smindex]+drsm*acsmooth_st[1+smindex];
			r3inveff=accsm*r3inveff;
			ax[j]=ax[j]+aa*pmass[i]*r3inveff;
			ay[j]=ay[j]+bb*pmass[i]*r3inveff;
			az[j]=az[j]+cc*pmass[i]*r3inveff;
			ax[i]=ax[i]-aa*pmass[j]*r3inveff;
			ay[i]=ay[i]-bb*pmass[j]*r3inveff;
			az[i]=az[i]-cc*pmass[j]*r3inveff;
		}
 
//	free_dvector(acsmooth,0,30001);
	printf("\nforced:  Disk accelerations computed <<\n");
}

local void forced_tree(void)
{
//  Subroutine to compute forces on disk particles from the disk's self-
//  gravity. Tree scheme (O(N Log N)).

	int i,l,j,smindex;
	real aa,bb,cc,sdrdotdr,rinveff,
			r3inveff,drdeldrg,drsm,accsm,drdotdr;
    bodyptr p, bodytab;
	int nbody;

	for (l=1; l<=gd.ndisk; l++) {
           ax[l]=0.;
           ay[l]=0.;
           az[l]=0.;
	}

// Init tree ...
	gdtreegrav.rsize=1.01;
	gdtreegrav.usequad = TRUE;
	gdtreegrav.eps = 1.0E-04;
	gdtreegrav.theta = 0.7;
	
	nbody = gd.ndisk;
	bodytab = (bodyptr) allocate(nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+nbody) {
		i = p - bodytab + 1;
		Mass(p) = pmass[i];
		Pos(p)[0] = x[i];
		Pos(p)[1] = y[i];
		Pos(p)[2] = z[i];
//		fprintf(stdout,"\n%d %g %g %g",i,x[i],y[i],z[i]);
		Update(p)=TRUE;
	}
	maketree_grav(bodytab, nbody, &gdtreegrav);
	fprintf(stdout,"\n\nncell, tdepth : %d %d\n",gdtreegrav.ncell,gdtreegrav.tdepth);

    gdtreegrav.nbbcalc = gdtreegrav.nbccalc = 0;   
	gdtreegrav.cpuforce=0.;

	normal_gravcalc(bodytab, nbody, &gdtreegrav);
//
 
	for (i=1; i<=gd.ndisk; i++) {
		p = bodytab + i - 1;
		ax[i] = Acc(p)[0];
		ay[i] = Acc(p)[1];
		az[i] = Acc(p)[2];
	}
 
	printf("\nforced (tree):  Disk accelerations computed <<\n");
}


local void forceds(void)
{
// Subroutine to compute forces on satellite center from disk particles.

//	int ninterp;
	int j,smindex,i;
//	real deldrg,xw,xw2,xw3,xw4;
	real aa,bb,cc,sdrdotdr,rinveff,
			r3inveff,drdeldrg,drsm,accsm,drdotdr,phsm;

//	realptr acsmooth, phsmooth;

//	acsmooth=nr_dvector(0,30001);
//	phsmooth=nr_dvector(0,30001);

//	ninterp=30000;

//	deldrg=2./ninterp;

/*	for (i=0; i<=1+ninterp; i++) {
		xw=i*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0) {
			phsmooth[i]=-2.*xw3*(ONE/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.;
			acsmooth[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
		} else {
			phsmooth[i]=-ONE/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.);
			acsmooth[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;
		}
		if (xw >= 2.0) {
			phsmooth[i]=ONE;
			acsmooth[i]=1.0;
		}
	}
*/
 
	gd.potsat=0.0;
	gd.axsat=0.0;
	gd.aysat=0.0;
	gd.azsat=0.0;
 
	for (j=1; j<=gd.ndisk; j++) {
		aa=cmd.xsat-x[j];
		bb=cmd.ysat-y[j];
		cc=cmd.zsat-z[j];
		drdotdr=aa*aa+bb*bb+cc*cc;
		sdrdotdr=rsqrt(drdotdr);
		rinveff=1./(sdrdotdr+1.e-10);
		r3inveff=rinveff/(drdotdr+1.e-10);
		drdeldrg=sdrdotdr*NINTERP/(cmd.epsdisk+cmd.epsdisk);
		smindex=drdeldrg;
		if (NINTERP < smindex) smindex=NINTERP;
		if (1.0 < drdeldrg-smindex)
			drsm=1.0;
		else
			drsm=drdeldrg-smindex;

		phsm=(1.-drsm)*phsmooth_st[smindex]+drsm*phsmooth_st[1+smindex];
		accsm=(1.-drsm)*acsmooth_st[smindex]+drsm*acsmooth_st[1+smindex];
		rinveff=phsm*rinveff;
		r3inveff=accsm*r3inveff;
		gd.potsat=gd.potsat-pmass[j]*rinveff;
		gd.axsat=gd.axsat-aa*pmass[j]*r3inveff;
		gd.aysat=gd.aysat-bb*pmass[j]*r3inveff;
		gd.azsat=gd.azsat-cc*pmass[j]*r3inveff;
	}

//	free_dvector(acsmooth,0,30001);
//	free_dvector(phsmooth,0,30001);
	printf("\nforceds:  Disk-satellite accelerations computed <<\n");
}

local void forceh(void)
{
// Subroutine to compute forces on disk particles from the halo.

	int i,ihalo;
	real xmhr,drhalo,arad;

	if (scanopt(cmd.halotype, "IS")) {

		drhalo=rhalo[5]-rhalo[4];
 
		for (i=1; i<=gd.ndisk; i++) {
			ihalo=radsph[i]/drhalo;
			ihalo=ihalo+1;

			if (radsph[i] > rhalo[gd.ntabhalo])
				xmhr=xmhalo[gd.ntabhalo];
			else
				if (radsph[i] <= rhalo[2])
					xmhr=xmhalo[2]*rqbe(radsph[i])/rqbe(rhalo[2]);
				else
                    xmhr=(radsph[i]-rhalo[ihalo])*xmhalo[ihalo+1]/
						drhalo-(radsph[i]-rhalo[ihalo+1])*
						xmhalo[ihalo]/drhalo;

			arad=-xmhr/(rsqr(radsph[i])+1.e-10);

			ax[i]=ax[i]+arad*x[i]/(radsph[i]+1.e-10);
			ay[i]=ay[i]+arad*y[i]/(radsph[i]+1.e-10);
			az[i]=az[i]+arad*z[i]/(radsph[i]+1.e-10);
		}
 
	} else 
		for (i=1; i<=gd.ndisk; i++) {
			arad= -cmd.halomass/rsqr(radsph[i]+cmd.ahalo);
			ax[i]=ax[i]+arad*x[i]/(radsph[i]+1.e-10);
			ay[i]=ay[i]+arad*y[i]/(radsph[i]+1.e-10);
			az[i]=az[i]+arad*z[i]/(radsph[i]+1.e-10);
		}

	printf("\nforceh:  Halo accelerations computed <<\n");
}

local void forcehs(void)
{
//  Subroutine to compute forces on satellite from the halo.

	int ihalo;
	real xmhr, drhalo, arad;

	gd.radsat=rsqrt(rsqr(cmd.xsat)+rsqr(cmd.ysat)+rsqr(cmd.zsat));
 
	if (scanopt(cmd.halotype, "IS")) {

		drhalo=rhalo[5]-rhalo[4];

		ihalo=gd.radsat/drhalo;
		ihalo=ihalo+1;

		if (gd.radsat > rhalo[gd.ntabhalo]) {
			xmhr=xmhalo[gd.ntabhalo];
			gd.potsat=gd.potsat+uhalo[gd.ntabhalo]*rhalo[gd.ntabhalo]/gd.radsat;
		} else {
			if (gd.radsat <= rhalo[2]) {
				xmhr=xmhalo[2]*rqbe(gd.radsat)/rqbe(rhalo[2]);
				gd.potsat=gd.potsat+uhalo[2];
			} else {
				xmhr=(gd.radsat-rhalo[ihalo])*xmhalo[ihalo+1]/
					drhalo-(gd.radsat-rhalo[ihalo+1])*
					xmhalo[ihalo]/drhalo;
				gd.potsat=gd.potsat+((gd.radsat-rhalo[ihalo])*uhalo[ihalo+1]/
					drhalo-(gd.radsat-rhalo[ihalo+1])*
					uhalo[ihalo]/drhalo);
			}
		}
 
		arad=-xmhr/(rsqr(gd.radsat)+1.e-10);

		gd.axsat=gd.axsat+arad*cmd.xsat/(gd.radsat+1.e-10);
		gd.aysat=gd.aysat+arad*cmd.ysat/(gd.radsat+1.e-10);
		gd.azsat=gd.azsat+arad*cmd.zsat/(gd.radsat+1.e-10);

	} else {

		gd.potsat= gd.potsat-cmd.halomass/(gd.radsat+cmd.ahalo);

		arad= -cmd.halomass/rsqr(gd.radsat+cmd.ahalo);
		gd.axsat=gd.axsat+arad*cmd.xsat/(gd.radsat+1.e-10);
		gd.aysat=gd.aysat+arad*cmd.ysat/(gd.radsat+1.e-10);
		gd.azsat=gd.azsat+arad*cmd.zsat/(gd.radsat+1.e-10);
	}
 
	printf("\nforcehs:  Halo-satellite accelerations computed <<\n");
}


local void setmesh(int nphi, realptr xmesh, realptr ymesh, realptr rmesh)
{
// Subroutine to set up a mesh of x and y points to evaluate the     
// acceleration to compute the acceleration gradient.

	int i,j,n;
    real delphi,delr,r,phi;

	delphi=2.0*PI/((real)(nphi));
	delr=MAX(cmd.rmax,cmd.rmaxgas)/((real)(maxtab));

	r=delr;
	n=0;

	for (i=1; i<=maxtab; i++) {

		phi=0.;
 
		for (j=1; j<=nphi; j++) {
              n=n+1;
              xmesh[n]=r*rcos(phi);
              ymesh[n]=r*rsin(phi);
              rmesh[n]=r;
              phi=phi+delphi;
		}
		r=r+delr;
	}
}

local void setrot(void)
{
// Subroutine to set up centripital, anticlockwise rotation according 
// to the mean rotation curve.

	int i;
	real smallnum;

	smallnum=1.e-07;

	for (i=1; i<=gd.ndisk; i++)
		if (radcyl[i] > smallnum) {
              vx[i]=vx[i]-rotmean[i]*y[i]/radcyl[i];
              vy[i]=vy[i]+rotmean[i]*x[i]/radcyl[i];
		} else
			fprintf(gd.outlog,"\nsetrot:  particle: %d had rad = %g\n",i,radcyl[i]);

	printf("\nsetrot:  Disk rotation established <<\n");
}


local void setsigma(void)
{
// Subroutine to set up the velocity dispersion, correcting for the 
// radial dispersion near the center.

	int i;
	real smallnum,vmean,gaussz,rgauss,gaussr,gaussphi;

	smallnum=1.e-07;

	vmean=0.;

	for (i=1; i<=gd.ndisk; i++) {
		vx[i]=0.0;
		vy[i]=0.0;
		vz[i]=0.0;
		if (sigz[i] != 0.0) {
			gaussz=grandom(vmean,sigz[i]);
			vz[i]=gaussz;
		} else
			if (i > cmd.ndgas)
				printf("\nsetv  Particle : %d had sigz = %g\n",i,sigz[i]);

		if (radcyl[i] > smallnum)
			if (sigr[i] != 0.0) {
                 gaussr=grandom(vmean,sigr[i]);
                 gaussphi=grandom(vmean,sigphi[i]);
                 vx[i]=gaussr*(x[i]/radcyl[i])-gaussphi*(y[i]/radcyl[i]);
                 vy[i]=gaussr*(y[i]/radcyl[i])+gaussphi*(x[i]/radcyl[i]);
			} else
				if (cmd.qsolar != 0.)
					if (i > cmd.ndgas)
                       printf("\nsetv  particle: %d had sigr = %g\n",i,sigr[i]);

		else
			printf("\nsetv   particle : %d had rad = %g\n",i,radcyl[i]);
 
	}
 
	printf("\nsetsigma:  Velocity dispersion set up <<\n");
}

local void sigalar(void)
{
/*
Subroutine to compute the critical Toomre radial velocity     
dispersion at radius R for a given Q from

                             pi Q g SIGMA(R)
           sigma[R]_crit ~ -----------------
                               KAPPA(R)
*/

	int i;
	real term;

	term=cmd.qsolar*PI;

	for (i=1; i<=gd.ndisk; i++)
		if (kappa[i] !=0.)
			sigt[i]=term*surfd[i]/kappa[i];
		else {
			fprintf(gd.outlog,"\nsigalar:  particle: %d had kappa= %g\n",i,kappa[i]);
			sigt[i]=0.;
		}
 
	fprintf(stdout,"\nsigalar:  Toomre velocity dispersion computed <<\n");
}

local void sigcheck(void)
{
/* 
Subroutine to find the radius at which the mean motion becomes 
imaginary. This radius "a" is used to correct the radial velocity 
dispersion profile.  The radius where THETA_mean^2 goes negative 
(R=a) is taken as the scale length of the correcting function:

        <pi^2 ^(1/2)_new = <pi^2 (r=0)_old exp(-sqrt(r^2+2a^2)/2h).

*/

	int nrings,i,nfound,j,isis;
	real delr,r,fnsum,omega,fnmean,sign,asign;

	delr=0.25;
 
	gd.acorrgas=0.0;

	r=cmd.rmax;
	nrings=cmd.rmax/delr;

	for (i=1; i<=nrings-1; i++) {
		fnsum=0.;
		nfound=0;
 
		for (j=cmd.ndgas+1; j<=gd.ndisk; j++)
			if (radcyl[j] <= r && radcyl[j] > r-delr) {
				nfound=nfound+1;
				omega=rotcirc[j]/radcyl[j];
				fnsum=fnsum+1.-2.*(radcyl[j]/gd.h)+(rotcirc[j]*
					rotcirc[j]/(sigr[j]*sigr[j]))-0.25*
					kappa[j]*kappa[j]/(omega*omega);
			}
 
		fnmean=0.;
 
		if (nfound > 0) fnmean=fnsum/((real)(nfound));
 
		if (fnmean >= 0.0 || nfound ==0)
			sign=1.;
		else
			sign=-1.;
 
		if (i == 1)
			asign=sign;
		else
			if (asign != sign) {
				gd.acorr=r-delr/2.;
				printf("\nsigcheck:  Critical radius acorr : %g\n",gd.acorr);
				goto CIENTOUNO;
			} else
				asign=sign;

		r=r-delr;
 
	}
 
	printf("\nsigcheck:  No critical radius found, acorr=0 <<\n");
	gd.acorr=r-delr/2.;

CIENTOUNO: 
 
	for (isis=cmd.ndgas+1; isis<=gd.ndisk; isis++) {
           sigr[isis]=gd.sigr0*rexp(-(rsqrt(radcyl[isis]*radcyl[isis]+
                     2.0*gd.acorr*gd.acorr)/(2.0*gd.h)));
           sigz[isis]=rsqrt(PI*gd.G*cmd.z0*gd.surfd0)*rexp(-(rsqrt(radcyl[isis]*		// g es G?
                     radcyl[isis]+2.0*gd.acorr*gd.acorr)/(2.0*gd.h)));
	}
}


local void meanrot(void)
{
/*
Subroutine to compute the mean circular velocity for each particle
according to the C function and Q under the van der Kruit and 
Searle radial dispersion profile, according to the collisionless 
Boltzmann moment equation.

              2       2      2                           2
         Vmean - Vcirc = < pi  ( 1. - (R/h) - R dln(< pi  ) - C(R/h))
                                                 -----------
                                                     dR
*/

	int i,nerror,nout,ierror,imax;
	real errormax,errorsum,rsuma,rsumb,radmax,xrad,delC,error,
		term1,term2,sigratio,errmean,rmeana,sigmax,t1max,
		t2max,square,rmeanb;

	for (i=1; i<=cmd.ndgas; i++)
		rotmean[i]=rsqrt(rsqr(rotcirc[i]));

	nerror=0;
	errormax=0.;
	errorsum=0.;
	nout=0;
	rsuma=0.;
	rsumb=0.;
	radmax=0.;
 
	for (i=cmd.ndgas+1; i<=gd.ndisk; i++) {
		xrad=radcyl[i]/gd.h;
		sigratio=cfunc(rotcirc[i],radcyl[i],kappa[i]);
 
		if (sigratio > 1.) {
			rsuma=rsuma+radcyl[i];
			nerror=nerror+1;
			delC=sigratio-1.;
			error=delC*100./sigratio;
			errorsum=errorsum+error;

			if (error >= errormax) {
				errormax=error;
				ierror=i;
			}
 
			sigratio=1.;
		}

		term1=rotcirc[i]*rotcirc[i]/(sigr[i]*sigr[i])+1-xrad-
			radcyl[i]*radcyl[i]/(gd.h*rsqrt(radcyl[i]*radcyl[i]+
			2.*gd.acorr*gd.acorr));
		term2=1-xrad-radcyl[i]*radcyl[i]/(gd.h*rsqrt(radcyl[i]*
			radcyl[i]+2.*gd.acorr*gd.acorr));

		if (term1 <= sigratio || term2 >= sigratio) {
			nout=nout+1;
			rsumb=rsumb+radcyl[i];

			if (radcyl[i] >= radmax) {
				radmax=radcyl[i];
				imax=i;
				sigmax=sigratio;
				t1max=term1;
				t2max=term2;
			}

			square=0.; 
		} else
			square=rotcirc[i]*rotcirc[i]-sigr[i]*sigr[i]*
				(sigratio+xrad-1.+radcyl[i]*radcyl[i]/(gd.h*
				rsqrt(radcyl[i]*radcyl[i]+2.*gd.acorr*gd.acorr)));

		rotmean[i]=rsqrt(square);
	}
 
	if (nerror != 0) errmean=errorsum/((real)(nerror));
	if (nerror != 0) rmeana=rsuma/((real)(nerror));

	printf("\nmeanrot:  "); 
	printf("\nmeanrot: %d particles had C > 1",nerror);
	printf("\nmeanrot: mean error  : %g",errmean); 
	printf("\nmeanrot: mean radius : %g",rmeana);
	printf("\nmeanrot:  ");
	printf("\nmeanrot: %d particles out of limits on C",nout);

	if (nout != 0) rmeanb=rsumb/((real)(nout));

	printf("\nmeanrot:  mean radius = %g\n",rmeanb);
}

local void meshacc(int nphi, realptr xmesh, realptr ymesh, realptr rmesh, realptr tabbuff)
{
/*
Subroutine to compute the acceleration at each mesh point by a 
direct sum over the disk particles and by evaluating the halo
contribution at each point.
*/

	int n,i,j,k,ihalo;
	real axm,aym,drhalo,sum,delx,dely,
		delz,term,term2,aradd,xmhr,aradh,aradb;
 
	n=0;
	drhalo=rhalo[5]-rhalo[4];

	for (i=1; i<=maxtab; i++) { 
		sum=0.;

		for (j=1; j<=nphi; j++) {
			n=n+1;
			axm=0.;
			aym=0.;

			for (k=1; k<=gd.ndisk; k++) {
				delx=xmesh[n]-x[k];
				dely=ymesh[n]-y[k];
				delz=z[k];
				term=delx*delx + dely*dely + delz*delz + gd.epsdisk2;
				term2=term*rsqrt(term);
				axm=axm-pmass[k]*delx/term2;
				aym=aym-pmass[k]*dely/term2;
			}

			aradd=(axm*xmesh[n]+aym*ymesh[n])/rmesh[n];

			if (cmd.usebulge && cmd.axibulge && cmd.selfgbul) {

				axm=0.0;
				aym=0.0;

				for (k=gd.ndisk+1; k<=gd.ndisk+cmd.nbulge; k++) {
                    delx=xmesh[n]-x[k];
                    dely=ymesh[n]-y[k];
                    delz=z[k];
                    term=delx*delx + dely*dely + delz*delz 
							+ cmd.epsbulge*cmd.epsbulge;
                    term2=term*rsqrt(term);
                    axm=axm-pmass[k]*delx/term2;
                    aym=aym-pmass[k]*dely/term2;
				}

				aradb=(axm*xmesh[n]+aym*ymesh[n])/rmesh[n];

				aradd=aradd+aradb;

			}

			aradh=0.0;

			if (cmd.usehalo)
				if (scanopt(cmd.halotype, "IS")) {
                    ihalo=rmesh[n]/drhalo;
                    ihalo=ihalo+1;

					if (rmesh[n] > rhalo[gd.ntabhalo])
						xmhr=xmhalo[gd.ntabhalo];
                    else
						if (rmesh[n] <= rhalo[2])
							xmhr=xmhalo[2]*rqbe(rmesh[n])/rqbe(rhalo[2]);
						else
							xmhr=(rmesh[n]-rhalo[ihalo])*xmhalo[ihalo+1]/
								drhalo-(rmesh[n]-rhalo[ihalo+1])*
								xmhalo[ihalo]/drhalo;

                    aradh=-xmhr/(rsqr(rmesh[n])+1.e-10);

				 } else
                    aradh= -cmd.halomass/rsqr(rmesh[n]+cmd.ahalo);
  
			sum=sum+aradd+aradh;
 
			if (cmd.usebulge && (!cmd.axibulge))
				sum=sum-cmd.bulgmass/rsqr(rmesh[n]+cmd.abulge);

		}

		tabbuff[i]=sum/((real)(nphi));

	}
}

local void dadr(void)
{
/*
Subroutine to compute the azimuthally averaged gradient of the 
radial acceleration and store it as a table.  The table is 
contructed by computing the accelerations on a uniform polar mesh 
(maxtab bins in radius and nphi bins in angle). The acceleration at
each mesh point is computed using the n**2 method with softening 
length epsdisk. The accelerations are the azimuthally averaged
and differenced to compute a mean acceleration gradient across 
each radial bin.
*/

	int maxbuff=25000;
	int maxmesh=25000;
	int nphi=100;
	int k;
	realptr xmesh, ymesh, rmesh, tabbuff;
	real delr;

	xmesh=nr_dvector(1,maxmesh);
	ymesh=nr_dvector(1,maxmesh);
	rmesh=nr_dvector(1,maxmesh);
	tabbuff=nr_dvector(1,maxbuff);			//tabbuff contains the azimuthally 
											// averaged radial accelerations

	if (maxbuff < maxtab*nphi || maxmesh < maxtab*nphi)
		error("\nDimension error in dadr\n");
 
	setmesh(nphi,xmesh,ymesh,rmesh);
	meshacc(nphi,xmesh,ymesh,rmesh,tabbuff);

	delr= MAX(cmd.rmax,cmd.rmaxgas)/((real)(maxtab));

	for (k=1; k<=maxtab; k++)
		if (k==1)
			dadrtab[k]=tabbuff[1]/delr;
		else
			dadrtab[k]=(tabbuff[k]-tabbuff[k-1])/delr; 

	free_dvector(xmesh,1,maxmesh);
	free_dvector(ymesh,1,maxmesh);
	free_dvector(rmesh,1,maxmesh);
	free_dvector(tabbuff,1,maxbuff);

	printf("\ndadr:  Radial acceleration gradient computed <<\n");
}

void DiskStat(void)
{
/*
   Subroutine to produce a table of disk characteristics. The 
   analytic values for the rotational velocity and epicyclic frequency
   (softening free) are computed for the particular disk and halo 
   characteristic and stored in vexact and exactk.       
*/

	int maxbin=1000;
	int kk,ilower,iupper,i,j,n;

	int *nbin;

	realptr rbin,surfbin,vcbin,
		sigrbin,sigzbin,vexact,
		exactk,sigpbin,binmean,
		sigtbin,qtoomre,c1bin,
		c2bin,c3bin,kapbin;

	real kapsum,drhalo,atemp,delr,r,va2,surfsum,vcsum,sigrsum,sigzsum,
		sigpsum,sigtsum,summean,csum1,csum2,csum3,omega,rs,bi0,
		bi1,bk0,bk1,vcirc2d,vcirc2h,halokapp,diskkapp,vcirc2b,
		bulgkapp;

    stream outstat_gas, outstat_stars, outtemp;
	char buf[200];

	printf("\nDisk_Stat...\n");

	sprintf(buf,"%s",cmd.filename);
	if(!(outtemp=fopen(buf,"w")))
    	{
      	  error("can't open file `%s`\n",buf);
	} 

	sprintf(buf,"%s-gas",cmd.statfile);
	if(!(outstat_gas=fopen(buf,"w")))
    	{
      	  error("can't open file `%s`\n",buf);
	}
	sprintf(buf,"%s-stars",cmd.statfile);
	if(!(outstat_stars=fopen(buf,"w")))
    	{
      	  error("can't open file `%s`\n",buf);
	} 

	rbin=nr_dvector(1,maxbin);
	surfbin=nr_dvector(1,maxbin);
	vcbin=nr_dvector(1,maxbin);
	sigrbin=nr_dvector(1,maxbin);
	sigzbin=nr_dvector(1,maxbin);
	vexact=nr_dvector(1,maxbin);
	exactk=nr_dvector(1,maxbin);
	sigpbin=nr_dvector(1,maxbin);
	binmean=nr_dvector(1,maxbin);
	sigtbin=nr_dvector(1,maxbin);
	qtoomre=nr_dvector(1,maxbin);
	c1bin=nr_dvector(1,maxbin);
	c2bin=nr_dvector(1,maxbin);
	c3bin=nr_dvector(1,maxbin);
	kapbin=nr_dvector(1,maxbin);

	nbin=nr_ivector(1,maxbin);

	drhalo=rhalo[5]-rhalo[4];

	for (kk=1; kk<=2; kk++) {
 
		if (kk == 1) {
			ilower=1;
			iupper=cmd.ndgas;
			atemp=gd.acorrgas;
		} else {
			ilower=cmd.ndgas+1;
			iupper=gd.ndisk;
			atemp=gd.acorr;
		}
 
		delr=0.5;
		r=0.;
 
// va2 is left from before; this needs to be fixed.
// ------------------------------------------------
 
		va2=1.;
//          va2=vhalo*vhalo
 
		for (i=1; i<=40; i++) {
			n=0;
			kapsum=0.;
			surfsum=0.;
			vcsum=0.;
			sigrsum=0.;
			sigzsum=0.;
			sigpsum=0.;
			sigtsum=0.;
			summean=0.;
			csum1=0.;
			csum2=0.;
			csum3=0.;
 
			for (j=ilower; j<=iupper; j++) {
				if (radcyl[j] > r && radcyl[j] <= r+delr) {
					n=n+1;
                    kapsum=kapsum+kappa[j];
                    surfsum=surfsum+surfd[j];
                    vcsum=vcsum+rotcirc[j];
                    summean=summean+rotmean[j];
                    sigrsum=sigrsum+sigr[j];
                    sigzsum=sigzsum+sigz[j];
                    sigpsum=sigpsum+sigphi[j];
                    sigtsum=sigtsum+sigt[j];
                    omega=rotcirc[j]/radcyl[j];
                    if (omega != 0 && kappa[j] != 0.) {
						csum1=csum1+0.25*kappa[j]*kappa[j]/(omega*omega);
						if (sigr[j] != 0.0)
							csum2=csum2+
								rotcirc[j]*rotcirc[j]/(sigr[j]*sigr[j])+
								1.-(radcyl[j]/gd.h)-(radcyl[j]*radcyl[j]/(
								gd.h*sqrt(radcyl[j]*radcyl[j]+2.*atemp*
								atemp)));

						csum3=csum3+
							1.-(radcyl[j]/gd.h)-(radcyl[j]*radcyl[j]/(
							gd.h*sqrt(radcyl[j]*radcyl[j]+2.*atemp*atemp))); 
					}
				}
			}
 
			nbin[i]=n;
			rbin[i]=r+delr/2.;
			rs=rbin[i]/(2.*gd.h);

			bessel(rs,&bi0,&bi1,&bk0,&bk1);

			vcirc2d=(rbin[i]*rbin[i]*gd.G*gd.diskmass/(2.*rqbe(gd.h)))*
					(bi0*bk0-bi1*bk1);

			vcirc2h=0.0;
			halokapp=0.0;
			vcirc2b=0.0;
			bulgkapp=0.0;

			if (cmd.usehalo)
				if (scanopt(cmd.halotype, "IS")) {
                    vcirc2h=va2*(1.-(cmd.gamhalo/rbin[i])*atan(rbin[i]/
							cmd.gamhalo));
                    halokapp=(va2/(rbin[i]*rbin[i]))*
							(2.-(cmd.gamhalo/rbin[i])*atan(rbin[i]/
							cmd.gamhalo)-1./(1.+(rbin[i]/cmd.gamhalo)*
							(rbin[i]/cmd.gamhalo)));
				} else {
					vcirc2h=(cmd.halomass*rbin[i])/rsqr(rbin[i]+cmd.ahalo);
                    halokapp=(cmd.halomass*(rbin[i]+3.*cmd.ahalo)/(rbin[i]*
						rqbe(rbin[i]+cmd.ahalo)));
				}
 
			if (cmd.usebulge) {
				vcirc2b=(cmd.bulgmass*rbin[i])/rsqr(rbin[i]+cmd.abulge);
				bulgkapp=(cmd.bulgmass*(rbin[i]+3.*cmd.abulge)/(rbin[i]*
							rqbe(rbin[i]+cmd.abulge)));
			}

			vexact[i]=rsqrt(vcirc2d+vcirc2h+vcirc2b);
			diskkapp=(gd.G*gd.diskmass/rqbe(gd.h))*(
						2.*bi0*bk0-bi1*bk1+(rbin[i]/(2.*gd.h))*(
						bi1*bk0-bi0*bk1));
			exactk[i]=rsqrt(halokapp+diskkapp+bulgkapp);

			if (n != 0) {
				kapbin[i]=kapsum/((real)(n));
				surfbin[i]=surfsum/((real)(n));
				vcbin[i]=vcsum/((real)(n));
				sigrbin[i]=sigrsum/((real)(n));
				sigzbin[i]=sigzsum/((real)(n));
				sigpbin[i]=sigpsum/((real)(n));
				sigtbin[i]=sigtsum/((real)(n));
				if (cmd.qsolar != 0.)
					qtoomre[i]=cmd.qsolar*sigrbin[i]/sigtbin[i];
				else
					qtoomre[i]=0.;

				binmean[i]=summean/((real)(n));
				c1bin[i]=csum1/((real)(n));
				c2bin[i]=csum2/((real)(n));
				c3bin[i]=csum3/((real)(n));
			} else {
				kapbin[i]=0.;
				surfbin[i]=0.;
				vcbin[i]=0.;
				sigrbin[i]=0.;
				sigzbin[i]=0.;
				sigpbin[i]=0.;
				sigtbin[i]=0.;
				binmean[i]=0.;
				c1bin[i]=0.;
				c2bin[i]=0.;
				c3bin[i]=0.;
			}
 
			r=r+delr;

		} 

		if (kk==1) {
		fprintf(outstat_gas,"%s %s\n#",
				"#R [kpc]  N  SIGMA  KAPPA  Kexact  Vcirc  Vexact ", 
				"Vmean  sigma r  sigma T  Q  sigma phi  sigma z");
 
		for (i=1; i<=40; i++)
			fprintf(outstat_gas,"\n%g %d %g %g %g %g %g %g %g %g %g %g %g",
				rbin[i],nbin[i],surfbin[i],kapbin[i],
				exactk[i],vcbin[i],vexact[i],
				binmean[i],sigrbin[i],sigtbin[i],
				qtoomre[i],sigpbin[i],sigzbin[i]);
		} else {
		fprintf(outstat_stars,"%s %s\n#",
				"#R [kpc]  N  SIGMA  KAPPA  Kexact  Vcirc  Vexact ",
				"Vmean  sigma r  sigma T  Q  sigma phi  sigma z");
 
		for (i=1; i<=40; i++)
			fprintf(outstat_stars,"\n%g %d %g %g %g %g %g %g %g %g %g %g %g",
				rbin[i],nbin[i],surfbin[i],kapbin[i],
				exactk[i],vcbin[i],vexact[i],
				binmean[i],sigrbin[i],sigtbin[i],
				qtoomre[i],sigpbin[i],sigzbin[i]);
		}

		for (i=1; i<=40; i++)
			fprintf(outtemp,"\n%g %g %g %g",rbin[i],c3bin[i],c1bin[i],c2bin[i]);

	}

	fclose(outtemp);
	fclose(outstat_gas);
	fclose(outstat_stars);

	free_dvector(rbin,1,maxbin);
	free_dvector(surfbin,1,maxbin);
	free_dvector(vcbin,1,maxbin);
	free_dvector(sigrbin,1,maxbin);
	free_dvector(sigzbin,1,maxbin);
	free_dvector(vexact,1,maxbin);
	free_dvector(exactk,1,maxbin);
	free_dvector(sigpbin,1,maxbin);
	free_dvector(binmean,1,maxbin);
	free_dvector(sigtbin,1,maxbin);
	free_dvector(qtoomre,1,maxbin);
	free_dvector(c1bin,1,maxbin);
	free_dvector(c2bin,1,maxbin);
	free_dvector(c3bin,1,maxbin);
	free_dvector(kapbin,1,maxbin);

	free_ivector(nbin,1,maxbin);
}

void DiskVel(void)
{
//  Subroutine to initialize velocities of disk particles.
 
	printf("\nDisk_Vel...\n");

	set_interpolation_force_table();
	forced();
	if (cmd.usebulge) forceb();
	if (cmd.usehalo) forceh();
	radacc();
	dadr();
	getkappa();
	sigalar();
	sigmar();
	sigmaz();
	circv();
	sigcheck();
	sigmap();
	setsigma();
	meanrot();
	setrot();
}

local void getkappa(void)
{
//   Subroutine to compute the epicyclic frequency distribution.  The 
//   epicyclic frequency is given by :
//
//                               3     d PHI         d^2  PHI
//           KAPPA(R) = sqrt [ -----  --------  +   ---------- ]
//                               R     d R           d  R^2
//
//   where PHI is the potential at R. 

	int i,iplace;
	real smallnum,brack,dadr;

	smallnum=1.e-07;

	for (i=1; i<=gd.ndisk; i++) {
		iplace=(radcyl[i]/MAX(cmd.rmax,cmd.rmaxgas))*((real)(maxtab))+1.;
		if (iplace == maxtab+1) {
			printf("\ngetkappa  iplace=%d for i=%d <<",iplace,i);
			iplace=maxtab;
		}

		dadr=dadrtab[iplace];

		if (radcyl[i] > smallnum) {
			brack=-3.*(aradcyl[i]/radcyl[i])-dadr;
			if (brack < 0.) {
				fprintf(gd.outlog,"\ngetkappa:  particle : %d kappa^2=%g",i,brack);
				fprintf(gd.outlog,"\ngetkappa:                          rad = %g\n",
							radcyl[i]);
				kappa[i]=0.;
			} else
				kappa[i]=rsqrt(brack);

		} else {
			fprintf(gd.outlog,"\ngetkappa:  particle: %d had rad= %g\n",i,radcyl[i]);
			kappa[i]=0.;
		}
 
	}
       
	printf("\ngetkappa  Kappa computed <<\n");
}


local void halopot(void)
{
// Subroutine to compute potentials of halo particles.

	int ninterp,j,smindex,i,ihalo;
	real deldrg,xw,xw2,xw3,xw4,aa,bb,cc,sdrdotdr,rinveff,drdeldrg,
		drsm,drdotdr,phsm,drhalo;

	realptr phsmooth;

    double cpustart;
    cpustart = cputime();

	phsmooth=nr_dvector(0,30001);

	ninterp=30000;

	deldrg=2./ninterp;
 
	for (i=0; i<=1+ninterp; i++) {
		xw=i*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0)
			phsmooth[i]=-2.*xw3*(ONE/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.;
		else
			phsmooth[i]=-ONE/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.);

		if (xw >= 2.0)
			phsmooth[i]=ONE;
	}
 
	for (i=gd.nbodies-cmd.nhalo+1-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++)
		pot[i]=0.0;

//   Compute halo-disk interaction.
//   ------------------------------
 
	for (i=gd.nbodies-cmd.nhalo+1-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++) {
		for (j=1; j<=gd.ndisk; j++) {
			aa=x[i]-x[j];
			bb=y[i]-y[j];
			cc=z[i]-z[j];
			drdotdr=aa*aa+bb*bb+cc*cc;
			sdrdotdr=rsqrt(drdotdr);
			rinveff=1./(sdrdotdr+1.e-10);
			drdeldrg=sdrdotdr*ninterp/(cmd.epsdisk+cmd.epsdisk);
			smindex=drdeldrg;								// Correcto?
			if (ninterp < smindex) smindex=ninterp;
			if (1.0 < drdeldrg-smindex)
				drsm=1.0;
			else
				drsm=drdeldrg-smindex;

			phsm=(1.-drsm)*phsmooth[smindex]+drsm*phsmooth[1+smindex];
			rinveff=phsm*rinveff;
			pot[i]=pot[i]-pmass[j]*rinveff;
		}

		if (cmd.usebulge && cmd.selfgbul && cmd.axibulge)
			for (j=gd.ndisk+1; j<=gd.ndisk+cmd.nbulge; j++) {
				aa=x[i]-x[j];
				bb=y[i]-y[j];
				cc=z[i]-z[j];
				drdotdr=aa*aa+bb*bb+cc*cc;
				sdrdotdr=rsqrt(drdotdr);
				rinveff=1./(sdrdotdr+1.e-10);
				drdeldrg=sdrdotdr*((real) ninterp)/(cmd.epsbulge+cmd.epsbulge);
				smindex=drdeldrg;
				if (ninterp < smindex) smindex=ninterp;
				if (1.0 < drdeldrg-smindex)
                    drsm=1.0;
				else
					drsm=drdeldrg-smindex;

				phsm=(1.-drsm)*phsmooth[smindex]+drsm*phsmooth[1+smindex];
				rinveff=phsm*rinveff;
				pot[i]=pot[i]-pmass[j]*rinveff;
			}
 
	}
 
//   Compute halo self-interaction.
//   ------------------------------
 
	if (scanopt(cmd.halotype, "IS")) {

		drhalo=rhalo[5]-rhalo[4];
 
		for (i=gd.nbodies-cmd.nhalo+1-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++) {
			ihalo=radsph[i]/drhalo;
			ihalo=ihalo+1;

			if (radsph[i] > rhalo[gd.ntabhalo])
                 pot[i]=pot[i]+uhalo[gd.ntabhalo]*rhalo[gd.ntabhalo]/radsph[i];
			else
				if (radsph[i] <= rhalo[2])
					pot[i]=pot[i]+uhalo[2];
				else
					pot[i]=pot[i]+((radsph[i]-rhalo[ihalo])*
						uhalo[ihalo+1]/drhalo-(radsph[i]-
						rhalo[ihalo+1])*uhalo[ihalo]/drhalo);
		}
 
	} else
		for (i=gd.nbodies-cmd.nhalo+1-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++)
              pot[i]= pot[i]-cmd.halomass/(radsph[i]+cmd.ahalo);

// Compute halo-bulge interaction
// ------------------------------

	if (cmd.usebulge && cmd.selfgbul && (!cmd.axibulge))
		for (i=gd.nbodies-cmd.nhalo+1-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++)
			pot[i]= pot[i]-cmd.bulgmass/(radsph[i]+cmd.abulge);

	printf("\nhalopot  Halo potentials computed <<\n");
 
	free_dvector(phsmooth,0,30001);

	fprintf(stdout,"\n\nhalopot CPU time: %g\n",cputime()-cpustart);
	fflush(stdout);
}

void HaloVel(void)
{
/*
   Subroutine to initialize velocities of halo particles, using the
   spherical Jeans equations.  That is, the radial velocity dispersion
   at radius r is:

                            infinity
                               /
                         1     |  G M(r)
	   <v_r ^2 (r) =  -------  |  ------ rho_h(r) dr
                      rho_h(r) /   r^2
                               r


   where rho_h(r) is the halo density at r.
*/
// Se llama a: cputime(), error(), halopot()

	int maxnbin=10000;

	int nbin,i,j,ir,ntest,irlower,irupper,irindex;
	int irtmp;

	realptr xintbin, xmbin, gridmass;

	real radmax,alphhalo,qhalo,ria,ri,rj,r,dr,rhoi,rhoj,rja,p,
		xv,gspeed,cth,sth,signs,phi,tmass,fracmass,epsilon,rlower,
		rupper,vescape,sigfix,vmaxsig,rkludge,vmax,potrmax,vkludge;

    double cpustart;
    cpustart = cputime();

	xintbin=nr_dvector(0,maxnbin+1000);
	xmbin=nr_dvector(0,maxnbin+1000);
	gridmass=nr_dvector(0,maxnbin+1000);

	qhalo=cmd.gamhalo/cmd.rthalo;
	alphhalo=1./(1.-rsqrt(PI)*qhalo*rexp(qhalo*qhalo)*erfcc(qhalo));

	nbin=1000;

	if (nbin >= maxnbin) error("\nnbin error in halovel\n");
 
	radmax=1.01*cmd.rmaxhalo;			//	LINEA ORIGINAL
//	radmax=1.02*rmaxhalo;

	dr=radmax/nbin;

	for (i=0; i<=nbin; i++) {
		xintbin[i]=0.0;
		xmbin[i]=0.0;
		gridmass[i]=0.0;
	}
 
	rkludge=0.75;

	for (i=1; i<=gd.nbodies-cmd.nsat; i++) {				// SATELLITE NOT INCLUDED...
		r=rsqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
		epsilon=cmd.epsdisk;
		if (i > gd.ndisk && i < gd.nbodies-cmd.nhalo+1-cmd.nsat) epsilon=cmd.epsbulge;
		if (i > gd.nbodies-cmd.nhalo-cmd.nsat) epsilon=cmd.epshalo;

		if (i <= gd.ndisk) r=r*rkludge;

		rlower=r-epsilon/2.0;
		rupper=r+epsilon/2.0;
		ir=r/dr;
		ir=ir+1;
		irlower=rlower/dr;
		irlower=irlower+1;

		if (rlower < 0.0) {
			irlower=rabs(rlower)/dr;
			irlower= -irlower;
		}
 
		irupper=rupper/dr;
		irupper=irupper+1;
		tmass=0.0;

		for (j=irlower; j<=irupper; j++) {
			irindex=j;
			if (j < ir) irindex=ir;
			if (j > maxnbin+1000) irindex=maxnbin+1000;

			if (irlower == irupper) {
                 tmass=tmass+pmass[i];
                 gridmass[irindex]=gridmass[irindex]+pmass[i];
			} else {
				if (j != irlower && j != irupper) {
                    fracmass=(dr/epsilon)*pmass[i];
                    tmass=tmass+fracmass;
                    gridmass[irindex]=gridmass[irindex]+fracmass;
				} else {
					if (j == irupper) {
                       fracmass=((rupper-(irupper-1)*dr)/epsilon)*pmass[i];
                       tmass=tmass+fracmass;
                       gridmass[irindex]=gridmass[irindex]+fracmass;
					} else {
                       fracmass=((irlower*dr-rlower)/epsilon)*pmass[i];
                       tmass=tmass+fracmass;
                       gridmass[irindex]=gridmass[irindex]+fracmass;
					}
				}
			}
 
		}

		if (rabs(tmass-pmass[i])/pmass[i] > 1.e-4)
			error("\nmass assignment error in halovel\n");
 
	}
 
	for (i=1; i<=nbin; i++)
		xmbin[i]=xmbin[i-1]+gridmass[i];

	xmbin[nbin+1]=xmbin[nbin];
 
	if (scanopt(cmd.halotype,"LH")) {
 
		for (i=1; i<=nbin; i++) {
 
              ri=i*dr;
              ria=ri/cmd.ahalo;
              rhoi=cmd.halomass/(2.*PI*rqbe(cmd.ahalo)*ria*rqbe(1.+ria));
 
			for (j=i; j<=nbin; j++) {
				rj=j*dr+0.5*dr;
				rja=rj/cmd.ahalo;
				rhoj=cmd.halomass/(2.*PI*rqbe(cmd.ahalo)*rja*rqbe(1.+rja));
				xintbin[i]=xintbin[i]+rhoj*0.5*(xmbin[j]+
						xmbin[j+1])*dr/(rj*rj);
		}

		xintbin[i]=xintbin[i]/rhoi;

		}
   
	} else {
 
		for (i=1; i<=nbin; i++) {
			ri=i*dr;
			rhoi=cmd.halomass*alphhalo*rexp(-ri*ri/rsqr(cmd.rthalo))/
				(2.*PI*rsqrt(PI)*cmd.rthalo*(ri*ri+rsqr(cmd.gamhalo)));
 
			for (j=i; j<=nbin; j++) {
				rj=j*dr+0.5*dr;
				rhoj=cmd.halomass*alphhalo*rexp(-rj*rj/rsqr(cmd.rthalo))/
					(2.*PI*rsqrt(PI)*cmd.rthalo*(rj*rj+rsqr(cmd.gamhalo)));
				xintbin[i]=xintbin[i]+rhoj*0.5*(xmbin[j]+
						xmbin[j+1])*dr/(rj*rj);
			}
 
			xintbin[i]=xintbin[i]/rhoi;
 
		}
 
	}

	xintbin[nbin+1]=0.0;

	for (i=gd.nbodies-cmd.nhalo+1-cmd.nsat; i<=gd.nbodies-cmd.nsat; i++) {
		r=rsqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
		ir=r/dr;
		irtmp=MIN(ir,nbin);
		ir=irtmp;
		sigr[i]=((xintbin[ir+1]-xintbin[ir])*r+xintbin[ir]*
				(ir+1)*dr-xintbin[ir+1]*ir*dr)/dr;
		sigr[i]=rsqrt(sigr[i]);
	}
 
//   Initialize velocities isotropically so that the distribution of
//   speeds is proportional to v^2 EXP[-v^2/(2*sigma_r^2)].  Limit
//   speed to the local escape speed.
//   ---------------------------------------------------------------
 
	halopot();
 
	ntest=gd.nbodies-cmd.nhalo+1-cmd.nsat;
 
	vkludge=0.95;

CIENTOSESENTA:

	vescape=rsqrt(2.0*rabs(pot[ntest]));

	potrmax= -(cmd.halomass+cmd.bulgmass+gd.diskmass)/cmd.rmaxhalo;
	vmax=vescape*vescape+2.0*potrmax;

	if (vmax > 0.0) {
		vmax=rsqrt(vmax);
 
		xv=xrandom(0.0,1.0)*vescape;

		vmaxsig=vkludge*vmax/(rsqrt(2.0)*sigr[ntest]);
		sigfix=1.0-erfcc(vmaxsig)-8.0*vmaxsig*(0.75+0.5*rsqr(vmaxsig))*
				rexp(-rsqr(vmaxsig))/(3.0*rsqrt(PI));
		sigfix=sigr[ntest]/rsqrt(sigfix);

		gspeed=xv*xv*rexp(-xv*xv/(2.0*rsqr(sigfix)));

		gspeed=gspeed/(2.0*sigfix*sigfix*rexp(-1.0));

		p=xrandom(0.0,1.0);
	} else {
           gspeed=0.0;
           p=0.0;
	}
 
	if (p <= gspeed || vmax <= 0.0) {
		cth=2.0*(xrandom(0.0,1.0)-0.5);
		sth=rsqrt(1.0-cth*cth);
		signs=2.0*(xrandom(0.0,1.0)-0.5);
		cth=signs*cth/rabs(signs);
		phi=TWOPI*xrandom(0.0,1.0);
		vx[ntest]=xv*sth*rcos(phi);
		vy[ntest]=xv*sth*rsin(phi);
		vz[ntest]=xv*cth;

		if (rsqrt(rsqr(vx[ntest])+rsqr(vy[ntest])+rsqr(vz[ntest])) >
			vkludge*vmax && vmax > 0.0) goto CIENTOSESENTA;

		ntest=ntest+1;
	} else
		goto CIENTOSESENTA;

	if (ntest <= gd.nbodies-cmd.nsat) goto CIENTOSESENTA;

	free_dvector(xintbin,0,maxnbin+1000);
	free_dvector(xmbin,0,maxnbin+1000);
	free_dvector(gridmass,0,maxnbin+1000);

	fprintf(stdout,"\n\nHaloVel CPU time: %g\n",cputime()-cpustart);
	fflush(stdout);
}


local void set_interpolation_force_table(void)
{
	real deldrg,xw,xw2,xw3,xw4;
	int i;

	deldrg=2./NINTERP;

	for (i=0; i<=1+NINTERP; i++) {
		xw=((real)i)*deldrg;
		xw2=xw*xw;
		xw3=xw2*xw;
		xw4=xw2*xw2;
		if (xw <= 1.0) {
			phsmooth_st[i]=-2.*xw3*(ONE/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.;
			acsmooth_st[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
		} else {
			phsmooth_st[i]=-ONE/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.);
			acsmooth_st[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;
		}

		if (xw >= 2.0) {
			phsmooth_st[i]=ONE;
			acsmooth_st[i]=1.0;
		}
	}
}

