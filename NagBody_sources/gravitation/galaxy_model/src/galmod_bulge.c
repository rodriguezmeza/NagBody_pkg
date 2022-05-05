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

local void inbmass(void);
local void setbulge(void);
local void cmbulge(void);

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


