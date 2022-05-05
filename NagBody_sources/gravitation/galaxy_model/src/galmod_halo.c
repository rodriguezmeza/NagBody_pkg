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

local void inhmass(void);
local void sethalo(void);
local void cmhalo(void);


// ---------------------------------------------


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


