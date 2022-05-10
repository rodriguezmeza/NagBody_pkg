/* =============================================================================
	MODULE: galmod_general.c		[galaxy_models]
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
	Copyright: (c) 1999-2006 Mar.  All Rights Reserved.
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "stdinc.h"
#include "constant.h"
#include "mathfns.h"
#include "numrec.h"

//#include "globaldefs.h"
#include "protodefs.h"

local void oblatep(real, real, realptr);
local void obsigt(real, real, real, realptr);
local void obsigp(real, real, real, realptr);
local void rhopsiR(real, real, realptr, realptr);
local void obdpotdz(real, real, realptr);
local void obdpotdR(real, real, realptr);

local real G, M, a, c;						// Tener cuidado con el uso interno 
											// de estas variables
local int nsimp;

static R, z;								// Se comparten con rhopsiR (abajo)

void obsigma(real xR, real xz, realptr sigR, realptr sigz, realptr sigp, 
					real abulge, real cbulge, real bmass, int ns, 
					real rmaxbulg, real zmaxbulg)
{
/*
Subroutine to compute velocity dispersions for an axisymmetric
mass distribution, according to:

                                          psi
                                           /
                                   1       |
       sigma_R^2 = sigma_z^2 = ----------  | rho(psi,R) d psi ,
                               rho(psi,R)  |
                                           /
                                           0

                                             psi
                                              /
                                      1       |   d(rho(psi,R))
       sigma_phi^2 = sigma_R^2 +  ----------  | R ------------- d psi,
                                  rho(psi,R)  |        dR
                                              /
                                              0


where psi = - potential.
*/
	real pot,psi,sigt,rho,mRz;
	real sigRtmp, sigptmp;

	R=xR;
	z=xz;

	a=abulge;
	c=cbulge;
	M=bmass;

	G=1.0;
 
	nsimp=ns;

	oblatep(xR,xz, &pot);

	psi=-pot;

	obsigt(psi,xR,xz, &sigt);
	obsigp(psi,xR,xz, &sigptmp);
 
	mRz=rsqrt(rsqr(xR)/rsqr(a)+rsqr(xz)/rsqr(c));
	rho=M/(2.*PI*a*a*c*mRz*rpow(1.+mRz,3));

	sigRtmp=sigt/rho;
	sigptmp=sigRtmp+sigptmp/rho;
	*sigR=rsqrt(sigRtmp);
 
	*sigz=*sigR;

	*sigp=rsqrt(sigptmp);
}


#define nscoef		1000

local void oblatep(real R, real z, realptr pot)
{
/*
Subroutine to determine the negative potential psi from the
integral:

                       infinity
                          /
                    G M   |              1
       psi (R,z) = -----  | du  -------------------- ,
                     2    |     delta(u) (1+m(u))**2
                          /
                          0

where:

                    R^2         z^2
       m^2(u) =  --------- + ---------
                  a^2 + u     c^2 + u


       delta (u) = (a^2 + u) SQRT(c^2 + u).
*/
	int i;
	real pot1;
	static bool firstc = TRUE;
	static realptr simpcoef, xarg, xarg2, xarg6;
	static a2, c2;
	real pottmp;

	if (firstc) {
		firstc=FALSE;

		simpcoef=nr_dvector(1,nscoef);
		xarg=nr_dvector(1,nscoef);
		xarg2=nr_dvector(1,nscoef);
		xarg6=nr_dvector(1,nscoef);
	
		if (2*nsimp+1>nscoef) error("\noverflow in oblatep\n");
 
		simpcoef[1]=1.0;

		for (i=2; i<=2*nsimp; i+=2) {
              simpcoef[i]=4.0;
              simpcoef[i+1]=2.0;
		}
 
		simpcoef[2*nsimp+1]=1.0;
 
		for (i=1; i<=2*nsimp+1; i++) {
              xarg[i]=( (real)i-1.0)/(2.0* (real)nsimp);
              xarg2[i]=rsqr(xarg[i]);
              xarg6[i]=rpow(xarg2[i],3);
		}

		a2=a*a;
		c2=c*c;
	}

	pottmp=0.0;
	pot1=0.0;

	for (i=1; i<=2*nsimp+1; i++) {
		pottmp=pottmp+simpcoef[i]/
			(
			rsqr(1.0+rsqrt(R*R/(a2+xarg[i])+z*z/(c2+xarg[i]))) 
				* (a2+xarg[i]) * rsqrt(c2+xarg[i])
			);

		pot1=pot1+simpcoef[i]*xarg2[i]/
			(
			rsqr(1.0
			+rsqrt(
				R*R*xarg6[i]/(xarg6[i]*a2+1)+z*z*xarg6[i]/(xarg6[i]*c2+1)
			)
			)
			* (a2*xarg6[i]+1.0)*rsqrt(c2*xarg6[i]+1.0)
			);
	}

	*pot = -0.5*G*M*(pottmp+6.*pot1)/(2.0*((real)nsimp)*3.0);
 }


static realptr rhosave, zsave;

local void obsigt(real psi, real R, real z, realptr sigt)
{
/*
Subroutine to determine the integral:

           psi
            /
            |
       I =  | d-psi rho(psi,R).
            |
            /
            0
*/
	int i;
	real rho,znew;
	static bool firstc = TRUE;
	static realptr simpcoef, xarg;
	real sigttmp;
 
	if (firstc) {
		firstc=FALSE;
 
		simpcoef=nr_dvector(1,nscoef);
		xarg=nr_dvector(1,nscoef);
		rhosave=nr_dvector(1,nscoef);
		zsave=nr_dvector(1,nscoef);

		if (2*nsimp+1>nscoef) error("\noverflow in oblatep\n");

		simpcoef[1]=1.0;

		for (i=2; i<=2*nsimp; i+=2) {
              simpcoef[i]=4.0;
              simpcoef[i+1]=2.0;
		}
 
		simpcoef[2*nsimp+1]=1.0;

		for (i=1; i<=2*nsimp+1; i++)
              xarg[i]=(((real)i)-1.0)/(2.0*((real)nsimp));
	}
 
	sigttmp=0.0;

	for (i=1; i<=2*nsimp+1; i++) {
		rhopsiR(R,psi*xarg[i],&rho,&znew);
		sigttmp=sigttmp+simpcoef[i]*rho;
		rhosave[i]=rho;
		zsave[i]=znew;
	}

	*sigt=sigttmp*psi/(2.0*((real)nsimp)*3.0);
}

local void obsigp(real psi, real R, real z, realptr sigp)
{
/*
Subroutine to determine the integral:

           psi
            /
            |         partial-rho(psi,R)
       I =  | d-psi R ------------------ .
            |             partial-R
            /
            0
*/
	int i;
	real mRz,drhodR,drhodz,dpotdz,dpotdR;
	static bool firstc = TRUE;
	static realptr simpcoef;
	real sigptmp;
 
	if (firstc) {
		firstc=FALSE;
 
		simpcoef=nr_dvector(1,nscoef);

		if (2*nsimp+1>nscoef) error("\noverflow in oblatep\n");
 
		simpcoef[1]=1.0;

		for (i=2; i<=2*nsimp; i+=2) {
              simpcoef[i]=4.0;
              simpcoef[i+1]=2.0;
		}
 
           simpcoef[2*nsimp+1]=1.0;

	}

	sigptmp=0.0;
 
	for (i=1; 2*nsimp+1; i++) {
		mRz=rsqrt(rsqr(R)/rsqr(a)+rsqr(zsave[i])/rsqr(c));
 
		drhodR=-rhosave[i]*(1.+4.*mRz)*R/(rsqr(mRz)*rsqr(a)*(1.+mRz));
		drhodz=-rhosave[i]*(1.+4.*mRz)*zsave[i]/(rsqr(mRz)*rsqr(c)*(1.+mRz));
 
		obdpotdz(R,zsave[i],&dpotdz);
		obdpotdR(R,zsave[i],&dpotdR);

		drhodR=drhodR-drhodz*dpotdR/dpotdz;

		sigptmp=sigptmp+simpcoef[i]*R*drhodR;
 
	}
 
	*sigp=sigptmp*psi/(2.0*nsimp*3.0);
}

local void rhopsiR(real xR, real xpsi, realptr rho, realptr znew)
{
/*
Subroutine to compute rho (psi,R), given psi = psi(R,z) and 
rho = rho(R,z).  From input values of R and psi, the relation
psi(R,z) is inverted numerically to provide z = z(psi,R).
The density is then computed from R and z.
*/

	int nit;
	real pot,zinit,zold,fz,
		dfdz,mRz,dpotdz,pottest,zlow,zhigh,potlow,pothigh,ztest,potmax;

	if (xpsi==0.0) {
		*rho=0.0;
		*znew=1.e10;
		return;
	}

	oblatep(xR, ZERO, &potmax);

	if (xpsi==-potmax) {
		*znew=ZERO;
		return;
	}

	if (xpsi>-potmax) {
		printf("\nrange error in rhopsiR\n");
		printf("psimax = %g\n",-potmax);
		error(" ");
	}
 
	zinit=rsqr(G*M/xpsi-a)-rsqr(xR);

	if (zinit<=1.e-4)
		zinit=c;

	zinit=rsqrt(zinit);

	zold=zinit;
	*znew=zinit;

	nit=0;

NDIEZ:

	oblatep(xR,zold, &pot);
	fz= -pot - xpsi;
	obdpotdz(xR,zold, &dpotdz);
	dfdz = -dpotdz;
	*znew = zold - fz/dfdz;

	if (rabs(zold-*znew)/zold > 1.e-6) { 
           zold=*znew; 
           nit=nit+1; 
           if (nit > 20) goto QUINCE;
           goto NDIEZ;
	}

QUINCE:

	oblatep(xR,*znew, &pottest);
 
	if (rabs(-pottest-xpsi)/xpsi > 1.e-6) {
 
		zlow=0.0;
		zhigh=1000.;
 
		oblatep(xR,zlow, &potlow);
		oblatep(xR,zhigh, &pothigh);

		if ((xpsi+potlow)*(xpsi+pothigh) >= 0.0) {
			printf("\nlimit error in rhopsiR\n");
			exit(1);
		}
 
		nit=0;

DIEZYSEIS:
 
		ztest=0.5*(zlow+zhigh);
 
		oblatep(xR,ztest, &pottest);
 
		if (rabs(pottest+xpsi)/xpsi <= 1.0e-6) goto NVEINTE;

		if (xpsi <= -pottest)
              zlow=ztest;
		else
              zhigh=ztest;
 
		nit=nit+1;

		if (nit > 100) {
			printf("\nnit1 error in rhopsiR\n");
			exit(1);
		}
 
		goto DIEZYSEIS;

NVEINTE:

		*znew=ztest;
 
	}

	mRz=rsqrt(rsqr(xR)/rsqr(a)+rsqr(*znew)/rsqr(c));
	*rho=M/(2.*PI*a*a*c*mRz*rpow(1.+mRz,3.0));
}


local void obdpotdz(real R, real z, realptr dpotdz)
{
//     Subroutine to compute partial-psi / partial-z.

	int i;
	real dpotdz1,mofu;
	static bool firstc = TRUE;
	static realptr simpcoef, xarg, xarg2, xarg3, xarg5, xarg6;
	static a2, c2;
	real dpotdztmp;

	if (firstc) {
		firstc=FALSE;

		simpcoef=nr_dvector(1,nscoef);
		xarg=nr_dvector(1,nscoef);
		xarg2=nr_dvector(1,nscoef);
		xarg3=nr_dvector(1,nscoef);
		xarg5=nr_dvector(1,nscoef);
		xarg6=nr_dvector(1,nscoef);
 
		if (2*nsimp+1 > nscoef) error("\noverflow in oblatep\n");
 
		simpcoef[1]=1.0;

		for (i=2; i<=2*nsimp; i+=2) {
			simpcoef[i]=4.0;
			simpcoef[i+1]=2.0;
		}
 
		simpcoef[2*nsimp+1]=1.0;

		for (i=1; i<=2*nsimp+1; i++) {
              xarg[i]=(((real) i)-1.0)/(2.0*((real) nsimp));
              xarg2[i]=rsqr(xarg[i]);
              xarg3[i]=xarg2[i]*xarg[i];
              xarg5[i]=xarg2[i]*xarg3[i];
              xarg6[i]=rpow(xarg2[i],3);
		}
 
		a2=a*a;
		c2=c*c;
	}

	dpotdztmp=0.0;
	dpotdz1=0.0;
 
	for (i=1; i<=2*nsimp+1; i++) {
		mofu=rsqrt(R*R/(a2+xarg[i])+z*z/(c2+xarg[i]));
		dpotdztmp=dpotdztmp+simpcoef[i]*z/(( rpow(1.0+mofu,3.0) * 
				(a2+xarg[i]) * rsqrt(c2+xarg[i]))*mofu*(c2+xarg[i]));
		mofu=rsqrt(R*R/(xarg6[i]*a2+1.0)+z*z/(xarg6[i]*c2+1.0));
		dpotdz1=dpotdz1+simpcoef[i]*xarg5[i]*z/( rpow(1.0+xarg3[i]*mofu,3.0)
			*(a2*xarg6[i]+1.0)*rsqrt(c2*xarg6[i]+1.0)*mofu*(c2*xarg6[i]+1.0));
	}
 
	*dpotdz = -0.5*G*M*(-2.*dpotdztmp-12.*dpotdz1)/(2.0*((real) nsimp)*3.0);
}

local void obdpotdR(real R, real z, realptr dpotdR)
{
//     Subroutine to compute partial-psi / partial-R.
 
	int i;
	real dpotdR1,mofu;
	static bool firstc = TRUE;
	static realptr simpcoef, xarg, xarg2, xarg3, xarg5, xarg6;
	static a2, c2;
	real dpotdRtmp;

	if (firstc) {
		firstc=FALSE;

		simpcoef=nr_dvector(1,nscoef);
		xarg=nr_dvector(1,nscoef);
		xarg2=nr_dvector(1,nscoef);
		xarg3=nr_dvector(1,nscoef);
		xarg5=nr_dvector(1,nscoef);
		xarg6=nr_dvector(1,nscoef);

		if (2*nsimp+1 > nscoef) error("\noverflow in oblatep\n");

           simpcoef[1]=1.0;

		for (i=2; i<=2*nsimp; i+=2) {
			simpcoef[i]=4.0;
			simpcoef[i+1]=2.0;
		}

		simpcoef[2*nsimp+1]=1.0;

		for (i=1; i<=2*nsimp+1; i++) {
			xarg[i]=(i-1.0)/(2.0*nsimp);
			xarg2[i]=rsqr(xarg[i]);
			xarg3[i]=xarg2[i]*xarg[i];
			xarg5[i]=xarg2[i]*xarg3[i];
			xarg6[i]=rpow(xarg2[i],3.0);
		}
 
		a2=a*a;
		c2=c*c;
	}
 
        dpotdRtmp=0.0;
        dpotdR1=0.0;
 
	for (i=1; i<=2*nsimp+1; i++) {
		mofu=rsqrt(R*R/(a2+xarg[i])+z*z/(c2+xarg[i]));
		dpotdRtmp=dpotdRtmp+simpcoef[i]*R/(( rpow(1.0+mofu,3.0) * 
				(a2+xarg[i]) * rsqrt(c2+xarg[i]))*mofu*(a2+xarg[i]));
		mofu=rsqrt(R*R/(xarg6[i]*a2+1.0)+z*z/(xarg6[i]*c2+1.0));
		dpotdR1=dpotdR1+simpcoef[i]*xarg5[i]*R/( rpow(1.0+xarg3[i]*mofu,3.0) *
			(a2*xarg6[i]+1.0)*rsqrt(c2*xarg6[i]+1.0)*
			mofu*(a2*xarg6[i]+1.0));
	}

	*dpotdR = -0.5*G*M*(-2.*dpotdRtmp-12.*dpotdR1)/(2.0*nsimp*3.0);
}


#undef nscoef

real cfunc(real xrotcirc, real xrad, real xkappa)
{
/*
Function defining the ratio of the squares of azimuthal and radial
velocity dispersion.                            
                                                                
OLD:                                                           
         A(x) = A(0) exp(-x) + 0.5                              
                                                                
where A(x)=<theta^2 /<phi^2  , A(0)=0.5 and x=R/h.
                                                                
NEW:                                                           
         cfunc(rotcirc,rad,kappa) = 0.25 * kappa^2/omega^2         
                                                                 
where  omega = rotcirc/rad                                  
                                                                
This function ensures a circular velocity ellipsoid in the center 
and an asymptotic behaviour like that predicted by the epicyclic 
approximation.                               
*/

	real xomega, cfunctmp;

	xomega=xrotcirc/xrad;

	if (xomega == 0.0 || xkappa == 0.0)
		cfunctmp=1.0;
	else
		cfunctmp=0.25*xkappa*xkappa/(xomega*xomega);
	return cfunctmp;
}

