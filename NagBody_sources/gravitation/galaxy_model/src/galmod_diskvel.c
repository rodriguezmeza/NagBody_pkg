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


local void sigmap(void);
local void sigmar(void);
local void sigmaz(void);

local void circv(void);

local void radacc(void);


local void forceb(void);
local void forced(void);
local void forced_direct_original(void);
local void forced_direct_new(void);
local void forced_tree(void);
local void forceh(void);

local void setmesh(int nphi, realptr, realptr, realptr);
local void setrot(void);

local void setsigma(void);
local void sigalar(void);
local void sigcheck(void);

local void meanrot(void);
local void meshacc(int nphi, realptr, realptr, realptr, realptr);

local void dadr(void);
local void getkappa(void);

// Structure data and routines to compute forces and potential

#define NINTERP 30000

static acsmooth_st[NINTERP+1];
static phsmooth_st[NINTERP+1];

local void set_interpolation_force_table(void);

// ---------------------------------------------

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
// LIBERAR LA MEMORIA (bodytab)
    free(bodytab);
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

