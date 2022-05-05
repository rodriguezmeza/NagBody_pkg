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


local void bulgepot(void);

// Structure data and routines to compute forces and potential

#define NINTERP 30000

static acsmooth_st[NINTERP+1];
static phsmooth_st[NINTERP+1];

// ---------------------------------------------


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

#undef NINTERP




