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



local void halopot(void);

// Structure data and routines to compute forces and potential

#define NINTERP 30000

static acsmooth_st[NINTERP+1];
static phsmooth_st[NINTERP+1];

//local void set_interpolation_force_table(void);

// ---------------------------------------------



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

// Mar debugging...
    int micontador; //BORRAR

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

    printf("\n Aqui voy [1]:: Antes de entrar en halopot\n");

	halopot();

    printf("\n Aqui voy [2]:: Saliendo de halopot\n");
 
	ntest=gd.nbodies-cmd.nhalo+1-cmd.nsat;
// Mar debugging...
fprintf(stdout,"\n Antes de CIENTOSESENTA:: %d, %d, %d, %d, %d, %d",
        ntest, gd.nbodies, cmd.nhalo, cmd.nbulge, cmd.ndstars, cmd.nsat);
//    pause;
 
	vkludge=0.95;

// Mar debugging...
    
    micontador = 0;

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
    } else {
// Mar debugging...
        fprintf(stdout,"\n Yendo a CIENTOSESENTA:: %d %d", ++micontador, ntest);
		goto CIENTOSESENTA;
    }

    if (ntest <= gd.nbodies-cmd.nsat) {
// Mar debugging...
        fprintf(stdout,"\n Yendo a CIENTOSESENTA:: %d %d", ++micontador, ntest);
        goto CIENTOSESENTA;
    }

// Mar debugging...
fprintf(stdout,"\n Antes de CIENTOSESENTA:: %d, %d, %d, %d, %d, %d",
            ntest, gd.nbodies, cmd.nhalo, cmd.nbulge, cmd.ndstars, cmd.nsat);

    free_dvector(xintbin,0,maxnbin+1000);
	free_dvector(xmbin,0,maxnbin+1000);
	free_dvector(gridmass,0,maxnbin+1000);

	fprintf(stdout,"\n\nHaloVel CPU time: %g\n",cputime()-cpustart);
	fflush(stdout);
}

#undef NINTERP
