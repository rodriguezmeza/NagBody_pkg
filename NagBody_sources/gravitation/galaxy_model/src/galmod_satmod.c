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


local void insmass(void);
local void satvel(void);

local void cmsat(void);
local void cmsatmod(void);

local void setsat(void);

local void forcebs(void);
local void forceds(void);
local void forcehs(void);

// Structure data and routines to compute forces and potential

#define NINTERP 30000

static acsmooth_st[NINTERP+1];
static phsmooth_st[NINTERP+1];

// ---------------------------------------------


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
 
//	for (i=0; i<=1+ninterp; i++) {
//		xw=i*deldrg;
//		xw2=xw*xw;
//		xw3=xw2*xw;
//		xw4=xw2*xw2;
//		if (xw <= 1.0) {
//			phsmooth[i]=-2.*xw3*(ONE/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.;
//			acsmooth[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
//		} else {
//			phsmooth[i]=-ONE/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.);
//			acsmooth[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;
//		}
//		if (xw >= 2.0) {
//			phsmooth[i]=ONE;
//			acsmooth[i]=1.0;
//		}
//	}


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

//	for (i=0; i<=1+ninterp; i++) {
//		xw=i*deldrg;
//		xw2=xw*xw;
//		xw3=xw2*xw;
//		xw4=xw2*xw2;
//		if (xw <= 1.0) {
//			phsmooth[i]=-2.*xw3*(ONE/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.;
//			acsmooth[i]=xw3*(4./3.-6.*xw2/5.+0.5*xw3);
//		} else {
//			phsmooth[i]=-ONE/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.);
//			acsmooth[i]=-1.0/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-xw4*xw2/6.;
//		}
//		if (xw >= 2.0) {
//			phsmooth[i]=ONE;
//			acsmooth[i]=1.0;
//		}
//	}
 
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


