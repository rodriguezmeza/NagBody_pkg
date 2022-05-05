/*==============================================================================
	MODULE: analysis_grav_io.c		[analysis_grav]
	Written by: Mario A. Rodriguez-Meza
	Starting date:	January, 2005
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:	stdinc.h, mathfns.h, vectmath
						vectmath.h, getparam.h
						types.h, stat.h, inout.h
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

local void PrintRhoTheta(string, string);
local void PrintRho(string, string);
local void PrintVcR(string, string);

//local void EvalVacf(void);
local void AccumVacf(string, string);
local void InitVacf(void);
local void ZeroVacf(void);
local void PrintVacf(string, string);

local void PrintSnap_Slab_long(string, string, int, char *);
local void PrintSnap_Slab_not_long(string, string, int, char *);

void startoutput(void)
{
	int nb;

    printf("\n%s\n%s: %s\n\t %s\n", gd.headline0, gd.headline1, gd.headline2, gd.headline3);

    if (! strnull(cmd.options))                     
        printf("\n\toptions: %s\n", cmd.options);

	gd.countRhoTheta = 0;
    gd.histRho0Theta = AllocVecR(cmd.sizeHistRhoTheta);
    gd.histRho1Theta = AllocVecR(cmd.sizeHistRhoTheta);
    gd.histRho2Theta = AllocVecR(cmd.sizeHistRhoTheta);

    gd.histRho0ThetaSave = AllocVecR(cmd.sizeHistRhoTheta);
    gd.histRho1ThetaSave = AllocVecR(cmd.sizeHistRhoTheta);
    gd.histRho2ThetaSave = AllocVecR(cmd.sizeHistRhoTheta);

    gd.histAccRTheta = AllocVecR(cmd.sizeHistRhoTheta);
    gd.histAccTTheta = AllocVecR(cmd.sizeHistRhoTheta);
    gd.histVelTTheta = AllocVecR(cmd.sizeHistRhoTheta);

	if (scanopt(cmd.data_analysis_type, "rho-anim")) {
		gd.countRho = 0;
		gd.histRho = AllocVecR(cmd.sizeHistRho);
		gd.histMass = AllocVecR(cmd.sizeHistRho);
	}

	if (scanopt(cmd.data_analysis_type, "vcr-anim")) {
		gd.countVcR = 0;
		gd.histVcR = AllocVecR(cmd.sizeHistVcR);
		gd.histMass = AllocVecR(cmd.sizeHistVcR);
		gd.histForce = AllocVecR(cmd.sizeHistVcR);
		gd.histNPart = AllocVecR(cmd.sizeHistVcR);
		gd.histRVelDisp = AllocVecR(cmd.sizeHistVcR);
		gd.histVrR = AllocVecR(cmd.sizeHistVcR);
	}

//	if (cmd.computeTransport) {
	if (scanopt(cmd.data_analysis_type, "acf-anim")) {
printf("\n\nAllocating memory for buffer,avAcfVel, ...\n");
		gd.nbody=400000;				// Reservar la cantidad apropiada de cuerpos... 
		AllocMem (gd.avAcfVel, cmd.nValAcf, real);
		AllocMem (gd.avAcfTherm, cmd.nValAcf, real);
		AllocMem (gd.avAcfVisc, cmd.nValAcf, real);
		AllocMem (tBuf, cmd.nBuffAcf, TBuf);
		for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
			AllocMem (tBuf[nb].acfVel, cmd.nValAcf, real);
			AllocMem (tBuf[nb].orgVel, gd.nbody, vector);
			AllocMem (tBuf[nb].acfTherm, cmd.nValAcf, real);
			AllocMem (tBuf[nb].acfVisc, cmd.nValAcf, real);
		}

//		if (strnull(cmd.restorefile))
		InitVacf();
	}

}

void EvalRhoTheta(string fname, string fnametmp)
{
    bodyptr p;
	real deltaTheta;
	real histSum;
	int i, l, j;
	double cpustart;
	real rangeTheta;
	real xi, yi, xf, yf, theta;
	real hdz;
	vector cmpos, tmpv;
	real tmass,r,accr,acct,velt;
	real RhoRMax;
	real histTheta0Sum, histTheta1Sum, histTheta2Sum;

	printf("\nEvalRhoTheta: Entrando ...");
	cpustart = cputime();                       

	hdz=0.5*cmd.RhoDeltaZ;
	rangeTheta = 2.0*PI;

	gd.countRhoTheta = gd.countRhoTheta + 1;
	if (gd.countRhoTheta == 1) {
		for (j = 1; j < cmd.sizeHistRhoTheta; j++) gd.histRho0Theta[j] = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++) gd.histRho1Theta[j] = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++) gd.histRho2Theta[j] = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++) gd.histAccRTheta[j] = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++) gd.histAccTTheta[j] = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++) gd.histVelTTheta[j] = 0.0;
	}
	deltaTheta = rangeTheta/cmd.sizeHistRhoTheta;

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		if ( Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0] ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nBulge CM: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	RhoRMax=0.4;
	xi=cmpos[0]; yi=cmpos[1]; 
	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) &&
			 r <= RhoRMax ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nDisk CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	gdtreegrav.rsize=1.01;
	gdtreegrav.usequad = TRUE;
	gdtreegrav.eps = 0.015;
	gdtreegrav.theta = 0.75;
	DO_BODY(p, bodytab, bodytab+gd.nbody) Update(p)=TRUE;
	maketree_grav(bodytab, gd.nbody, &gdtreegrav);
	fprintf(stdout,"\n\nncell, tdepth : %d %d\n",gdtreegrav.ncell,gdtreegrav.tdepth);

    gdtreegrav.nbbcalc = gdtreegrav.nbccalc = 0;   
	gdtreegrav.cpuforce=0.;
	xi=cmpos[0]; yi=cmpos[1];
fprintf(stdout,"\n\nhdz, RhoR, RhoDeltaR : %g %g %g\n",hdz,cmd.RhoR,cmd.RhoDeltaR);
fprintf(stdout,"ID's bulge: %d %d\n",gd.bodyIDMin[0],gd.bodyIDMax[0]);
fprintf(stdout,"ID's disk: %d %d\n",gd.bodyIDMin[1],gd.bodyIDMax[1]);
fprintf(stdout,"ID's halo: %d %d\n",gd.bodyIDMin[2],gd.bodyIDMax[2]);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( theta/deltaTheta ) +1;
			if (i>cmd.sizeHistRhoTheta) i=cmd.sizeHistRhoTheta;
				gd.histRho0Theta[i] = gd.histRho0Theta[i]+1.0;
			ind_normal_gravcalc(bodytab, gd.nbody, p, &gdtreegrav);
			acct = Acc(p)[0]*(-rsin(theta)) + Acc(p)[1]*rcos(theta);
			accr = Acc(p)[0]*rcos(theta) + Acc(p)[1]*rsin(theta);
			gd.histAccRTheta[i] = gd.histAccRTheta[i] + accr;
			gd.histAccTTheta[i] = gd.histAccTTheta[i] + acct;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVelTTheta[i] = gd.histVelTTheta[i] + velt;
		}
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( theta/deltaTheta ) +1;
			if (i>cmd.sizeHistRhoTheta) i=cmd.sizeHistRhoTheta;
				gd.histRho1Theta[i] = gd.histRho1Theta[i]+1.0;
			ind_normal_gravcalc(bodytab, gd.nbody, p, &gdtreegrav);
			acct = Acc(p)[0]*(-rsin(theta)) + Acc(p)[1]*rcos(theta);
			accr = Acc(p)[0]*rcos(theta) + Acc(p)[1]*rsin(theta);
			gd.histAccRTheta[i] = gd.histAccRTheta[i] + accr;
			gd.histAccTTheta[i] = gd.histAccTTheta[i] + acct;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVelTTheta[i] = gd.histVelTTheta[i] + velt;
		}
		if ( (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( theta/deltaTheta ) +1;
			if (i>cmd.sizeHistRhoTheta) i=cmd.sizeHistRhoTheta;
				gd.histRho2Theta[i] = gd.histRho2Theta[i]+1.0;
			ind_normal_gravcalc(bodytab, gd.nbody, p, &gdtreegrav);
			acct = Acc(p)[0]*(-rsin(theta)) + Acc(p)[1]*rcos(theta);
			accr = Acc(p)[0]*rcos(theta) + Acc(p)[1]*rsin(theta);
			gd.histAccRTheta[i] = gd.histAccRTheta[i] + accr;
			gd.histAccTTheta[i] = gd.histAccTTheta[i] + acct;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVelTTheta[i] = gd.histVelTTheta[i] + velt;
		}
	}

	if (gd.countRhoTheta == cmd.stepAvgRhoTheta) {
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++)
			histSum= histSum+gd.histRho0Theta[j];
		histTheta0Sum = histSum;
		for (j=1; j<=cmd.sizeHistRhoTheta; j++) {
			gd.histRho0ThetaSave[j]= gd.histRho0Theta[j];
			gd.histRho0Theta[j]= gd.histRho0Theta[j]/histSum;
		}

		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++)
			histSum= histSum+gd.histRho1Theta[j];
		histTheta1Sum = histSum;
		for (j=1; j<=cmd.sizeHistRhoTheta; j++) {
			gd.histRho1ThetaSave[j]= gd.histRho1Theta[j];
			gd.histRho1Theta[j]= gd.histRho1Theta[j]/histSum;
		}

		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++)
			histSum= histSum+gd.histRho2Theta[j];
		histTheta2Sum = histSum;
		for (j=1; j<=cmd.sizeHistRhoTheta; j++) {
			gd.histRho2ThetaSave[j]= gd.histRho2Theta[j];
			gd.histRho2Theta[j]= gd.histRho2Theta[j]/histSum;
		}

		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoTheta; j++)
			histSum= histSum+gd.histAccRTheta[j];
		for (j=1; j<=cmd.sizeHistRhoTheta; j++)
			gd.histAccRTheta[j]= gd.histAccRTheta[j]/histSum;

//		histSum = 0.0;
//		for (j = 1; j < cmd.sizeHistRhoTheta; j++) {
//			histSum= histSum+gd.histAccTTheta[j];
		histSum= histTheta0Sum + histTheta1Sum + histTheta0Sum;
		for (j=1; j<=cmd.sizeHistRhoTheta; j++)
//			gd.histAccTTheta[j]= gd.histAccTTheta[j]/histSum;
			gd.histAccTTheta[j]= gd.histAccTTheta[j]
			/(gd.histRho0ThetaSave[j]+gd.histRho1ThetaSave[j]+gd.histRho2ThetaSave[j]);

//		histSum = 0.0;
//		for (j = 1; j < cmd.sizeHistRhoTheta; j++)
//			histSum= histSum+gd.histVelTTheta[j];
		histSum= histTheta0Sum + histTheta1Sum + histTheta0Sum;
		for (j=1; j<=cmd.sizeHistRhoTheta; j++)
//			gd.histVelTTheta[j]= gd.histVelTTheta[j]/histSum;
			gd.histVelTTheta[j]= gd.histVelTTheta[j]/
			(gd.histRho0ThetaSave[j]+gd.histRho1ThetaSave[j]+gd.histRho2ThetaSave[j]);

		PrintRhoTheta(fname, fnametmp);
		gd.RhoTheta_flag = 1;
		gd.countRhoTheta = 0;
	}

	fprintf(stdout,"\ncputree, cpuforce : %g %g\n",gdtreegrav.cputree,gdtreegrav.cpuforce);
	fprintf(stdout,"nbbcalc, nbccalc : %d %d\n",gdtreegrav.nbbcalc,gdtreegrav.nbccalc);
	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

local void PrintRhoTheta(string fname, string fnametmp)
{
	real thetaBin; 
	int n;
    stream outstr;
    char   buf[200];
	real rangeTheta;

    outstr = stropen(fnametmp, "w!");

	rangeTheta = 2.0*PI;
	for (n=1; n<=cmd.sizeHistRhoTheta; n++) {
		thetaBin = (180.0/PI)*(n-0.5)*rangeTheta/cmd.sizeHistRhoTheta;
//		fprintf(outstr,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
		fprintf(outstr,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
//			thetaBin,gd.histRho0Theta[n],
			thetaBin,gd.histRho1Theta[n],
			thetaBin,gd.histRho2Theta[n],
//			thetaBin,gd.histAccRTheta[n],
			thetaBin,gd.histAccTTheta[n],
			thetaBin,gd.histVelTTheta[n]);
	}
    fclose(outstr);
	sprintf(buf,"mv %s %s",fnametmp,fname);
	printf("\nsystem: %s ...\n",buf);
	system(buf);
}


void EvalRhoDist(string fname, string fnametmp)
{
    bodyptr p;
	real deltaR;
	real histSum, normFac;
	int i, l, j;
	double cpustart;
	real xi, yi, zi, r;
	vector cmpos, tmpv;
	real tmass;

	printf("\nEvalRho: Entrando ...");
	cpustart = cputime();                       

	gd.countRho = gd.countRho + 1;
	if (gd.countRho == 1)
		for (j = 1; j <= cmd.sizeHistRho; j++) {
			gd.histRho[j] = 0.0;
			gd.histMass[j] = 0.0;
		}

	deltaR = cmd.rangeR/cmd.sizeHistRho;

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		tmass += Mass(p);
		MULVS(tmpv, Pos(p), Mass(p));           
		ADDV(cmpos, cmpos, tmpv);
	}
	DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"Initial CM: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi));
		if (r<=cmd.rangeR) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
	DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"Intermediate CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi));
		if (r<=cmd.rangeR/2.0) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
	DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"Final CM and range R : %g %g %g %g\n",cmpos[0],cmpos[1],cmpos[2],cmd.rangeR);
	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi));
		if ( r <= cmd.rangeR ) {
			i = (int) ( r/deltaR ) +1;
			if (i>cmd.sizeHistRho) i=cmd.sizeHistRho;
			gd.histRho[i] = gd.histRho[i]+1.0;
			gd.histMass[i] += Mass(p);
		}
	}

	if (gd.countRho == cmd.stepAvgRho) {
		normFac = 1. /	(4.0 * PI * rpow(deltaR, 3.0) );
		histSum = 0.0;
		for (j = 1; j <= cmd.sizeHistRho; j++)
			histSum= histSum+gd.histRho[j];
		for (j=1; j<=cmd.sizeHistRho; j++) {
			gd.histRho[j] = gd.histRho[j] * normFac / ( rsqr((double)j-0.5) * histSum);
			gd.histMass[j] = gd.histMass[j] * normFac / ( rsqr((double)j-0.5) );
		}

		PrintRho(fname, fnametmp);
		gd.Rho_flag = 1;
		gd.countRho = 0;
	}

	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

local void PrintRho(string fname, string fnametmp)
{
	real rBin; 
	int n;
    stream outstr;
    char   buf[200];

    outstr = stropen(fnametmp, "w!");

	for (n=1; n<=cmd.sizeHistRho; n++) {
		rBin = (n-0.5)*cmd.rangeR/cmd.sizeHistRho;
		fprintf(outstr,"%g %g %g\n",
			rBin,gd.histMass[n],gd.histRho[n]);
	}
    fclose(outstr);
	sprintf(buf,"mv %s %s",fnametmp,fname);
	printf("\nsystem: %s ...\n",buf);
	system(buf);
}

void EvalVcR(string fname, string fnametmp)
{
    bodyptr p;
	real deltaR;
	real histSum;
	int i, j;
	double cpustart;
	real xi, yi, zi, xf, yf;
	real theta;
	vector cmpos, tmpv;
	real tmass,r,velt,velr,amag;
	vector cmvel;
	real rxy, x, y, z, vx, vy, vz;
	real vxi, vyi, vzi, veltheta, velphi;

	printf("\nEvalVcR: Entrando ...");
	cpustart = cputime();                       

	gd.countVcR = gd.countVcR + 1;
	if (gd.countVcR == 1)
		for (j = 1; j <= cmd.sizeHistVcR; j++) {
			gd.histVcR[j] = 0.0;
			gd.histMass[j] = 0.0;
			gd.histForce[j] = 0.0;
			gd.histNPart[j] = 0.0;
			gd.histRVelDisp[j] = 0.0;
			gd.histVrR[j] = 0.0;
		}

	deltaR = cmd.rangeR/cmd.sizeHistVcR;


/*
     CLRV(cmvel);                                
...
       MULVS(tmpv, Pos(p), Mass(p));           
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           
        ADDV(cmvel, cmvel, tmpv);
    }
    etot[0] = etot[1] + etot[2];                
    DIVVS(cmpos, cmpos, mtot);                  
    DIVVS(cmvel, cmvel, mtot);
*/

	CLRV(cmpos);
	CLRV(cmvel);                                
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		tmass += Mass(p);
		MULVS(tmpv, Pos(p), Mass(p));           
		ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           
        ADDV(cmvel, cmvel, tmpv);
	}
	DIVVS(cmpos, cmpos, tmass);
    DIVVS(cmvel, cmvel, tmass);
	fprintf(stdout,"Initial CM: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
	fprintf(stdout,"Initial VCM: %g %g %g\n",cmvel[0],cmvel[1],cmvel[2]);
	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi));
		if (r<=cmd.rangeR) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
			MULVS(tmpv, Vel(p), Mass(p));           
			ADDV(cmvel, cmvel, tmpv);
		}
	}
	DIVVS(cmpos, cmpos, tmass);
    DIVVS(cmvel, cmvel, tmass);
	fprintf(stdout,"Intermediate CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
	fprintf(stdout,"Intermediate VCM : %g %g %g\n",cmvel[0],cmvel[1],cmvel[2]);
	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi));
		if (r<=cmd.rangeR/2.0) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
			MULVS(tmpv, Vel(p), Mass(p));           
			ADDV(cmvel, cmvel, tmpv);
		}
	}
	DIVVS(cmpos, cmpos, tmass);
    DIVVS(cmvel, cmvel, tmass);
	fprintf(stdout,"Final CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
	fprintf(stdout,"Final VCM : %g %g %g\n",cmvel[0],cmvel[1],cmvel[2]);
	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];
	vxi=cmvel[0]; vyi=cmvel[1]; vzi=cmvel[2];

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) 
				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi) );
		if ( r <= cmd.rangeR ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( r/deltaR ) +1;
			if (i>cmd.sizeHistVcR) i=cmd.sizeHistVcR;
//
			x = Pos(p)[0]-xi;
			y = Pos(p)[1]-yi;
			z = Pos(p)[2]-zi;
			vx = Vel(p)[0] - vxi;
			vy = Vel(p)[1] - vyi;
			vz = Vel(p)[2] - vzi;
			rxy = rsqrt( rsqr(x) + rsqr(y) + rsqr(z));
			velr = (vx*x/r) + (vy*y/r) + (vz*z/r);
			veltheta = (vx*z*x/(r*rxy)) + (vy*z*y/(r*rxy)) - (vz*rxy/r);
			velphi = -vx*y/rxy + vy*x/rxy;
			velt = rsqrt(rsqr(veltheta) + rsqr(velphi));
//
//			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
//			velr = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVcR[i] = gd.histVcR[i] + velt;
//			gd.histVcR[i] = gd.histVcR[i] + rabs(velt);
			ABSV(amag, Acc(p));
//			gd.histForce[i] = gd.histForce[i] + Mass(p)*amag;
			gd.histForce[i] = gd.histForce[i] + amag;
			gd.histNPart[i] += 1.0;
			gd.histMass[i] += Mass(p);
			gd.histRVelDisp[i] = gd.histRVelDisp[i] + velr*velr;
			gd.histVrR[i] = gd.histVrR[i] + velr;
		}
	}

	if (gd.countVcR == cmd.stepAvgVcR) {
//		histSum = 0.0;
//		for (j = 1; j <= cmd.sizeHistVcR; j++)
//			histSum= histSum+gd.histVcR[j];
		for (j=1; j<=cmd.sizeHistVcR; j++) {
			if (gd.histNPart != 0)
				gd.histVcR[j]= gd.histVcR[j]/gd.histNPart[j];
			else
				gd.histVcR[j]= 0.0;
		}
//		histSum = 0.0;
//		for (j = 1; j <= cmd.sizeHistVcR; j++)
//			histSum= histSum+gd.histForce[j];
		for (j=1; j<=cmd.sizeHistVcR; j++) {
			if (gd.histNPart != 0)
				gd.histForce[j]= gd.histForce[j]/gd.histNPart[j];
			else
				gd.histForce[j]= 0.0;
		}
		for (j=1; j<=cmd.sizeHistVcR; j++) {
			if (gd.histNPart != 0)
				gd.histRVelDisp[j]= gd.histRVelDisp[j]/gd.histNPart[j];
			else
				gd.histRVelDisp[j]= 0.0;
		}
		for (j=1; j<=cmd.sizeHistVcR; j++) {
			if (gd.histNPart != 0)
				gd.histVrR[j]= gd.histVrR[j]/gd.histNPart[j];
			else
				gd.histVrR[j]= 0.0;
		}
		PrintVcR(fname, fnametmp);
		gd.VcR_flag = 1;
		gd.countVcR = 0;
	}

	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

local void PrintVcR(string fname, string fnametmp)
{
	real RBin; 
	int n;
    stream outstr;
    char   buf[200];
	realptr histAcumMass;

    outstr = stropen(fnametmp, "w!");
    histAcumMass = AllocVecR(cmd.sizeHistVcR);

	histAcumMass[1]=gd.histMass[1];
	for (n=2; n<=cmd.sizeHistVcR; n++)
		histAcumMass[n]=histAcumMass[n-1]+gd.histMass[n];

	fprintf(outstr,"%s\n%s\n",
			"# r   | GM/r | <Vc> | rF  | <vr^2> | sigmar | NPart Hist | Mass Hist |",
			"# <1> | <2>  | <3>  | <4> | <5>    | <6>    | <7>        | <8>");
	for (n=1; n<=cmd.sizeHistVcR; n++) {
		RBin = (n-0.5)*cmd.rangeR/cmd.sizeHistVcR;
		fprintf(outstr,"%g %g %g %g %g %g %g %g\n",
			RBin,
			rsqrt(cmd.gravityConstant*histAcumMass[n]/RBin),
			gd.histVcR[n], 
			rsqrt(RBin*gd.histForce[n]), 
			gd.histRVelDisp[n],
			gd.histRVelDisp[n], //- rsqr(gd.histVrR[n]), 
			gd.histNPart[n], 
			gd.histMass[n] );
	}
    fclose(outstr);
	sprintf(buf,"mv %s %s",fnametmp,fname);
	printf("\nsystem: %s ...\n",buf);
	system(buf);
}

// COMIENZO DE LA EVALUACION DE LAS PROPIEDADES DE TRANSPORTE ------------------

void EvalAcf (string fname, string fnametmp)
{
	vector vecTherm, vecVisc;
	int n, nb, ni;
	bodyptr p;
	double cpustart;
	real xi, yi;
//	real xf, yf, theta;
	real hdz;
	vector cmpos, tmpv;
	real tmass,r,accr,acct,velt;
	real RhoRMax;
	int i;

	printf("EvalVacf: Entrando ... ");

	cpustart = cputime();                       

	hdz=0.5*cmd.RhoDeltaZ;


// COMIENZO EVALUACION DE LA POSICION DEL CENTRO DE MASA -----------------------
	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		if ( Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0] ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nBulge CM: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	RhoRMax=0.5;
	xi=cmpos[0]; yi=cmpos[1]; 
	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) &&
			 r <= RhoRMax ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nDisk CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
// FIN EVALUACION DE LA POSICION DEL CENTRO DE MASA ----------------------------

	gdtreegrav.rsize=1.01;
	gdtreegrav.usequad = TRUE;
	gdtreegrav.eps = 0.015;
	gdtreegrav.theta = 0.75;
	DO_BODY(p, bodytab, bodytab+gd.nbody) Update(p)=TRUE;
	maketree_grav(bodytab, gd.nbody, &gdtreegrav);
	fprintf(stdout,"\n\nncell, tdepth : %d %d\n",gdtreegrav.ncell,gdtreegrav.tdepth);

    gdtreegrav.nbbcalc = gdtreegrav.nbccalc = 0;   
	gdtreegrav.cpuforce=0.;
	xi=cmpos[0]; yi=cmpos[1];

fprintf(stdout,"\n\nhdz, RhoR, RhoDeltaR : %g %g %g\n",hdz,cmd.RhoR,cmd.RhoDeltaR);
fprintf(stdout,"ID's bulge: %d %d\n",gd.bodyIDMin[0],gd.bodyIDMax[0]);
fprintf(stdout,"ID's disk: %d %d\n",gd.bodyIDMin[1],gd.bodyIDMax[1]);
fprintf(stdout,"ID's halo: %d %d\n",gd.bodyIDMin[2],gd.bodyIDMax[2]);

// COMIENZO EVALUACION DE LA FUERZA, ENERGIA Y TENSORES DE PRESION -------------
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
//			xf=Pos(p)[0]; yf=Pos(p)[1];
//			theta = angle(xi, yi, xf, yf);
//			i = (int) ( theta/deltaTheta ) +1;
//			if (i>cmd.sizeHistRhoTheta) i=cmd.sizeHistRhoTheta;
//				gd.histRho0Theta[i] = gd.histRho0Theta[i]+1.0;
			ind_normal_gravcalc(bodytab, gd.nbody, p, &gdtreegrav);
//			acct = Acc(p)[0]*(-rsin(theta)) + Acc(p)[1]*rcos(theta);
//			accr = Acc(p)[0]*rcos(theta) + Acc(p)[1]*rsin(theta);
//			gd.histAccRTheta[i] = gd.histAccRTheta[i] + accr;
//			gd.histAccTTheta[i] = gd.histAccTTheta[i] + acct;
//			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
//			gd.histVelTTheta[i] = gd.histVelTTheta[i] + velt;
		}
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
/*
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( theta/deltaTheta ) +1;
			if (i>cmd.sizeHistRhoTheta) i=cmd.sizeHistRhoTheta;
				gd.histRho1Theta[i] = gd.histRho1Theta[i]+1.0;
*/
			ind_normal_gravcalc(bodytab, gd.nbody, p, &gdtreegrav);
/*
			acct = Acc(p)[0]*(-rsin(theta)) + Acc(p)[1]*rcos(theta);
			accr = Acc(p)[0]*rcos(theta) + Acc(p)[1]*rsin(theta);
			gd.histAccRTheta[i] = gd.histAccRTheta[i] + accr;
			gd.histAccTTheta[i] = gd.histAccTTheta[i] + acct;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVelTTheta[i] = gd.histVelTTheta[i] + velt;
*/
		}
		if ( (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
/*
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( theta/deltaTheta ) +1;
			if (i>cmd.sizeHistRhoTheta) i=cmd.sizeHistRhoTheta;
				gd.histRho2Theta[i] = gd.histRho2Theta[i]+1.0;
*/
			ind_normal_gravcalc(bodytab, gd.nbody, p, &gdtreegrav);
/*
			acct = Acc(p)[0]*(-rsin(theta)) + Acc(p)[1]*rcos(theta);
			accr = Acc(p)[0]*rcos(theta) + Acc(p)[1]*rsin(theta);
			gd.histAccRTheta[i] = gd.histAccRTheta[i] + accr;
			gd.histAccTTheta[i] = gd.histAccTTheta[i] + acct;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVelTTheta[i] = gd.histVelTTheta[i] + velt;
*/
		}
	}
	fprintf(stdout,"\ncputree, cpuforce : %g %g\n",gdtreegrav.cputree,gdtreegrav.cpuforce);
	fprintf(stdout,"nbbcalc, nbccalc : %d %d\n",gdtreegrav.nbbcalc,gdtreegrav.nbccalc);
// FIN EVALUACION DE LA FUERZA, ENERGIA Y TENSORES DE PRESION ------------------

	VZero (vecVisc);
	VZero (vecTherm);
//printf("\n\nAqui voy [1]\n");
	DO_BODY(p, bodytab, bodytab+gd.nbody) {	// HACER LA VERSION 2D ...
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			vecVisc[0] += Vel(p)[1] * Vel(p)[2] + 0.5 * rf(p)[1][2];
			vecVisc[1] += Vel(p)[2] * Vel(p)[0] + 0.5 * rf(p)[2][0];
			vecVisc[2] += Vel(p)[0] * Vel(p)[1] + 0.5 * rf(p)[0][1];
			en(p) += VLenSq (Vel(p));
			VVSAdd (vecTherm, 0.5 * en(p), Vel(p));
			vecTherm[0] += 0.5 * VDot(Vel(p), rf(p)[0]);
			vecTherm[1] += 0.5 * VDot(Vel(p), rf(p)[1]);
			vecTherm[2] += 0.5 * VDot(Vel(p), rf(p)[2]);
		}
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			vecVisc[0] += Vel(p)[1] * Vel(p)[2] + 0.5 * rf(p)[1][2];
			vecVisc[1] += Vel(p)[2] * Vel(p)[0] + 0.5 * rf(p)[2][0];
			vecVisc[2] += Vel(p)[0] * Vel(p)[1] + 0.5 * rf(p)[0][1];
			en(p) += VLenSq (Vel(p));
			VVSAdd (vecTherm, 0.5 * en(p), Vel(p));
			vecTherm[0] += 0.5 * VDot(Vel(p), rf(p)[0]);
			vecTherm[1] += 0.5 * VDot(Vel(p), rf(p)[1]);
			vecTherm[2] += 0.5 * VDot(Vel(p), rf(p)[2]);
		}
		if ( (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			vecVisc[0] += Vel(p)[1] * Vel(p)[2] + 0.5 * rf(p)[1][2];
			vecVisc[1] += Vel(p)[2] * Vel(p)[0] + 0.5 * rf(p)[2][0];
			vecVisc[2] += Vel(p)[0] * Vel(p)[1] + 0.5 * rf(p)[0][1];
			en(p) += VLenSq (Vel(p));
			VVSAdd (vecTherm, 0.5 * en(p), Vel(p));
			vecTherm[0] += 0.5 * VDot(Vel(p), rf(p)[0]);
			vecTherm[1] += 0.5 * VDot(Vel(p), rf(p)[1]);
			vecTherm[2] += 0.5 * VDot(Vel(p), rf(p)[2]);
		}
	}
//printf("\n\nAqui voy [2]\n");

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBuf[nb].count == 0) {
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
				if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) 
						&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
						&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
//printf("Aqui voy [0] : %d\n",p-bodytab+1);
					n = p-bodytab;
					SETV(tBuf[nb].orgVel[n], Vel(p));
				}
				if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) 
						&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
						&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
//printf("Aqui voy [1] : %d\n",p-bodytab+1);
					n = p-bodytab;
					SETV(tBuf[nb].orgVel[n], Vel(p));
				}
				if ( (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2]) 
						&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
						&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
//printf("Aqui voy [2] : %d\n",p-bodytab+1);
					n = p-bodytab;
					SETV(tBuf[nb].orgVel[n], Vel(p));
				}
			}
		}
//printf("\n\nAqui voy [3]\n");
		if (tBuf[nb].count >= 0) {
			ni = tBuf[nb].count;
			tBuf[nb].acfVel[ni] = 0.;
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
				if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) 
						&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
						&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
					n = p-bodytab;
					tBuf[nb].acfVel[ni] += VDot(tBuf[nb].orgVel[n], Vel(p));
				}
				if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) 
						&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
						&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
					n = p-bodytab;
					tBuf[nb].acfVel[ni] += VDot(tBuf[nb].orgVel[n], Vel(p));
				}
				if ( (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2]) 
						&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
						&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
					n = p-bodytab;
					tBuf[nb].acfVel[ni] += VDot(tBuf[nb].orgVel[n], Vel(p));
				}
			}
		}
//printf("\n\nAqui voy [4]\n");
		if (tBuf[nb].count == 0) {
			SETV(tBuf[nb].orgVisc, vecVisc);
			SETV(tBuf[nb].orgTherm, vecTherm);
		}
		tBuf[nb].acfVisc[ni] = VDot(tBuf[nb].orgVisc, vecVisc);
		tBuf[nb].acfTherm[ni] = VDot(tBuf[nb].orgTherm, vecTherm);
		++ tBuf[nb].count;
	}
//printf("\n\nAqui voy [5]\n");
	AccumVacf(fname, fnametmp);
//printf("\n\nAqui voy [6]\n");

	printf("Saliendo : EvalACF CPU time = %lf\n",cputime()-cpustart); 
}

local void AccumVacf(string fname, string fnametmp)
{
	real fac;
	int j, nb;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBuf[nb].count == cmd.nValAcf) {
			for (j = 0; j < cmd.nValAcf; j ++) gd.avAcfVel[j] += tBuf[nb].acfVel[j];
			for (j = 0; j < cmd.nValAcf; j ++) {
				gd.avAcfVisc[j] += tBuf[nb].acfVisc[j];
				gd.avAcfTherm[j] += tBuf[nb].acfTherm[j];
			}
			tBuf[nb].count = 0;
			++ gd.countAcfAv;
			if (gd.countAcfAv == cmd.limitAcfAv) {
//				fac = cmd.stepAcf * gd.dtime / (NDIM * gd.nbody * cmd.limitAcfAv);
				fac = 1.0 / (NDIM * cmd.limitAcfAv);
				gd.intAcfVel = fac * Integrate (gd.avAcfVel, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) gd.avAcfVel[j] /= gd.avAcfVel[0];
				gd.avAcfVel[0] = 1.;
//				fac = cmd.density * cmd.stepAcf * gd.dtime /
//						(3. * cmd.temperature * gd.nbody * cmd.limitAcfAv);
				gd.intAcfVisc = fac * Integrate (gd.avAcfVisc, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) gd.avAcfVisc[j] /= gd.avAcfVisc[0];
				gd.avAcfVisc[0] = 1.;
//				fac = cmd.density * cmd.stepAcf * gd.dtime /
//						(3. * rsqr (cmd.temperature) * gd.nbody * cmd.limitAcfAv);
				gd.intAcfTherm = fac * Integrate (gd.avAcfTherm, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) gd.avAcfTherm[j] /= gd.avAcfTherm[0];
				gd.avAcfTherm[0] = 1.;
				PrintVacf (fname, fnametmp);
				ZeroVacf ();
			}
		}
	}
}

local void InitVacf(void)
{
	int nb;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++)
		tBuf[nb].count = - nb * cmd.nValAcf / cmd.nBuffAcf;
	ZeroVacf ();
}

local void ZeroVacf (void)
{
	int j;

	gd.countAcfAv = 0;
	for (j = 0; j < cmd.nValAcf; j ++) gd.avAcfVel[j] = 0.;
	for (j = 0; j < cmd.nValAcf; j ++) {
		gd.avAcfVisc[j] = 0.;
		gd.avAcfTherm[j] = 0.;
	}
}

//#define saveacftmp "acf.tmp"
//#define saveacf "acf.dat"

local void PrintVacf(string fname, string fnametmp)
{
	real tVal;
	int j;
    stream outstr;
    char   buf[200];

    outstr = stropen(fnametmp, "w!");

	fprintf (outstr, "acf\n");
	for (j = 0; j < cmd.nValAcf; j ++) {
//		tVal = j * cmd.stepAcf * gd.dtime;
		tVal = j;
		fprintf (outstr, "%8.4f %8.4f %8.4f %8.4f %8.4f\n", tVal, gd.tnow,
		gd.avAcfVel[j], gd.avAcfVisc[j], gd.avAcfTherm[j]);
	}
	fprintf (outstr, "acf integrals: %8.3f %8.3f %8.3f\n",
		gd.intAcfVel, gd.intAcfVisc, gd.intAcfTherm);

    fclose(outstr);
	sprintf(buf,"mv %s %s",fnametmp,fname);
//	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);
}

//#undef saveacftmp
//#undef saveacf

// FIN DE LA EVALUACION DE LAS PROPIEDADES DE TRANSPORTE -----------------------

void LocateBodiesID(string fname, string fnametmp)
{
    bodyptr p;
	double cpustart;
	real xi, yi, zi, xf, yf, theta;
	real hdz;
	vector cmpos, tmpv;
	real tmass,r;
	real RhoRMax;
    stream outstr;
    char   buf[200];

	printf("\nLocateBodiesID: Entrando ...");
	cpustart = cputime();

    outstr = stropen(fnametmp, "w!");

	hdz=0.5*cmd.RhoDeltaZ;

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		if ( Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0] ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nBulge CM: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	RhoRMax=0.4;
	xi=cmpos[0]; yi=cmpos[1]; 
	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) &&
			 r <= RhoRMax ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nDisk CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];
	fprintf(stdout,"\n\nhdz, RhoR, RhoDeltaR : %g %g %g\n",hdz,cmd.RhoR,cmd.RhoDeltaR);
	fprintf(stdout,"ID's bulge: %d %d\n",gd.bodyIDMin[0],gd.bodyIDMax[0]);
	fprintf(stdout,"ID's disk: %d %d\n",gd.bodyIDMin[1],gd.bodyIDMax[1]);
	fprintf(stdout,"ID's halo: %d %d\n",gd.bodyIDMin[2],gd.bodyIDMax[2]);

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			if (theta>=cmd.ThetaMin && theta<=cmd.ThetaMax)
				fprintf(outstr,"%d %g %g %g\n",Id(p),
					Pos(p)[0]-xi,Pos(p)[1]-yi,Pos(p)[2]-zi);
		}
	}
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			if (theta>=cmd.ThetaMin && theta<=cmd.ThetaMax)
				fprintf(outstr,"%d %g %g %g\n",Id(p),
					Pos(p)[0]-xi,Pos(p)[1]-yi,Pos(p)[2]-zi);
		}
	}
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz)
				&& (r>=cmd.RhoR && r<=cmd.RhoR+cmd.RhoDeltaR) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			if (theta>=cmd.ThetaMin && theta<=cmd.ThetaMax)
				fprintf(outstr,"%d %g %g %g\n",Id(p),
					Pos(p)[0]-xi,Pos(p)[1]-yi,Pos(p)[2]-zi);
		}
	}

    fclose(outstr);
	sprintf(buf,"mv %s %s",fnametmp,fname);
	printf("\nsystem: %s ...\n",buf);
	system(buf);
	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define massesfile	"./groups_catalog/masses_dist.dat"
#define outsnapgroup	"./groups_catalog/snapgroup%03d"

void GroupsCatalog(int snapcount)
{
    char namebuf[256];
//    stream outstr;
//    bodyptr p;
    struct stat buf;
    stream instr, outstr;
    bodyptr p;
	int i, ngroups, *tpgroup, *seq, *ig;
	realptr gmass;

	printf("\n\nGroups catalog process ...");

//printf("\n\nAqui voy [2]\n");
    if (stat(cmd.foffile, &buf) != 0) {
		error("\n\nfof file %s does not exist\n",cmd.foffile);
    } else {
        instr = stropen(cmd.foffile, "r");
	}

	in_int(instr, &ngroups);
    gmass = (realptr) allocate(ngroups*sizeof(real));
    tpgroup = (int *) allocate(ngroups*sizeof(int));
    seq = (int *) allocate(ngroups*sizeof(int));
    ig = (int *) allocate(ngroups*sizeof(int));
	DO_BODY(p, bodytab, bodytab+gd.nbody)
		in_int(instr, &IdG(p));
    fclose(instr);

	for (i=1; i<=ngroups; i++) {
		gmass[i-1] = 0.;
		tpgroup[i-1] = 0;
		ig[i-1] = i;
		DO_BODY(p, bodytab, bodytab+gd.nbody) {
			if (i == IdG(p)) {
				gmass[i-1] += Mass(p);
				tpgroup[i-1] += 1;
			}
		}
	}

//printf("\n\nAqui voy [3]\n");
/*
	outstr = stropen(massesfile, "w!");
	for (i=1; i<=ngroups; i++)
		fprintf(outstr,"%d\t%d\t%g\n",i,tpgroup[i-1],gmass[i-1]);
	fclose(outstr);
*/

//    sprintf(namebuf, out, snapcount);           
//    outstr = stropen(namebuf, "w!");
	outstr = stropen(massesfile, "w!");

//	HeapSortDescend(&gmass[0], &seq[0], ngroups);
//	for (i=1; i<=ngroups; i++) {
	HeapSort(&gmass[0], &seq[0], ngroups);
	for (i=ngroups; i>=1; i--) {
		fprintf(outstr,"%d\t%d\t%d\t%g\n",ngroups-i+1, ig[seq[i-1]], tpgroup[seq[i-1]],gmass[seq[i-1]]);
	}
	fclose(outstr);

//printf("\n\nAqui voy [4]\n");

	for (i=ngroups; i>=1; i--) {
		sprintf(namebuf, outsnapgroup, ngroups-i+1);
		outstr = stropen(namebuf, "w!");
//		out_int(outstr, tpgroup[seq[i-1]]);
//		out_int(outstr, NDIM);
//		out_real(outstr, tnow);
		fprintf(outstr,"#   nbody NDIM time\n# %d %d %lf\n",tpgroup[seq[i-1]],NDIM,gd.tnow);
		DO_BODY(p, bodytab, bodytab+gd.nbody) {
			if (ig[seq[i-1]] == IdG(p)) {
//				out_real(outstr, Mass(p));              
//				out_vector(outstr, Pos(p));
//				out_vector(outstr, Vel(p));
				out_int_mar(outstr, Id(p));
				out_real_mar(outstr, Mass(p));
// Apply boundary condition -------------------
				out_vector_mar(outstr, Pos(p));             
// --------------------------------------------
				out_vector_mar(outstr, Vel(p));
				fprintf(outstr,"\n");
			}
		}
		fclose(outstr);
	}

//printf("\n\nAqui voy [5]\n");

	printf(" done.\n");
}

#undef massesfile
#undef outsnapgroup

void PrintBodiesSets_Spherical(string savesnap, string savesnaptmp, int nbody, 
	int nbodiesSets, int *bodyIDMin, int *bodyIDMax, real RMax, char *options)
{
    stream outstr;
    bodyptr p;
	int i, id;
    char   buf[200];
// Comienzo de inclusion para ver disco...
	vector cmpos, tmpv;
	real xi, yi, zi;
	real tmass,r;
//	real RhoRMax;
// Fin de inclusion para ver disco...

//printf("\nAqui voy dentro de PrintBodiesSets... (1)\n");

	for (i=0; i<nbodiesSets; i++)
		if (bodyIDMin[i]>nbody || bodyIDMax[i]>nbody)
			error("\n\nPrintBodiesSets: IDMin or IDMax out of range\n");

    outstr = stropen(savesnaptmp, "w!");

	CLRV(cmpos);
	tmass=0.;

    printf("IDMin, IDMax: %d %d\n",bodyIDMin[0],bodyIDMax[0]);

	DO_BODY(p, bodytab, bodytab+nbody) {
		if ( Id(p)>=bodyIDMin[0] && Id(p)<=bodyIDMax[0] ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nCM first set: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

/*
	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];
	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi) );
		if ( (Id(p)>=bodyIDMin[1] && Id(p)<=bodyIDMax[1]) ) {
//        if ( (Id(p)>=bodyIDMin[1] && Id(p)<=bodyIDMax[1]) &&
//            r <= RMax ) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	fprintf(stdout,"\n\nDisk CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];
*/

printf(">>>Selecting nbodiesSets: %d\n", nbodiesSets);

	DO_BODY(p, bodytab, bodytab+nbody) {
//		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
//				+(Pos(p)[2]-zi)*(Pos(p)[2]-zi) );
//		if (r<=RMax) {
		for (i=0; i<nbodiesSets; i++) {
			if (Id(p)>=bodyIDMin[i] && Id(p)<=bodyIDMax[i]) {
//			if (Id(p)>=bodyIDMin[1] && Id(p)<=bodyIDMax[1]) { // Ver disco...
				out_vector_mar(outstr, Pos(p));
				out_vector_mar(outstr, Vel(p));
				out_int_mar(outstr, Id(p));
				out_int_mar(outstr, Type(p));
				if (scanopt(options, "out-phi"))
					out_real_mar(outstr, Phi(p));
				if (scanopt(options, "out-acc"))
					out_vector_mar(outstr, Acc(p));
				fprintf(outstr,"\n");
			}
		}
//		}
    }

    
/*
// ALTERNATIVE ADDED CODE LINES ************************************************

	DO_BODY(p, bodytab, bodytab+nbody) {
		if (scanopt(options, "wr-cm")) {
			r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
					+(Pos(p)[2]-zi)*(Pos(p)[2]-zi) );
			if (r<=RMax) {
			for (i=0; i<nbodiesSets; i++) {
//			if (Id(p)>=bodyIDMin[i] && Id(p)<=bodyIDMax[i]) {
			if (Id(p)>=bodyIDMin[1] && Id(p)<=bodyIDMax[1]) { // Ver disco...
				tmpv[0]=Pos(p)[0]-xi;
				tmpv[1]=Pos(p)[1]-yi;
				tmpv[2]=Pos(p)[2]-zi;
				out_vector_mar(outstr, &tmpv[0]);
				out_vector_mar(outstr, Vel(p));
				out_int_mar(outstr, Id(p));
				out_int_mar(outstr, Type(p));
				if (scanopt(options, "out-phi"))
					out_real_mar(outstr, Phi(p));
				if (scanopt(options, "out-acc"))
					out_vector_mar(outstr, Acc(p));
				fprintf(outstr,"\n");
			}
			}
			}
		} else {
			out_vector_mar(outstr, Pos(p));             
			out_vector_mar(outstr, Vel(p));
			out_int_mar(outstr, Id(p));
			out_int_mar(outstr, Type(p));
			if (scanopt(options, "out-phi"))
				out_real_mar(outstr, Phi(p));
			if (scanopt(options, "out-acc"))
				out_vector_mar(outstr, Acc(p));
			fprintf(outstr,"\n");
		}
    }

// END ALTERNATIVE ADDED CODE LINES ********************************************
*/
    fclose(outstr);
	sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
	printf("\nsystem: %s",buf);
	system(buf);

//printf("\nAqui voy dentro de PrintBodiesSets... (2)\n");

}

void PrintSnap_Spherical(string savesnap, string savesnaptmp, int nbody, char *options)
{
    stream outstr_snap;
    bodyptr p;
    char   buf[200];
	vector cmpos, tmpv;
	real xi, yi, zi;
	real tmass,r;

    outstr_snap = stropen(savesnaptmp, "w!");

	fprintf(stdout,"In PrintSnap_Spherical\n");

	if (scanopt(options, "wr-cm")) {
		CLRV(cmpos);
		tmass=0.;
		DO_BODY(p, bodytab, bodytab+nbody) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
		DIVVS(cmpos, cmpos, tmass);
		fprintf(stdout,"Initial CM: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
		xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];

		CLRV(cmpos);
		tmass=0.;
		DO_BODY(p, bodytab, bodytab+nbody) {
			r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
					+(Pos(p)[2]-zi)*(Pos(p)[2]-zi));
			if (r<=cmd.rangeR) {
				tmass += Mass(p);
				MULVS(tmpv, Pos(p), Mass(p));           
				ADDV(cmpos, cmpos, tmpv);
			}
		}
		DIVVS(cmpos, cmpos, tmass);
		fprintf(stdout,"Intermediate CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
		xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];

		CLRV(cmpos);
		tmass=0.;
		DO_BODY(p, bodytab, bodytab+nbody) {
			r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
					+(Pos(p)[2]-zi)*(Pos(p)[2]-zi));
			if (r<=cmd.rangeR/2.0) {
				tmass += Mass(p);
				MULVS(tmpv, Pos(p), Mass(p));           
				ADDV(cmpos, cmpos, tmpv);
			}
		}
		DIVVS(cmpos, cmpos, tmass);
		fprintf(stdout,"Final CM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
		xi=cmpos[0]; yi=cmpos[1]; zi=cmpos[2];
	}

	DO_BODY(p, bodytab, bodytab+nbody) {
		if (scanopt(options, "wr-cm")) {
			r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi)
					+(Pos(p)[2]-zi)*(Pos(p)[2]-zi) );
			if (r<=cmd.rangeR) {
				tmpv[0]=Pos(p)[0]-xi;
				tmpv[1]=Pos(p)[1]-yi;
				tmpv[2]=Pos(p)[2]-zi;
				out_vector_mar(outstr_snap, &tmpv[0]);
				out_vector_mar(outstr_snap, Vel(p));
				out_int_mar(outstr_snap, Id(p));
				out_int_mar(outstr_snap, Type(p));
				if (scanopt(options, "out-phi"))
					out_real_mar(outstr_snap, Phi(p));
				if (scanopt(options, "out-acc"))
					out_vector_mar(outstr_snap, Acc(p));
				fprintf(outstr_snap,"\n");
			}
		} else {
			out_vector_mar(outstr_snap, Pos(p));             
			out_vector_mar(outstr_snap, Vel(p));
			out_int_mar(outstr_snap, Id(p));
			out_int_mar(outstr_snap, Type(p));
			if (scanopt(options, "out-phi"))
				out_real_mar(outstr_snap, Phi(p));
			if (scanopt(options, "out-acc"))
				out_vector_mar(outstr_snap, Acc(p));
			fprintf(outstr_snap,"\n");
		}
    }
    fclose(outstr_snap);
	sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
	system(buf);
}

void PrintSnap_Slab(string savesnap, string savesnaptmp, int nbody, char *options)
{
	if (gd.in_long_fmt==1) {
		PrintSnap_Slab_long(savesnap, savesnaptmp, nbody, options);
	} else {
		PrintSnap_Slab_not_long(savesnap, savesnaptmp, nbody, options);
	}	
}

local void PrintSnap_Slab_long(string savesnap, string savesnaptmp, int nbody, char *options)
{
    stream outstr_snap;
    bodyptr_long p;
    char   buf[200];
	vector cmpos, tmpv;
	real xi, yi, zi;
	real tmass,r;

    outstr_snap = stropen(savesnaptmp, "w!");

	fprintf(stdout,"In PrintSnap_Slab_long\n");

	DO_BODY(p, bodytab_long, bodytab_long+nbody) {
		if ( Pos_long(p)[0] >= cmd.xmin && Pos_long(p)[0] <= cmd.xmax &&
			 Pos_long(p)[1] >= cmd.ymin && Pos_long(p)[1] <= cmd.ymax &&
			 Pos_long(p)[2] >= cmd.zmin && Pos_long(p)[2] <= cmd.zmax ) {
			out_vector_mar(outstr_snap, Pos_long(p));             
			out_vector_mar(outstr_snap, Vel_long(p));
			out_int_mar(outstr_snap, Id_long(p));
			out_int_mar(outstr_snap, Type_long(p));
//			if (scanopt(options, "out-phi"))
//				out_real_mar(outstr_snap, Phi_long(p));
//			if (scanopt(options, "out-acc"))
//				out_vector_mar(outstr_snap, Acc_long(p));
			fprintf(outstr_snap,"\n");
		}
    }
    fclose(outstr_snap);
	sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
	system(buf);
}

local void PrintSnap_Slab_not_long(string savesnap, string savesnaptmp, int nbody, char *options)
{
    stream outstr_snap;
    bodyptr p;
    char   buf[200];
	vector cmpos, tmpv;
	real xi, yi, zi;
	real tmass,r;
	
    outstr_snap = stropen(savesnaptmp, "w!");
	
	fprintf(stdout,"In PrintSnap_Slab_not_long\n");
	
	DO_BODY(p, bodytab, bodytab+nbody) {
		if ( Pos(p)[0] >= cmd.xmin && Pos(p)[0] <= cmd.xmax &&
			Pos(p)[1] >= cmd.ymin && Pos(p)[1] <= cmd.ymax &&
			Pos(p)[2] >= cmd.zmin && Pos(p)[2] <= cmd.zmax ) {
			out_vector_mar(outstr_snap, Pos(p));             
			out_vector_mar(outstr_snap, Vel(p));
			out_int_mar(outstr_snap, Id(p));
			out_int_mar(outstr_snap, Type(p));
			if (scanopt(options, "out-phi"))
				out_real_mar(outstr_snap, Phi(p));
			if (scanopt(options, "out-acc"))
				out_vector_mar(outstr_snap, Acc(p));
			fprintf(outstr_snap,"\n");
		}
    }
    fclose(outstr_snap);
	sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
	system(buf);
}

