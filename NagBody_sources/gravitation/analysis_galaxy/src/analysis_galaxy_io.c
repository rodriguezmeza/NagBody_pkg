/*==============================================================================
	MODULE: analysis_galaxy_io.c		[datanaly_galaxy]
	Written by: M.A. Rodriguez-Meza
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
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/mathfns.h"
//#include "../../../General_libs/math/mathutil.h"
//#include "../../../General_libs/math/vectmath.h"
//#include "../../../General_libs/general/getparam.h"
//#include "../../../General_libs/io/inout.h"
//#include "../../../General_libs/general/constant.h"

//#include "../../../General_libs/NagBody/nagbody.h"

#include "globaldefs.h"
#include "protodefs.h"

//#include <sys/types.h>
//#include <sys/stat.h>

/* #include <strings.h> */							// For unix
//#include "../../../General_libs/general/strings.h"	// For Visual C

local void PrintRhoTheta(string, string);
local void PrintVcR(string, string);

//local void EvalVacf(void);
local void AccumVacf(string, string);
local void InitVacf(void);
local void ZeroVacf(void);
local void PrintVacf(string, string);

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

	gd.countVcR = 0;
    gd.histVcR = AllocVecR(cmd.sizeHistVcR);

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

void EvalVcR(string fname, string fnametmp)
{
    bodyptr p;
	real deltaR;
	real histSum;
	int i, j;
	double cpustart;
	real xi, yi, xf, yf;
	real theta;
	real hdz;
	vector cmpos, tmpv;
	real tmass,r,velt;

	printf("\nEvalVcR: Entrando ...");
	cpustart = cputime();                       

	hdz=0.5*cmd.RhoDeltaZ;

	gd.countVcR = gd.countVcR + 1;
	if (gd.countVcR == 1) {
		for (j = 1; j < cmd.sizeHistVcR; j++) gd.histVcR[j] = 0.0;
	}
	deltaR = cmd.rangeR/cmd.sizeHistVcR;

	CLRV(cmpos);
	tmass=0.;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) ||
			 (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) ||
			 (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2])) {
			tmass += Mass(p);
			MULVS(tmpv, Pos(p), Mass(p));           
			ADDV(cmpos, cmpos, tmpv);
		}
	}
    DIVVS(cmpos, cmpos, tmass);
	printf("\n\nCM : %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);

	xi=cmpos[0]; yi=cmpos[1];

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		r=rsqrt( (Pos(p)[0]-xi)*(Pos(p)[0]-xi)+(Pos(p)[1]-yi)*(Pos(p)[1]-yi) );
		if ( (Id(p)>=gd.bodyIDMin[0] && Id(p)<=gd.bodyIDMax[0]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( r/deltaR ) +1;
			if (i>cmd.sizeHistVcR) i=cmd.sizeHistVcR;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVcR[i] = gd.histVcR[i] + velt;
		}
		if ( (Id(p)>=gd.bodyIDMin[1] && Id(p)<=gd.bodyIDMax[1]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( r/deltaR ) +1;
			if (i>cmd.sizeHistVcR) i=cmd.sizeHistVcR;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVcR[i] = gd.histVcR[i] + velt;
		}
		if ( (Id(p)>=gd.bodyIDMin[2] && Id(p)<=gd.bodyIDMax[2]) 
				&& (Pos(p)[2]>=-hdz && Pos(p)[2]<=hdz) ) {
			xf=Pos(p)[0]; yf=Pos(p)[1];
			theta = angle(xi, yi, xf, yf);
			i = (int) ( r/deltaR ) +1;
			if (i>cmd.sizeHistVcR) i=cmd.sizeHistVcR;
			velt = Vel(p)[0]*(-rsin(theta)) + Vel(p)[1]*rcos(theta);
			gd.histVcR[i] = gd.histVcR[i] + velt;
		}
	}

	if (gd.countVcR == cmd.stepAvgVcR) {
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistVcR; j++)
			histSum= histSum+gd.histVcR[j];
		for (j=1; j<=cmd.sizeHistVcR; j++)
			gd.histVcR[j]= gd.histVcR[j]/histSum;
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

    outstr = stropen(fnametmp, "w!");

	for (n=1; n<=cmd.sizeHistRhoTheta; n++) {
		RBin = (n-0.5)*cmd.rangeR/cmd.sizeHistVcR;
		fprintf(outstr,"%8.3f %8.3f\n",RBin,gd.histVcR[n]);
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

