/*==============================================================================
	MODULE: analysis_md_io.c		[analysis_md]
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
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"


local void PrintRhoAxes(string, string);
local void PrintVelDist(string, string);
local void PrintRdf(string, string);

local void PrintNFrecAxes(string, string);
local void PrintNFrec(int *, int, char *, real);
local void ComputeFrequencies(int *, int *, int);


void startoutput(void)									// CHECK 2D --- OK!!!
{
//    printf("\n%s\n%s: %s\n\t %s\n", gd.headline0, gd.headline1, 
//									gd.headline2, gd.headline3);

    if (! strnull(cmd.options))                     
        printf("\n\toptions: %s\n", cmd.options);

	gd.count = 0;
    gd.histRhoX = AllocVecR(cmd.sizeHist);
    gd.histRhoY = AllocVecR(cmd.sizeHist);

    gd.histRhoX1 = AllocVecR(cmd.sizeHist);
    gd.histRhoY1 = AllocVecR(cmd.sizeHist);

    gd.histRhoX2 = AllocVecR(cmd.sizeHist);
    gd.histRhoY2 = AllocVecR(cmd.sizeHist);

#if (NDIM==3)
    gd.histRhoZ = AllocVecR(cmd.sizeHist);
    gd.histRhoZ1 = AllocVecR(cmd.sizeHist);
    gd.histRhoZ2 = AllocVecR(cmd.sizeHist);
#endif

//	gd.countRhoAxes = 0;
//    gd.histRhoX = AllocVecR(cmd.sizeHistRhoAxes);
//    gd.histRhoY = AllocVecR(cmd.sizeHistRhoAxes);
//    gd.histRhoX = AllocVecR(cmd.sizeHistRhoAxes);
//    gd.histRhoY = AllocVecR(cmd.sizeHistRhoAxes);

//    gd.histRhoX1 = AllocVecR(cmd.sizeHistRhoAxes);
//    gd.histRhoY1 = AllocVecR(cmd.sizeHistRhoAxes);

//    gd.histRhoX2 = AllocVecR(cmd.sizeHistRhoAxes);
//    gd.histRhoY2 = AllocVecR(cmd.sizeHistRhoAxes);

//#if (NDIM==3)
//    gd.histRhoZ = AllocVecR(cmd.sizeHistRhoAxes);
//    gd.histRhoZ1 = AllocVecR(cmd.sizeHistRhoAxes);
//    gd.histRhoZ2 = AllocVecR(cmd.sizeHistRhoAxes);
//#endif

//	gd.countNFrecAxes = 0;
    gd.histNFrecX = AllocVecI(cmd.sizeHist);
    gd.histNFrecY = AllocVecI(cmd.sizeHist);

    gd.histNFrecX1 = AllocVecI(cmd.sizeHist);
    gd.histNFrecY1 = AllocVecI(cmd.sizeHist);

    gd.histNFrecX2 = AllocVecI(cmd.sizeHist);
    gd.histNFrecY2 = AllocVecI(cmd.sizeHist);

	gd.histNFrecXD = AllocVecI(cmd.sizeHist);
	gd.histNFrecYD = AllocVecI(cmd.sizeHist);

#if (NDIM==3)
    gd.histNFrecZ = AllocVecI(cmd.sizeHist);
    gd.histNFrecZ1 = AllocVecI(cmd.sizeHist);
    gd.histNFrecZ2 = AllocVecI(cmd.sizeHist);
	gd.histNFrecZD = AllocVecI(cmd.sizeHist);
#endif

//	gd.countVel = 0;
    gd.histVel = AllocVecR(cmd.sizeHist);

	fprintf(gd.outlog,"stepAvg=%d\n",cmd.stepAvg);
	fprintf(gd.outlog,"sizeHist=%d\n",cmd.sizeHist);
	fprintf(gd.outlog,"rangeVal=%g\n",cmd.rangeVal);

//	fprintf(gd.outlog,"stepAvgRhoAxes=%d\n",cmd.stepAvgRhoAxes);
//	fprintf(gd.outlog,"sizeHistRhoAxes=%d\n",cmd.sizeHistRhoAxes);

//	fprintf(gd.outlog,"stepAvgVel=%d\n",cmd.stepAvgVel);
//	fprintf(gd.outlog,"sizeHistVel=%d\n",cmd.sizeHistVel);
//	fprintf(gd.outlog,"rangeVel=%g\n",cmd.rangeVel);

//	gd.countRdf = 0;
	gd.histRdf = AllocVecR(cmd.sizeHist);
	
//	fprintf(gd.outlog,"stepAvgRdf=%d\n",cmd.stepAvgRdf);
//	fprintf(gd.outlog,"sizeHistRdf=%d\n",cmd.sizeHistRdf);
//	fprintf(gd.outlog,"rangeRdf=%g\n",cmd.rangeRdf);
}

void EvalRhoAxes(string fname, string fnametmp)		// In TWODIM RhoZ is along y-axis...
{
    bodyptr p;
	real deltaX, deltaY, deltaZ, histSum;
	int i, l, j;
	double cpustart;
	real rangeX, rangeY, rangeZ;
	real vol;

	printf("\nEvalRhoAxes: Entrando ... ");
	cpustart = cputime();                       

	rangeX = gd.Box[0];
	rangeY = gd.Box[1];
#if (NDIM==3)
	rangeZ = gd.Box[2];
#endif

//	gd.countRhoAxes = gd.countRhoAxes + 1;
	gd.count = gd.count + 1;
	if (gd.count == 1) {
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoX[j] = 0.0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoY[j] = 0.0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoX1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoY1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoX2[j] = 0.0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoY2[j] = 0.0;
#if (NDIM==3)
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoZ[j] = 0.0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoZ1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histRhoZ2[j] = 0.0;
#endif
	}
	deltaX = rangeX/cmd.sizeHist;
	deltaY = rangeY/cmd.sizeHist;
#if (NDIM==3)
	deltaZ = rangeZ/cmd.sizeHist;
#endif

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		i = (int) ( (Pos(p)[0]+0.5*gd.Box[0] )/deltaX)+1;
		l = (int) ( (Pos(p)[1]+0.5*gd.Box[1] )/deltaY)+1;
#if (NDIM==3)
		j = (int) ( (Pos(p)[2]+0.5*gd.Box[2] )/deltaZ)+1;
#endif
		if (i>cmd.sizeHist) i=cmd.sizeHist;
			gd.histRhoX[i] = gd.histRhoX[i]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoX1[i] = gd.histRhoX1[i]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoX2[i] = gd.histRhoX2[i]+1.0;

		if (l>cmd.sizeHist) l=cmd.sizeHist;
			gd.histRhoY[l] = gd.histRhoY[l]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoY1[l] = gd.histRhoY1[l]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoY2[l] = gd.histRhoY2[l]+1.0;

#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		if (j>cmd.sizeHist) j=cmd.sizeHist;
			gd.histRhoZ[j] = gd.histRhoZ[j]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoZ1[j] = gd.histRhoZ1[j]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoZ2[j] = gd.histRhoZ2[j]+1.0;
#endif
	}

	if (gd.count == cmd.stepAvg) {
#if (NDIM==3)
		vol = gd.Box[1]*gd.Box[2]*deltaX;
#else
		vol = gd.Box[1]*deltaX;
#endif
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHist; j++)
			histSum= histSum+gd.histRhoX[j];
		for (j=1; j<=cmd.sizeHist; j++) {
//			gd.histRhoX[j]= gd.histRhoX[j]/histSum;
			gd.histRhoX[j]= gd.histRhoX[j]/((real)cmd.stepAvg*vol);
			gd.histRhoX1[j]= gd.histRhoX1[j]/((real)cmd.stepAvg*vol);
			gd.histRhoX2[j]= gd.histRhoX2[j]/((real)cmd.stepAvg*vol);
		}

#if (NDIM==3)
		vol = gd.Box[0]*gd.Box[2]*deltaY;
#else
		vol = gd.Box[0]*deltaY;
#endif
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHist; j++)
			histSum= histSum+gd.histRhoY[j];
		for (j=1; j<=cmd.sizeHist; j++) {
//			gd.histRhoY[j]= gd.histRhoY[j]/histSum;
			gd.histRhoY[j]= gd.histRhoY[j]/((real)cmd.stepAvg*vol);
			gd.histRhoY1[j]= gd.histRhoY1[j]/((real)cmd.stepAvg*vol);
			gd.histRhoY2[j]= gd.histRhoY2[j]/((real)cmd.stepAvg*vol);
		}

#if (NDIM==3)
		vol = gd.Box[0]*gd.Box[1]*deltaZ;
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHist; j++)
			histSum= histSum+gd.histRhoZ[j];
		for (j=1; j<=cmd.sizeHist; j++) {
//			gd.histRhoZ[j]= gd.histRhoZ[j]/(histSum*vol);
			gd.histRhoZ[j]= gd.histRhoZ[j]/((real)cmd.stepAvg*vol);
			gd.histRhoZ1[j]= gd.histRhoZ1[j]/((real)cmd.stepAvg*vol);
			gd.histRhoZ2[j]= gd.histRhoZ2[j]/((real)cmd.stepAvg*vol);
		}
#endif
//		PrintRhoAxes(stdout);
		PrintRhoAxes(fname, fnametmp);
		gd.RhoAxes_flag = 1;
		gd.count = 0;
	}

	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define RHOFMT	\
"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n"

local void PrintRhoAxes(string fname, string fnametmp)	// CHECK 2D --- OK!!!
{
	real xBin, yBin, zBin;
	int n;
    stream outstr_rhoz;
    char   buf[200];
	real rangeX, rangeY, rangeZ;

    outstr_rhoz = stropen(fnametmp, "w!");

	rangeX = gd.Box[0];
	rangeY = gd.Box[1];
#if (NDIM==3)
	rangeZ = gd.Box[2];
#endif

	for (n=1; n<=cmd.sizeHist; n++) {
		xBin = (n-0.5)*rangeX/cmd.sizeHist-0.5*gd.Box[0];
		yBin = (n-0.5)*rangeY/cmd.sizeHist-0.5*gd.Box[1];
#if (NDIM==3)
		zBin = (n-0.5)*rangeZ/cmd.sizeHist-0.5*gd.Box[2];
		fprintf(outstr_rhoz, RHOFMT,
			xBin,gd.histRhoX[n],gd.histRhoX1[n],gd.histRhoX2[n],
			yBin,gd.histRhoY[n],gd.histRhoY1[n],gd.histRhoY2[n],
			zBin,gd.histRhoZ[n],gd.histRhoZ1[n],gd.histRhoZ2[n]);
#else
		fprintf(outstr_rhoz,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
			xBin,gd.histRhoX[n],gd.histRhoX1[n],gd.histRhoX2[n],
			yBin,gd.histRhoY[n],gd.histRhoY1[n],gd.histRhoY2[n]);
#endif
	}
    fclose(outstr_rhoz);
	sprintf(buf,"mv %s %s",fnametmp,fname);
	printf("\nsystem: %s ...\n",buf);
	system(buf);
}

#undef RHOFMT


void EvalNFrecAxes(string fname, string fnametmp)
{
    bodyptr p;
	real deltaX, deltaY;
	int histSum;
	int i, l, j;
	double cpustart;
	real rangeX, rangeY;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
	real deltaZ, rangeZ;
#endif

	fprintf(stdout,"EvalNFrecAxes: Entrando ... ");
	cpustart = cputime();                       

	rangeX = gd.Box[0];
	rangeY = gd.Box[1];
#if (NDIM==3)
	rangeZ = gd.Box[2];
#endif

	gd.count = gd.count + 1;
	if (gd.count == 1) {
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecX[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecX1[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecX2[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecXD[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecY[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecY1[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecY2[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecYD[j] = 0;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecZ[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecZ1[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecZ2[j] = 0;
		for (j = 1; j <= cmd.sizeHist; j++) gd.histNFrecZD[j] = 0;
#endif
	}
	deltaX = rangeX/cmd.sizeHist;
	deltaY = rangeY/cmd.sizeHist;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
	deltaZ = rangeZ/cmd.sizeHist;
#endif

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		i = (int) ( (Pos(p)[0]+0.5*gd.Box[0] )/deltaX)+1;
		l = (int) ( (Pos(p)[1]+0.5*gd.Box[1] )/deltaY)+1;
#if (NDIM==3)
		j = (int) ( (Pos(p)[2]+0.5*gd.Box[2] )/deltaZ)+1;
#endif
		if (i>cmd.sizeHist) i=cmd.sizeHist;
			gd.histNFrecX[i] = gd.histNFrecX[i]+1;
		if ( Type(p) == BODY1 )
			gd.histNFrecX1[i] = gd.histNFrecX1[i]+1;
		if ( Type(p) == BODY2 )
			gd.histNFrecX2[i] = gd.histNFrecX2[i]+1;

		if (l>cmd.sizeHist) l=cmd.sizeHist;
			gd.histNFrecY[l] = gd.histNFrecY[l]+1;
		if ( Type(p) == BODY1 )
			gd.histNFrecY1[l] = gd.histNFrecY1[l]+1;
		if ( Type(p) == BODY2 )
			gd.histNFrecY2[l] = gd.histNFrecY2[l]+1;

#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		if (j>cmd.sizeHist) j=cmd.sizeHist;
			gd.histNFrecZ[j] = gd.histNFrecZ[j]+1;
		if ( Type(p) == BODY1 )
			gd.histNFrecZ1[j] = gd.histNFrecZ1[j]+1;
		if ( Type(p) == BODY2 )
			gd.histNFrecZ2[j] = gd.histNFrecZ2[j]+1;
#endif
	}

	if (gd.count == cmd.stepAvg) {
//		histSum = 0;
//		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++)
//			histSum= histSum+gd.histNFrecX[j];
//		for (j=1; j<=cmd.sizeHistNFrecAxes; j++)
//			gd.histNFrecX[j]= gd.nbody*gd.histNFrecX[j]/histSum;

//		histSum = 0;
//		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++)
//			histSum= histSum+gd.histNFrecY[j];
//		for (j=1; j<=cmd.sizeHistNFrecAxes; j++)
//			gd.histNFrecY[j]= gd.nbody*gd.histNFrecY[j]/histSum;

#if (NDIM==3)
//		histSum = 0;
//		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++)
//			histSum= histSum+gd.histNFrecZ[j];
//		for (j=1; j<=cmd.sizeHistNFrecAxes; j++)
//			gd.histNFrecZ[j]= gd.nbody*gd.histNFrecZ[j]/histSum;
#endif

//		PrintNFrecAxes(stdout);
		for (j = 1; j <= cmd.sizeHist; j++) {
			gd.histNFrecXD[j] = gd.histNFrecX1[j]-gd.histNFrecX2[j];
			gd.histNFrecYD[j] = gd.histNFrecY1[j]-gd.histNFrecY2[j];
#if (NDIM==3)
			gd.histNFrecZD[j] = gd.histNFrecZ1[j]-gd.histNFrecZ2[j];
#endif
		}
		PrintNFrecAxes(fname, fnametmp);
		gd.NFrecAxes_flag = 1;
		gd.count = 0;
	}

	printf("saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define savenfrecaxestmp "nfrecaxes.tmp"

local void PrintNFrecAxes(string fname, string fnametmp) // CHECK 2D --- OK!!!
{
	char fnamex[200], fnamey[200];
	char fnamex1[200], fnamey1[200];
	char fnamex2[200], fnamey2[200];
	char fnamexd[200], fnameyd[200];
#if (NDIM==3)
	char fnamez[200];
	char fnamez1[200];
	char fnamez2[200];
	char fnamezd[200];
#endif

	real xBin, yBin;
#if (NDIM==3)
	real zBin;
#endif
	int n;
    stream outstr_rhoz;
    char   buf[200];
	real rangeX, rangeY;
#if (NDIM==3)
	real rangeZ;
#endif
	real vol;

    outstr_rhoz = stropen(fnametmp, "w!");

	rangeX = gd.Box[0];
	rangeY = gd.Box[1];
#if (NDIM==3)
	rangeZ = gd.Box[2];
#endif

	for (n=1; n<=cmd.sizeHist; n++) {
		xBin = (n-0.5)*rangeX/cmd.sizeHist-0.5*gd.Box[0];
		yBin = (n-0.5)*rangeY/cmd.sizeHist-0.5*gd.Box[1];
#if (NDIM==3)
		zBin = (n-0.5)*rangeZ/cmd.sizeHist-0.5*gd.Box[2];
		fprintf(outstr_rhoz,"%8.3f %d %8.3f %d %8.3f %d\n",
			xBin,gd.histNFrecX[n], yBin,gd.histNFrecY[n], zBin,gd.histNFrecZ[n]);
#else
		fprintf(outstr_rhoz,"%8.3f %d %8.3f %d\n",
			xBin,gd.histNFrecX[n], yBin,gd.histNFrecY[n]);
#endif
	}
    fclose(outstr_rhoz);
	sprintf(buf,"mv %s %s",fnametmp,fname);
	printf("\nPrintNFrecAxes %d %d : system: %s ...\n",cmd.stepAvg, gd.count, buf);
	system(buf);


#if (NDIM==3)
	sprintf(fnamez, "%s-z", fname);
	sprintf(fnamez1, "%s-z1", fname);
	sprintf(fnamez2, "%s-z2", fname);
	sprintf(fnamezd, "%s-zd", fname);
#endif

	sprintf(fnamex, "%s-x", fname);
	sprintf(fnamey, "%s-y", fname);
	sprintf(fnamex1, "%s-x1", fname);
	sprintf(fnamey1, "%s-y1", fname);
	sprintf(fnamex2, "%s-x2", fname);
	sprintf(fnamey2, "%s-y2", fname);
	sprintf(fnamexd, "%s-xd", fname);
	sprintf(fnameyd, "%s-yd", fname);

//	PrintNFrec(&gd.histNFrecX[1], cmd.sizeHistNFrecAxes, savenfrecx);
//	PrintNFrec(&gd.histNFrecY[1], cmd.sizeHistNFrecAxes, savenfrecy);
#if (NDIM==3)
	vol = gd.Box[0]*gd.Box[1]*gd.Box[2]/cmd.sizeHist;
#else
	vol = gd.Box[0]*gd.Box[1]/cmd.sizeHist;
#endif
	PrintNFrec(&gd.histNFrecX[1], cmd.sizeHist, fnamex, vol);
	PrintNFrec(&gd.histNFrecX1[1], cmd.sizeHist, fnamex1, vol);
	PrintNFrec(&gd.histNFrecX2[1], cmd.sizeHist, fnamex2, vol);
	PrintNFrec(&gd.histNFrecXD[1], cmd.sizeHist, fnamexd, vol);
#if (NDIM==3)
	vol = gd.Box[0]*gd.Box[1]*gd.Box[2]/cmd.sizeHist;
#else
	vol = gd.Box[0]*gd.Box[1]/cmd.sizeHist;
#endif
	PrintNFrec(&gd.histNFrecY[1], cmd.sizeHist, fnamey, vol);
	PrintNFrec(&gd.histNFrecY1[1], cmd.sizeHist, fnamey1, vol);
	PrintNFrec(&gd.histNFrecY2[1], cmd.sizeHist, fnamey2, vol);
	PrintNFrec(&gd.histNFrecYD[1], cmd.sizeHist, fnameyd, vol);
#if (NDIM==3)
//	PrintNFrec(&gd.histNFrecZ[1], cmd.sizeHistNFrecAxes, savenfrecz);
	vol = gd.Box[0]*gd.Box[1]*gd.Box[2]/cmd.sizeHist;
	PrintNFrec(&gd.histNFrecZ[1], cmd.sizeHist, fnamez, vol);
	PrintNFrec(&gd.histNFrecZ1[1], cmd.sizeHist, fnamez1, vol);
	PrintNFrec(&gd.histNFrecZ2[1], cmd.sizeHist, fnamez2, vol);
	PrintNFrec(&gd.histNFrecZD[1], cmd.sizeHist, fnamezd, vol);
#endif
}

														// CHECK 2D --- OK!!!
local void PrintNFrec(int *rn, int nval, char *filename, real vol)
{
    stream outstr_nfrec;
    char   buf[200];
	int *n, *fn, *p, *s, i;

	n = (int *) allocate(nval*sizeof(int));
	s = (int *) allocate(nval*sizeof(int));
	fn = (int *) allocate(nval*sizeof(int));
	for (p=n; p<n+nval; p++)
		*p = rn[p-n];

	ComputeFrequencies(n, fn, nval);
	HeapSortInt(n, s, nval);

    outstr_nfrec = stropen(savenfrecaxestmp, "w!");
	for (i=0; i<nval; i++) {
		if (n[s[i]]!=-1)
			fprintf(outstr_nfrec,"%g %d\n",	
				((real)n[s[i]])/( ((real)cmd.stepAvg)*vol),fn[s[i]]);
	}
    fclose(outstr_nfrec);
	sprintf(buf,"mv %s %s",savenfrecaxestmp,filename);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);

}

//#undef savenfrecx
//#undef savenfrecy
//#if (NDIM==3)
//#undef savenfrecz
//#endif
#undef savenfrecaxestmp

local void ComputeFrequencies(int *n, int *fn, int nval) // CHECK 2D --- OK!!!
{
	int *p, *q;

	for (p=n; p<n+nval; p++) {
		if (*p!=-1)
			fn[p-n] = 1;
		else
			continue;
		for (q=p+1; q<n+nval; q++)
			if (*q != -1 && *p==*q) {
				*q=-1;
				fn[p-n] += 1;
			}
	}
}


void EvalVelDist(string savevel, string saveveltmp)		// CHECK 2D --- OK!!!
{
    bodyptr p;
	real deltaV, histSum, vv;
	int j, k;

	gd.count = gd.count + 1;
	if (gd.count == 1) {
		for (j = 1; j < cmd.sizeHist; j++) gd.histVel[j] = 0.0;
	}
	deltaV = cmd.rangeVal/cmd.sizeHist;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		vv =0.0;
		DO_COORD(k)
			vv += rsqr(Vel(p)[k]);
		j = (int) (rsqrt(vv)/deltaV)+1;
		if (j>cmd.sizeHist) j=cmd.sizeHist;
		gd.histVel[j] = gd.histVel[j]+1.0;
	}
	if (gd.count == cmd.stepAvg) {
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHist; j++)
			histSum= histSum+gd.histVel[j];
		for (j=1; j<=cmd.sizeHist; j++)
			gd.histVel[j]= gd.histVel[j]/histSum;
//		PrintVelDist(stdout);
		PrintVelDist(savevel, saveveltmp);
		gd.Vel_flag = 1;
		gd.count = 0.0;
	}
}

local void PrintVelDist(string savevel, string saveveltmp)// CHECK 2D --- OK!!!
{
	real vBin;
	int n;
    stream outstr_vel;
    char   buf[200];

    outstr_vel = stropen(saveveltmp, "w!");
	
	for (n=1; n<=cmd.sizeHist; n++) {
		vBin = (n-0.5)*cmd.rangeVal/cmd.sizeHist;
		fprintf(outstr_vel,"%8.3f %8.3f\n",vBin,gd.histVel[n]);
	}
    fclose(outstr_vel);
	sprintf(buf,"mv %s %s",saveveltmp,savevel);
	printf("\nsystem: %s",buf);
	system(buf);
}


void EvalRdf(string saverdf, string saverdftmp)			// CHECK 2D --- OK!!!
{
    bodyptr j1, j2;
	real deltaR, normFac, rr, rrRange, Vol;
	int k, n;
	vector dr;

	gd.count = gd.count + 1;
	if (gd.count == 1) {
		for (n = 1; n <= cmd.sizeHist; n++) gd.histRdf[n] = 0.0;
	}
	rrRange = rsqr(cmd.rangeVal);
	deltaR = cmd.rangeVal/cmd.sizeHist;
	DO_BODY(j1, bodytab, bodytab+gd.nbody-1) {
		DO_BODY(j2, j1+1, bodytab+gd.nbody) {
			DO_COORD(k) {
				dr[k] = Pos(j1)[k] - Pos(j2)[k];
				dr[k]=dr[k] - ( (real) (nint(dr[k]/gd.Box[k])) )*gd.Box[k];
			}
			rr=0.0;
			DO_COORD(k)
				rr += rsqr(dr[k]);
			if (rr < rrRange) {
				n = (int) (rsqrt(rr) / deltaR) + 1;
				gd.histRdf[n] = gd.histRdf[n] + 1.;
			}
		}
	}
		if (gd.count == cmd.stepAvg) {
		Vol = 1.0;
		DO_COORD(k)
			Vol = Vol*gd.Box[k];
		if (NDIM == 3)
			normFac = Vol /	(2.0 * PI * rpow(deltaR, 3.0) * gd.nbody * gd.nbody * gd.count);
		else if (NDIM == 2)
			normFac = Vol /	(PI * rpow(deltaR, 2.0) * gd.nbody * gd.nbody * gd.count);
		else error("\n\nWrong NDIM!\n\n");
		
		for (n = 1; n <= cmd.sizeHist; n++)
			if (NDIM == 3)
				gd.histRdf[n] = gd.histRdf[n] * normFac / rsqr(n-0.5);
			else if (NDIM == 2)
				gd.histRdf[n] = gd.histRdf[n] * normFac / ((int)n-0.5);
			else error("\n\nWrong NDIM!\n\n");

		PrintRdf(saverdf, saverdftmp);
		gd.Rdf_flag = 1;
		gd.count = 0;
	}
}


local void PrintRdf(string saverdf, string saverdftmp)
{
	real rBin;
	int n;
    stream outstr_rdf;
    char   buf[200];
	real ssq,fac,dr2,rri,rri3;
//	real UN,Pressure;
//	realptr UNIntegrand, PressureIntegrand;

    outstr_rdf = stropen(saverdftmp, "w!");
	
//	AllocMem (UNIntegrand, cmd.sizeHist, real);
//	AllocMem (PressureIntegrand, cmd.sizeHist, real);

	ssq = cmd.sigma11*cmd.sigma11;						// Una sola especie...

	for (n=1; n<=cmd.sizeHist; n++) {
		rBin = ((int)n-0.5)*cmd.rangeVal/cmd.sizeHist;
		dr2 = rBin*rBin;
		rri=ssq/dr2; rri3=rri*rri*rri;
//		UNIntegrand[n] = dr2*(rri3-1.0)*rri3*gd.histRdf[n];
//		PressureIntegrand[n] = dr2*rBin*rri3*(rri3-0.5)*rri*gd.histRdf[n];
	}
	fac = 2.0*PI*cmd.density*gdforce.fphi11;
//	UN = fac * Integrate(UNIntegrand, cmd.sizeHist);
	fac = (2.0/3.0)*PI*cmd.density*cmd.density*gdforce.fa11;
//	Pressure = cmd.density*cmd.temperature 
//				- fac * Integrate(PressureIntegrand, cmd.sizeHist);

//	fprintf(outstr_rdf,"# U/N, Pressure : %g %g\n",UN,Pressure);

	for (n=1; n<=cmd.sizeHist; n++) {
		rBin = ((int)n-0.5)*cmd.rangeVal/cmd.sizeHist;
		fprintf(outstr_rdf,"%8.4f %8.4f\n",rBin,gd.histRdf[n]);
	}
    fclose(outstr_rdf);
	sprintf(buf,"mv %s %s",saverdftmp,saverdf);
	printf("\nsystem: %s",buf);
	system(buf);

//printf("\n\nAntes ... \n\n");
//	free(UNIntegrand);				// Tiene una problema de reservacion de memoria!!
//	free(PressureIntegrand);		// Checar!!!! descomentando estas y las anteriores lineas...
//printf("\n\nDespues ... \n\n");
}

#ifdef THREEDIM
void PrintSnap_Slab(string savesnap, string savesnaptmp, int nbody, char *options,
                    char mode[2])
{
    stream outstr_snap;
    bodyptr p;
    char   buf[200];
	vector cmpos, tmpv;
	real xi, yi, zi;
	real tmass,r;
	
	fprintf(stdout,"\nPrintSnap_Slab: mode : %s\n",mode);
    outstr_snap = stropen(savesnaptmp, mode);
//    outstr_snap = stropen(savesnaptmp, "w!");
	
	fprintf(stdout,"In PrintSnap_Slab\n");
	
	DO_BODY(p, bodytab, bodytab+nbody) {
		if ( Pos(p)[0] >= cmd.xmin && Pos(p)[0] <= cmd.xmax &&
			Pos(p)[1] >= cmd.ymin && Pos(p)[1] <= cmd.ymax &&
			Pos(p)[2] >= cmd.zmin && Pos(p)[2] <= cmd.zmax ) {
			out_vector_mar(outstr_snap, Pos(p));             
			out_vector_mar(outstr_snap, Vel(p));
			out_int_mar(outstr_snap, Id(p));
			out_int_mar(outstr_snap, Type(p));
//			if (scanopt(options, "out-phi"))
//				out_real_mar(outstr_snap, Phi(p));
//			if (scanopt(options, "out-acc"))
//				out_vector_mar(outstr_snap, Acc(p));
			fprintf(outstr_snap,"\n");
		}
    }
    fclose(outstr_snap);
//	sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
	if (scanopt(mode, "a")) {
		sprintf(buf,"cp %s %s",savesnaptmp,savesnap);
		printf("\nsystem: %s",buf);
	} else {
		sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
		printf("\nsystem: %s",buf);
	}
	system(buf);
}
#endif

