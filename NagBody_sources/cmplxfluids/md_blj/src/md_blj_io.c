/*==============================================================================
	MODULE: md_blj_io.c			[md_blj]
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

	Major revisions: November 2008;
	Copyright: (c) 2005-2011 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#include "globaldefs.h"


local void EvalProps(void);
local void AccumProps(int);
local real EvalChemPot(void);
local void EvalLatticeCorr(void);
local void PrintSummary(FILE *);

local void EvalRhoAxes(void);
local void PrintRhoAxes(FILE *);
local void EvalNFrecAxes(void);
local void PrintNFrecAxes(FILE *);
local void PrintNFrec(int *, int, char *);
local void ComputeFrequencies(int *, int *, int);

local void EvalVelDist(void);
local void PrintVelDist(FILE *);
local void EvalRdf(void);
local void PrintRdf(FILE *);

local void EvalPressAxes(void);
local void PrintPressAxes(FILE *);

local void EvalDiffusion (void);
local void AccumDiffusion (void);
local void InitDiffusion (void);
local void ZeroDiffusion (void);
local void PrintDiffusion (FILE *);

local void EvalVelAcf(void);
local void AccumVelAcf(void);
local void InitVelAcf(void);
local void ZeroVelAcf(void);
local void PrintVelAcf(FILE *);
local real IntegrateVelAcf (real *, int);

local void EvaldPAcf(void);
local void AccumdPAcf(void);
local void InitdPAcf(void);
local void ZerodPAcf(void);
local void PrintdPAcf(FILE *);
local real IntegratedPAcf (real *, int);

local void EvalVacf(void);
local void AccumVacf(void);
local void InitVacf(void);
local void ZeroVacf(void);
local void PrintVacf(FILE *);
local real IntegrateVacf(real *, int);

local void EvalSpacetimeCorr (void);
local void AccumSpacetimeCorr (void);
local void InitSpacetimeCorr (void);
local void ZeroSpacetimeCorr (void);
local void PrintSpacetimeCorr (FILE *);


local void diagnostics(void);

local void Global_to_Header(void);

local void savestate(string);                 

local real mtot;                                
local real etot[3];                             
local matrix keten;                            
local matrix peten;                             
local vector cmpos;                             
local vector cmvel;
#if (NDIM==3)
local vector amvec;
#else
#if (NDIM==2)
local real amvec;
#endif
#endif

#define outfile_thermo	"thermo.dat"
#define outfile_transc	"transportc.dat"
#define savestatetmp	"savestate-tmp"

void StartOutput(void)									// CHECK 2D --- OK!!!
{
    stream outstr_thermo;
    stream outstr_transc;
	char buf[200];
	int k, nb;

	if (strnull(cmd.restorefile)) {
		sprintf(buf,"%s",outfile_thermo);
		if(!(outstr_thermo=fopen(buf,gd.mode)))
		{
			error("can't open file `%s`\n",buf);
		}
		fprintf(outstr_thermo,"%1s%4s%8s%8s%8s%8s%8s%8s%7s%8s%7s%7s%9s",
			"#","nstep","time","vSum","TotE","ssTotE",
			"PotE","ssPotE","kE","sskE","CV","p","ssp");
		fprintf(outstr_thermo,"%9s%7s%9s","p/rho*T","T","rho");
		fprintf(outstr_thermo,"%9s%8s%8s","Rcut11","Rcut12","Rcut22");
		fprintf(outstr_thermo,"%9s","lattCorr");
		if (cmd.computeChemPot)
			fprintf(outstr_thermo,"%8s%10s","ChemPot", "ssChemPot");
		fprintf(outstr_thermo,"\n");
		fprintf(outstr_thermo,
			"%1s%4s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%9s%8s%8s",
			"#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>",
			"<10>","<11>","<12>","<13>","<14>","<15>","<16>");
		fprintf(outstr_thermo,"%8s%8s%8s","<17>","<18>","<19>");
		if (cmd.computeChemPot)
			fprintf(outstr_thermo,"%9s%9s","<20>","<21>");
		fprintf(outstr_thermo,"\n");
		fclose(outstr_thermo);
	}

	if (cmd.computeVelAcf || cmd.computeBulkViscosity || cmd.computeTransport) {
	if (strnull(cmd.restorefile)) {
		sprintf(buf,"%s",outfile_transc);
		if(!(outstr_transc=fopen(buf,gd.mode)))
		{
			error("can't open file `%s`\n",buf);
		}
		fprintf(outstr_transc,"%1s%4s%8s%12s%16s%16s%20s%20s%20s%20s%20s%20s%21s\n",
			"#","nstep","time","Diffusion",
			"Bulk-viscosity","Bulk-viscosity-KK","Bulk-viscosity-KP","Bulk-viscosity-PP",
			"Shear-viscosity","Shear-viscosity-KK","Shear-viscosity-KP","Shear-viscosity-PP",
			"Thermal-conductivity");
		fclose(outstr_transc);
	}
	}

    fprintf(stdout,"  \t -- %s --\n", gd.model_comment);
	fprintf(stdout,"\n  \t running ... %s force calculation method ...\n",
		cmd.forcecalc_method);

    fprintf(stdout,"\n%12s%9s%8s%7s","temperature","density","nbody","dtime");
    fprintf(stdout,"%9s%10s%8s\n", "dtout", "dtoutinfo", "tstop");
    fprintf(stdout,"%8.2f     %8.4f%8d%8.5f", 
		cmd.temperature, cmd.density, gd.nbody, gd.dtime);
    fprintf(stdout,"%8.2f%8.2f %9.4f\n", gd.dtout, gd.dtoutinfo, cmd.tstop);

    if (! strnull(cmd.options))
        fprintf(stdout,"\n\toptions: %s\n", cmd.options);

	if (strnull(cmd.restorefile)) {

	if (cmd.computeRhoAxes) {
		gd.countRhoAxes = 0;
		gd.histRhoX = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoY = AllocVecR(cmd.sizeHistRhoAxes);

		gd.histRhoX1 = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoY1 = AllocVecR(cmd.sizeHistRhoAxes);

		gd.histRhoX2 = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoY2 = AllocVecR(cmd.sizeHistRhoAxes);
#if (NDIM==3)
		gd.histRhoZ = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoZ1 = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoZ2 = AllocVecR(cmd.sizeHistRhoAxes);
#endif
		fprintf(gd.outlog,"stepAvgRhoAxes=%d\n",cmd.stepAvgRhoAxes);
		fprintf(gd.outlog,"stepRhoAxes=%d\n",cmd.stepRhoAxes);
		fprintf(gd.outlog,"sizeHistRhoAxes=%d\n",cmd.sizeHistRhoAxes);
	}

	if (cmd.computeNFrecAxes) {
		gd.countNFrecAxes = 0;
		gd.histNFrecX = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecY = AllocVecI(cmd.sizeHistNFrecAxes);

		gd.histNFrecX1 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecY1 = AllocVecI(cmd.sizeHistNFrecAxes);
	
		gd.histNFrecX2 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecY2 = AllocVecI(cmd.sizeHistNFrecAxes);

		gd.histNFrecXD = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecYD = AllocVecI(cmd.sizeHistNFrecAxes);
#if (NDIM==3)
		gd.histNFrecZ = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecZ1 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecZ2 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecZD = AllocVecI(cmd.sizeHistNFrecAxes);
#endif

		fprintf(gd.outlog,"stepAvgNFrecAxes=%d\n",cmd.stepAvgNFrecAxes);
		fprintf(gd.outlog,"stepNFrecAxes=%d\n",cmd.stepNFrecAxes);
		fprintf(gd.outlog,"sizeHistNFrecAxes=%d\n",cmd.sizeHistNFrecAxes);
	}

    gd.stateVel = FALSE;
	if (cmd.computeVelDist) {
        gd.stateVel = TRUE;
		gd.countVel = 0;
		gd.histVel = AllocVecR(cmd.sizeHistVel);
		fprintf(gd.outlog,"stepAvgVel=%d\n",cmd.stepAvgVel);
		fprintf(gd.outlog,"stepVel=%d\n",cmd.stepVel);
		fprintf(gd.outlog,"sizeHistVel=%d\n",cmd.sizeHistVel);
		fprintf(gd.outlog,"rangeVel=%g\n",cmd.rangeVel);
	}

    gd.stateRdf = FALSE;
	if (cmd.computeRdf) {
        gd.stateRdf = TRUE;
		gd.countRdf = 0;
		gd.histRdf = AllocVecR(cmd.sizeHistRdf);
		gd.histRdf11 = AllocVecR(cmd.sizeHistRdf);
		gd.histRdf12 = AllocVecR(cmd.sizeHistRdf);
		gd.histRdf22 = AllocVecR(cmd.sizeHistRdf);
		fprintf(gd.outlog,"stepAvgRdf=%d\n",cmd.stepAvgRdf);
		fprintf(gd.outlog,"stepRdf=%d\n",cmd.stepRdf);
		fprintf(gd.outlog,"sizeHistRdf=%d\n",cmd.sizeHistRdf);
		fprintf(gd.outlog,"rangeRdf=%g\n",cmd.rangeRdf);
	}

	if (cmd.computePressAxes) {
		gd.countPressAxes = 0;
		gd.histPressX = AllocVecR(cmd.sizeHistPressAxes);
		gd.histPressY = AllocVecR(cmd.sizeHistPressAxes);
#if (NDIM==3)
		gd.histPressZ = AllocVecR(cmd.sizeHistPressAxes);
#endif

		fprintf(gd.outlog,"computePressAxes=%d\n",cmd.computePressAxes);
		fprintf(gd.outlog,"stepAvgPressAxes=%d\n",cmd.stepAvgPressAxes);
		fprintf(gd.outlog,"stepPressAxes=%d\n",cmd.stepPressAxes);
		fprintf(gd.outlog,"sizeHistPressAxes=%d\n",cmd.sizeHistPressAxes);
	}

	if (cmd.computeChemPot) {
		gd.countChemPot = 0;
		gd.histChemPotX = AllocVecR(cmd.sizeHistChemPot);
		gd.histChemPotY = AllocVecR(cmd.sizeHistChemPot);
#if (NDIM==3)
		gd.histChemPotZ = AllocVecR(cmd.sizeHistChemPot);
#endif

		fprintf(gd.outlog,"computeChemPot=%d\n",cmd.computeChemPot);
		fprintf(gd.outlog,"stepAvgChemPot=%d\n",cmd.stepAvgChemPot);
		fprintf(gd.outlog,"stepChemPot=%d\n",cmd.stepChemPot);
		fprintf(gd.outlog,"sizeHistChemPot=%d\n",cmd.sizeHistChemPot);
		fprintf(gd.outlog,"numTestBodies=%d\n",cmd.numTestBodies);
	}

	if (cmd.computeDiffusion) {
		AllocMem (gd.rrDiffuseAv, cmd.nValDiffuse, real);
		AllocMem (tBufD, cmd.nBuffDiffuse, TBufDiffusion);
		for (nb = 0; nb < cmd.nBuffDiffuse; nb ++) {
			AllocMem (tBufD[nb].orgR, gd.nbody, vector);
			AllocMem (tBufD[nb].rTrue, gd.nbody, vector);
			AllocMem (tBufD[nb].rrDiffuse, cmd.nValDiffuse, real);
		}
		if (strnull(cmd.restorefile))
			InitDiffusion();
	}

	if (cmd.computeVelAcf) {
		AllocMem (gd.avAcfVel, cmd.nValAcf, real);
		AllocMem (tBufVAcf, cmd.nBuffAcf, TBufVelAcf);
		for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
			AllocMem (tBufVAcf[nb].acfVel, cmd.nValAcf, real);
			AllocMem (tBufVAcf[nb].orgVel, gd.nbody, vector);
		}

		if (strnull(cmd.restorefile))
			InitVelAcf();
	}

	if (cmd.computeBulkViscosity) {
		AllocMem (gd.avAcfdP, cmd.nValAcf, real);
		AllocMem (tBufdPAcf, cmd.nBuffAcf, TBufdPressAcf);
		for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
			AllocMem (tBufdPAcf[nb].acfdP, cmd.nValAcf, real);
			AllocMem (tBufdPAcf[nb].orgdP, gd.nbody, real);
		}

		if (strnull(cmd.restorefile))
			InitdPAcf();
	}

	if (cmd.computeTransport) {
		AllocMem (gd.avAcfVel, cmd.nValAcf, real);
		AllocMem (gd.avAcfTherm, cmd.nValAcf, real);
// - Bulk Viscosity
		AllocMem (gd.avAcfBVisc_KK, cmd.nValAcf, real);
		AllocMem (gd.avAcfBVisc_KP, cmd.nValAcf, real);
		AllocMem (gd.avAcfBVisc_PP, cmd.nValAcf, real);
// - Shear Viscosity
		AllocMem (gd.avAcfVisc_KK, cmd.nValAcf, real);
		AllocMem (gd.avAcfVisc_KP, cmd.nValAcf, real);
		AllocMem (gd.avAcfVisc_PP, cmd.nValAcf, real);
//
		AllocMem (gd.avAcfVisc, cmd.nValAcf, real);
		AllocMem (tBuf, cmd.nBuffAcf, TBuf);
		for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
			AllocMem (tBuf[nb].acfVel, cmd.nValAcf, real);
			AllocMem (tBuf[nb].orgVel, gd.nbody, vector);
			AllocMem (tBuf[nb].acfTherm, cmd.nValAcf, real);
// - Bulk Viscosity
			AllocMem (tBuf[nb].acfBVisc_KK, cmd.nValAcf, real);
			AllocMem (tBuf[nb].acfBVisc_KP, cmd.nValAcf, real);
			AllocMem (tBuf[nb].acfBVisc_PP, cmd.nValAcf, real);
// - Shear Viscosity
			AllocMem (tBuf[nb].acfVisc_KK, cmd.nValAcf, real);
			AllocMem (tBuf[nb].acfVisc_KP, cmd.nValAcf, real);
			AllocMem (tBuf[nb].acfVisc_PP, cmd.nValAcf, real);

			AllocMem (tBuf[nb].acfVisc, cmd.nValAcf, real);
		}

		if (strnull(cmd.restorefile))
			InitVacf();
	}

	if (cmd.computeSTCorr) {
#if (NDIM==3)
		AllocMem (gd.valST, 24 * cmd.nFunCorr, real);	// 24 = 2 * 4 * 3
														// 2 - parte real e imaginaria
														// 4 = densidad + tres componentes
														//    del vector corriente
														// 3 = tres direcciones posibles
														//		del vector k : (k,0,0),
														//			(0,k,0), y (0,0,k).
#else
		AllocMem (gd.valST, 12 * cmd.nFunCorr, real);	// 12 = 2 * 3 * 2
#endif
		AllocMem2 (gd.avAcfST, 3 * cmd.nFunCorr, cmd.nValCorr, real);
		AllocMem (tBufC, cmd.nBuffCorr, TBufCorr);
		for (nb = 0; nb < cmd.nBuffCorr; nb ++) {
#if (NDIM==3)
			AllocMem (tBufC[nb].orgST, 24 * cmd.nFunCorr, real);
#else
			AllocMem (tBufC[nb].orgST, 12 * cmd.nFunCorr, real);
#endif
			AllocMem2 (tBufC[nb].acfST, 3 * cmd.nFunCorr, cmd.nValCorr, real);
		}

		InitSpacetimeCorr();
	}

	} else {
        if (!gd.stateVel) {
            if (cmd.computeVelDist) {
                gd.stateVel = TRUE;
                gd.countVel = 0;
                gd.histVel = AllocVecR(cmd.sizeHistVel);
                fprintf(gd.outlog,"stepAvgVel=%d\n",cmd.stepAvgVel);
                fprintf(gd.outlog,"stepVel=%d\n",cmd.stepVel);
                fprintf(gd.outlog,"sizeHistVel=%d\n",cmd.sizeHistVel);
                fprintf(gd.outlog,"rangeVel=%g\n",cmd.rangeVel);
            }
        }

        // Una vez fijos rangeRdf y sizeHistRdf no se podran modificar ...
        // Falta tener algo parecido para rangeRdf y sizeHistRdf
        // rangeRdf si funciona ... pero hay que dejar pasar la escritura
        // del archivo rdf.dat varias veces para que se estabilize el nuevo
        // valor de rangeRdf ...
        if (!gd.stateRdf) {
            if (cmd.computeRdf) {
                gd.stateRdf = TRUE;
                gd.countRdf = 0;
                gd.histRdf = AllocVecR(cmd.sizeHistRdf);
                gd.histRdf11 = AllocVecR(cmd.sizeHistRdf);
                gd.histRdf12 = AllocVecR(cmd.sizeHistRdf);
                gd.histRdf22 = AllocVecR(cmd.sizeHistRdf);
                fprintf(gd.outlog,"stepAvgRdf=%d\n",cmd.stepAvgRdf);
                fprintf(gd.outlog,"stepRdf=%d\n",cmd.stepRdf);
                fprintf(gd.outlog,"sizeHistRdf=%d\n",cmd.sizeHistRdf);
                fprintf(gd.outlog,"rangeRdf=%g\n",cmd.rangeRdf);
            }
        }
    }

	fprintf(gd.outlog,"\n\nstepEquil=%d\n",cmd.stepEquil);

// borrar cuando ya este seguro de que no sirven...
//	fprintf(gd.outlog,"stepSnapInit=%d\n",gd.stepSnapInit);
//	fprintf(gd.outlog,"stepSnap=%d\n\n",cmd.stepSnap);

	fflush(gd.outlog);

	if (! strnull(cmd.statefile)) {
		savestate(savestatetmp);
		sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
		printf("system: %s\n",buf);
		system(buf);
	}

	if (strnull(cmd.restorefile))
		AccumProps(0);

}

void output(void)										// CHECK 2D --- OK!!!
{
    real cmabs, amabs, teff;
    char   buf[200];
	real drmin, vmax, tmp, amax, rri, rri3;
	bodyptr p, q;
    double cpustart, cpudt;
	int flaginfo=0;
	char mode[4];

    cpustart = cputime();                       

	strcpy(mode,"w!");

	teff = gd.tnow + gd.dtime/8.0;

	EvalProps(); AccumProps(1);

    if (teff >= gd.toutinfo) {
		flaginfo=1;
		diagnostics();                              
		ABSV(cmabs, cmvel);
#if (NDIM==3)
		ABSV(amabs, amvec);
#else
#if (NDIM==2)
		amabs = rabs(amvec);
#endif
#endif
		printf("\n\n%8s %8s %8s %10s %10s %7s %10s\n", "nstep",
           "time", "|Vcom|", "|Jtot|", "nbbtot", "CPUfc", "CPUTotFc");
		gd.cputotforce += gdforce.cpuforce;
		printf("%8d %8.3f %8.5f  %9.5f %9d %8.3f %9.3f\n", gd.nstep,
           gd.tnow,cmabs,amabs,gdforce.nbbcalc,gdforce.cpuforce,gd.cputotforce);
		if (! scanopt(cmd.forcecalc_method, "direct")  &&
			! scanopt(cmd.forcecalc_method, "direct2") &&
			! scanopt(cmd.forcecalc_method, "cells")	) {
			printf("\n\t%8s%8s%8s%8s%10s%10s%8s\n",
				"rsize", "tdepth", "ftree",
				"actmax", "nbbtot", "nbctot", "CPUtree");
			printf("\t%8.1f%8d%8.3f%8d%10d%10d%8.3f\n",
			gdtree.rsize, gdtree.tdepth, 
			(gd.nbody + gdtree.ncell - 1) / ((real) gdtree.ncell),
			gdforce.actmax, gdforce.nbbcalc, gdforce.nbccalc, gdtree.cputree);
		}
		printf("\n%9s %6s %9s %9s %10s %8s %8s\n", "Averages:",
           "Etot/N", "kinE/N", "PotE/N", "pressure", "vSum", "vvSum");
		printf("       %9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n\n",
           gd.totEnergy,gd.kinEnergy,gd.potEnergy,gd.pressure,gd.vSum,gd.vvSum);

#if (NDIM==3)
printf("\n\nCM pos, vel : %g %g %g %g %g %g\n\n",
	cmpos[0],cmpos[1],cmpos[2],cmvel[0],cmvel[1],cmvel[2]);
#else
printf("\n\nCM pos, vel : %g %g %g %g\n\n",cmpos[0],cmpos[1],cmvel[0],cmvel[1]);
#endif

		gd.toutinfo += gd.dtoutinfo;
	}

	if (cmd.computeRhoAxes) {
		if (gd.nstep>=cmd.stepEquil 
				&& (gd.nstep-cmd.stepEquil)%cmd.stepRhoAxes==0)
			EvalRhoAxes();
	}

	if (cmd.computeNFrecAxes) {
		if (gd.nstep>=cmd.stepEquil 
				&& (gd.nstep-cmd.stepEquil)%cmd.stepNFrecAxes==0)
			EvalNFrecAxes();
	}

	if (cmd.computeVelDist) {
		if (gd.nstep>=cmd.stepEquil && (gd.nstep-cmd.stepEquil)%cmd.stepVel==0)
			EvalVelDist();
	}

	if (cmd.computeRdf) {
		if (gd.nstep>=cmd.stepEquil && (gd.nstep-cmd.stepEquil)%cmd.stepRdf==0)
			EvalRdf();
	}

//	if (cmd.printSnap) { 	// borrar cuando ya este seguro de que no sirven...
//		if (gd.nstep>=gd.stepSnapInit 
//				&& (gd.nstep-gd.stepSnapInit)%cmd.stepSnap==0)
//			PrintSnap("snap.dat", "snap.tmp", gd.nbody, cmd.options, mode);
//	}

	if (gd.stepAvgFlag) {
		if ((gd.nstepNew)%gd.stepAvgOld == 0) {		// Restore behaviour ...
			AccumProps(2);
			PrintSummary(stdout); AccumProps(0);
		}
	} else {
		if (!strnull(cmd.restorefile) && gd.stepAvgStatus) {
			if ((gd.nstepNew-gd.nstepOld)%cmd.stepAvg == 0) {	// Restore behaviour ...
				AccumProps(2);
				PrintSummary(stdout); AccumProps(0);
			}
		} else {
			if (gd.nstepNew%cmd.stepAvg == 0) {				// Normal behaviour ...
				AccumProps(2);
				PrintSummary(stdout); AccumProps(0);
			}
		}
	}

	if (cmd.computePressAxes) {
		if (gd.nstep>=cmd.stepEquil
				&&(gd.nstep - cmd.stepEquil)%cmd.stepPressAxes==0)
			EvalPressAxes();
	}

	if (cmd.computeDiffusion) {
		if (gd.nstep>=cmd.stepEquil&&(gd.nstep - cmd.stepEquil)%cmd.stepDiffuse==0)
			EvalDiffusion();
	}

	if (cmd.computeVelAcf) {
		if (gd.nstep>=cmd.stepEquil&&(gd.nstep - cmd.stepEquil)%cmd.stepAcf==0) 
			EvalVelAcf();
	}

	if (cmd.computeBulkViscosity) {
		if (gd.nstep>=cmd.stepEquil&&(gd.nstep - cmd.stepEquil)%cmd.stepAcf==0) 
			EvaldPAcf();
	}

	if (cmd.computeTransport) {
		if (gd.nstep>=cmd.stepEquil&&(gd.nstep - cmd.stepEquil)%cmd.stepAcf==0) 
			EvalVacf();
	}

	if (cmd.computeSTCorr) {
		if (gd.nstep>=cmd.stepEquil&&(gd.nstep - cmd.stepEquil)%cmd.stepCorr==0) 
			EvalSpacetimeCorr();
	}

    if (! strnull(cmd.snapoutfile) && teff >= gd.tout) {
		Global_to_Header();
		outputdata(cmd.snapoutfile, cmd.snapoutfilefmt, gd.headerfmt, gd.nstep, 
					gd.nbody, gd.tnow, &hdr, cmd.options);
		gd.tout += gd.dtout;
    }

	if (gd.nstep%cmd.stepState == 0) {				// stepSnap = dtout/dtime
        if (! strnull(cmd.statefile)) {
			cpudt = cputime()-gd.cpuinit;
			gd.cputotal += cpudt;
			Global_to_Header();
            savestate(savestatetmp);
			gd.cputotal -= cpudt;
            sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
            printf("system: %s\n",buf);
            system(buf);
        }
	}

	gd.cputotout += cputime()-cpustart;
	if (flaginfo)
		printf("Individual and Accum. output cpu time : %g %g\n",
			cputime()-cpustart, gd.cputotout);
}
											// CHECK 2D --- OK!!!
local void diagnostics(void)				// se mantiene para chequeo cruzado
{
    register bodyptr p;
    real velsq;
	vector tmpv;
#if (NDIM==3)
    vector tmpav;
#else
#if (NDIM==2)
    real tmpav;
#endif
#endif
    matrix tmpt;

    mtot = 0.0;                                 
    etot[1] = etot[2] = 0.0;                    
    CLRM(keten);                                
    CLRM(peten);
#if (NDIM==3)
		CLRV(amvec);                                
#else
#if (NDIM==2)
		amvec=0.0;
#endif
#endif
    CLRV(cmpos);                                
    CLRV(cmvel);                                
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        mtot += Mass(p);                        
        DOTVP(velsq, Vel(p), Vel(p));           
        etot[1] += 0.5 * Mass(p) * velsq;       
        etot[2] += 0.5 * Mass(p) * Phi(p);      
        MULVS(tmpv, Vel(p), 0.5 * Mass(p));     
        OUTVP(tmpt, tmpv, Vel(p));
        ADDM(keten, keten, tmpt);
        MULVS(tmpv, Pos(p), Mass(p));           
        OUTVP(tmpt, tmpv, Acc(p));
        ADDM(peten, peten, tmpt);
#if (NDIM==3)
			CROSSVP(tmpav, Vel(p), Pos(p));          
			MULVS(tmpav, tmpav, Mass(p));
			ADDV(amvec, amvec, tmpav);
#else
#if (NDIM==2)
			CROSSVP(tmpav, Vel(p), Pos(p));          
			tmpav=tmpav*Mass(p);
			amvec=amvec+tmpav;
#endif
#endif
        MULVS(tmpv, Pos(p), Mass(p));           
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           
        ADDV(cmvel, cmvel, tmpv);
    }
    etot[0] = etot[1] + etot[2];                
    DIVVS(cmpos, cmpos, mtot);                  
    DIVVS(cmvel, cmvel, mtot);
}

local void EvalProps(void)								// CHECK 2D --- OK!!!
{
    bodyptr p;
	real v, vv;
	int k;
	
	gd.vSum = gd.vvSum = 0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		vv = 0.0;
		DO_COORD(k) {							// No se sincronizan vels
			v = Vel(p)[k];						// - 0.5*Acc(p)[k]*cmd.dtime;
			gd.vSum += v; vv += rsqr(v);
		}
		gd.vvSum += Mass(p)*vv;
	}
	gd.kinEnergy = 0.5*gd.vvSum/((real) (gd.nbody));
	gd.potEnergy = gdforce.uSum/((real) (gd.nbody));
	gd.totEnergy = gd.kinEnergy + gd.potEnergy;
	gd.pressure = cmd.density*(gd.vvSum+gdforce.virSum)/((real)(gd.nbody*NDIM));
//	gd.PV_K = cmd.density*(gd.vvSum)/((real)(gd.nbody*NDIM));
//	gd.PV_P = cmd.density*(gdforce.virSum)/((real)(gd.nbody*NDIM));
	gd.PV_K = (gd.vvSum);
	gd.PV_P = (gdforce.virSum);
	EvalLatticeCorr();
	if (cmd.computeChemPot) gd.ChemPot = EvalChemPot();
}

local void AccumProps(int icode)						// CHECK 2D --- OK!!!
{
	real stepAvg;

	if (icode == 0) {
		gd.sTotEnergy = gd.ssTotEnergy = 0.0;
		gd.sKinEnergy = gd.ssKinEnergy = 0.0;
		gd.sPotEnergy = gd.ssPotEnergy = 0.0;
		gd.sPressure = gd.ssPressure = 0.0;

		gd.sPV_K = gd.ssPV_K = 0.0;
		gd.sPV_P = gd.ssPV_P = 0.0;

		if (cmd.computeChemPot) {
			gd.sChemPot = gd.ssChemPot = 0.0;
		}
	} if (icode == 1) {
		gd.sTotEnergy += gd.totEnergy;
		gd.ssTotEnergy += rsqr(gd.totEnergy);
		gd.sKinEnergy += gd.kinEnergy;
		gd.ssKinEnergy += rsqr(gd.kinEnergy);
		gd.sPotEnergy += gd.potEnergy;
		gd.ssPotEnergy += rsqr(gd.potEnergy);
		gd.sPressure += gd.pressure;
		gd.ssPressure += rsqr(gd.pressure);

		gd.sPV_K += gd.PV_K;
		gd.ssPV_K += rsqr(gd.PV_K);
		gd.sPV_P += gd.PV_P;
		gd.ssPV_P += rsqr(gd.PV_P);

		if (cmd.computeChemPot) {
			gd.sChemPot += gd.ChemPot;
			gd.ssChemPot += rsqr(gd.ChemPot);
		}
	} if (icode == 2) {
		if (gd.tnow==0) 
			stepAvg = 1.0;
		else {
			if (gd.stepAvgFlag) {				// In order to restore work
				gd.stepAvgFlag = FALSE;
				stepAvg = gd.stepAvgOld;
				gd.nstepOld = gd.nstepNew;
			} else
				stepAvg = cmd.stepAvg;
		}
		gd.sTotEnergy = gd.sTotEnergy/stepAvg;
		gd.ssTotEnergy = rsqrt((gd.ssTotEnergy/stepAvg>rsqr(gd.sTotEnergy) ? 
			gd.ssTotEnergy/stepAvg - rsqr(gd.sTotEnergy) : 0) );
		gd.sKinEnergy = gd.sKinEnergy/stepAvg;
		gd.ssKinEnergy = rsqrt((gd.ssKinEnergy/stepAvg>rsqr(gd.sKinEnergy) ?
			gd.ssKinEnergy/stepAvg - rsqr(gd.sKinEnergy) : 0) );

		gd.sPotEnergy = gd.sPotEnergy/stepAvg;
		gd.ssPotEnergy = rsqrt((gd.ssPotEnergy/stepAvg>rsqr(gd.sPotEnergy) ?
			gd.ssPotEnergy/stepAvg - rsqr(gd.sPotEnergy) : 0) );

		gd.sPressure = gd.sPressure/stepAvg;
		gd.ssPressure = rsqrt((gd.ssPressure/stepAvg>rsqr(gd.sPressure) ?
			gd.ssPressure/stepAvg - rsqr(gd.sPressure) : 0) );

		gd.sPV_K = gd.sPV_K/stepAvg;
		gd.ssPV_K = rsqrt((gd.ssPV_K/stepAvg>rsqr(gd.sPV_K) ?
			gd.ssPV_K/stepAvg - rsqr(gd.sPV_K) : 0) );
		gd.sPV_P = gd.sPV_P/stepAvg;
		gd.ssPV_P = rsqrt((gd.ssPV_P/stepAvg>rsqr(gd.sPV_P) ?
			gd.ssPV_P/stepAvg - rsqr(gd.sPV_P) : 0) );

		if (cmd.computeBulkViscosity) gd.sPressureSave=0.5*(gd.sPressureSave+gd.sPressure);

		if (cmd.computeTransport) {
			gd.PV_KAvgSave = gd.sPV_K;
			gd.PV_PAvgSave = gd.sPV_P;
		}

		if (cmd.computeChemPot) {
			gd.sChemPot = gd.sChemPot/stepAvg;
			gd.ssChemPot = rsqrt((gd.ssChemPot/stepAvg>rsqr(gd.sChemPot) ?
			gd.ssChemPot/stepAvg - rsqr(gd.sChemPot) : 0) );
		}
	}
}

local real EvalChemPot(void)
{
	real chempot;
	real ustmp, virtmp, BSum, up, virp, xran, yran;
	short flag;
#if (NDIM==3)
	real zran;
#endif
	bodyptr p;
	double cpustart;
	int i, iter, maxiter=10000;

	cpustart = cputime();                       

	BSum = 0.;
	ustmp=gdforce.uSum;
	virtmp=gdforce.virSum;

	p = bodytab+gd.nbody;
	i = 1, iter=1;
	while (i <= cmd.numTestBodies) {
		Mass(p) = (xrandom(0.,1.) < 0.5 ? gd.mass1 : gd.mass2);
		Type(p) = (xrandom(0.,1.) < 0.5 ? BODY1 : BODY2);
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		xran = xrandom(-0.5*gdforce.Box[0], 0.5*gdforce.Box[0]);
		yran = xrandom(-0.5*gdforce.Box[1], 0.5*gdforce.Box[1]);
		zran = xrandom(-0.5*gdforce.Box[2], 0.5*gdforce.Box[2]);
		VSet(Pos(p), xran, yran, zran);
#else
		xran = xrandom(-0.5*gdforce.Box[0], 0.5*gdforce.Box[0]);
		yran = xrandom(-0.5*gdforce.Box[1], 0.5*gdforce.Box[1]);
		VSet(Pos(p), xran, yran);
#endif
		Update(p) = TRUE;
		Rcut(p) = (Type(p)==BODY1 ? gdforce.Rcut11Max : gdforce.Rcut22Max);

		ind_ljforcecalc_direct(bodytab, gd.nbody, &gdforce, 
								p, &up, &virp, &flag);
		if (flag == 1) {
			++i;
			BSum += rexp(-up/cmd.temperature);
		}
		++iter;
		if (iter>maxiter) error("\n\nEvalChemPot: max iter reached\n");
	}
	BSum /= cmd.numTestBodies;
	chempot = -cmd.temperature * rlog(BSum)/LNE;

	gdforce.uSum = ustmp;
	gdforce.virSum = virtmp;

	Mass(p) = Mass(p-1);
	Id(p) = p-bodytab+1;
	Type(p) = TESTBODYMU;
	CLRV(Pos(p));
	CLRV(Vel(p));
	CLRV(Acc(p));
	Update(p) = FALSE;
	Rcut(p) = gd.RcutAllMax;

	gd.cpuchempot = cputime() - cpustart;
	return (chempot);
}

local void EvalLatticeCorr(void)
{
  vector kVec;
  real si, sr, t;
  int n;
  bodyptr p;

//  kVec[0] = 2. * PI * gd.numUcell[0] / gdforce.Box[0];
//  kVec[1] = - kVec[0];
//#if (NDIM==3)										// CORRECCION 2-3D CHECAR
//  kVec[2] = kVec[0];
//#endif

  kVec[0] = 2. * PI * cmd.lattCorr_kx * gd.numUcell[0] / gdforce.Box[0];
  kVec[1] = 2. * PI * cmd.lattCorr_ky * gd.numUcell[1] / gdforce.Box[1];
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
  kVec[2] = 2. * PI * cmd.lattCorr_kz * gd.numUcell[2] / gdforce.Box[2];
#endif

  sr = 0.;
  si = 0.;
  DO_BODY(p, bodytab, bodytab+gd.nbody) {
    t = VDot(kVec, Pos(p));
    sr += rcos(t);
    si += rsin(t);
  }
  gd.latticeCorr = rsqrt(rsqr(sr) + rsqr(si)) / gd.nbody;
}

#define THERMOFMT	\
"%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f \
%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f"

local void PrintSummary(FILE *fp)
{
    stream outstr_thermo;
	real CV;

	CV=1.5/(1.-(2.*gd.nbody*rsqr(gd.ssPotEnergy))/(3.*rsqr(cmd.temperature)));
    outstr_thermo = stropen(outfile_thermo, "a");
	fprintf(outstr_thermo, THERMOFMT,
		gd.nstep, gd.tnow, gd.vSum, gd.sTotEnergy, gd.ssTotEnergy,
		gd.sPotEnergy, gd.ssPotEnergy,
		gd.sKinEnergy, gd.ssKinEnergy, CV, 
		gd.sPressure, gd.ssPressure,
		gd.sPressure/(cmd.density * cmd.temperature), 
		cmd.temperature, cmd.density, cmd.Rcut11, cmd.Rcut12, cmd.Rcut22,
		gd.latticeCorr);
	if (cmd.computeChemPot) fprintf(outstr_thermo," %7.4f %7.4f",
		gd.sChemPot, gd.ssChemPot);
	fprintf(outstr_thermo,"\n");
    fclose(outstr_thermo);
}

#undef THERMOFMT
#undef outfile_thermo

local void EvalRhoAxes(void)
{
    bodyptr p;
	real deltaX, deltaY, histSum;
	int i, l, j;
	double cpustart;
	real rangeX, rangeY;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
	real deltaZ, rangeZ;
#endif
	real vol;

	printf("EvalRhoAxes: Entrando ... ");
	cpustart = cputime();                       

	rangeX = gdforce.Box[0];
	rangeY = gdforce.Box[1];
#if (NDIM==3)
	rangeZ = gdforce.Box[2];
#endif

	gd.countRhoAxes = gd.countRhoAxes + 1;
	if (gd.countRhoAxes == 1) {
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoX[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoY[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoX1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoY1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoX2[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoY2[j] = 0.0;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoZ[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoZ1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoZ2[j] = 0.0;
#endif
	}
	deltaX = rangeX/cmd.sizeHistRhoAxes;
	deltaY = rangeY/cmd.sizeHistRhoAxes;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
	deltaZ = rangeZ/cmd.sizeHistRhoAxes;
#endif

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		i = (int) ( (Pos(p)[0]+0.5*gdforce.Box[0] )/deltaX)+1;
		l = (int) ( (Pos(p)[1]+0.5*gdforce.Box[1] )/deltaY)+1;
#if (NDIM==3)
		j = (int) ( (Pos(p)[2]+0.5*gdforce.Box[2] )/deltaZ)+1;
#endif
		if (i>cmd.sizeHistRhoAxes) i=cmd.sizeHistRhoAxes;
			gd.histRhoX[i] = gd.histRhoX[i]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoX1[i] = gd.histRhoX1[i]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoX2[i] = gd.histRhoX2[i]+1.0;

		if (l>cmd.sizeHistRhoAxes) l=cmd.sizeHistRhoAxes;
			gd.histRhoY[l] = gd.histRhoY[l]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoY1[l] = gd.histRhoY1[l]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoY2[l] = gd.histRhoY2[l]+1.0;

#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		if (j>cmd.sizeHistRhoAxes) j=cmd.sizeHistRhoAxes;
			gd.histRhoZ[j] = gd.histRhoZ[j]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoZ1[j] = gd.histRhoZ1[j]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoZ2[j] = gd.histRhoZ2[j]+1.0;
#endif
	}

	if (gd.countRhoAxes == cmd.stepAvgRhoAxes) {
		vol = gdforce.Box[0]*gdforce.Box[1]*deltaX;
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoAxes; j++)
			histSum= histSum+gd.histRhoX[j];
		for (j=1; j<=cmd.sizeHistRhoAxes; j++) {
			gd.histRhoX[j]= gd.histRhoX[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoX1[j]= gd.histRhoX1[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoX2[j]= gd.histRhoX2[j]/((real)cmd.stepAvgRhoAxes*vol);
		}

		vol = gdforce.Box[0]*gdforce.Box[1]*deltaY;
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoAxes; j++)
			histSum= histSum+gd.histRhoY[j];
		for (j=1; j<=cmd.sizeHistRhoAxes; j++) {
			gd.histRhoY[j]= gd.histRhoY[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoY1[j]= gd.histRhoY1[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoY2[j]= gd.histRhoY2[j]/((real)cmd.stepAvgRhoAxes*vol);
		}

#if (NDIM==3)
		vol = gdforce.Box[0]*gdforce.Box[1]*deltaZ;
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoAxes; j++)
			histSum= histSum+gd.histRhoZ[j];
		for (j=1; j<=cmd.sizeHistRhoAxes; j++) {
			gd.histRhoZ[j]= gd.histRhoZ[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoZ1[j]= gd.histRhoZ1[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoZ2[j]= gd.histRhoZ2[j]/((real)cmd.stepAvgRhoAxes*vol);
		}
#endif
		PrintRhoAxes(stdout);
		gd.countRhoAxes = 0.0;
	}

	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define saverhoaxestmp "rhoaxes.tmp"
#define saverhoaxes "rhoaxes.dat"
#define RHOFMT	\
"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n"

local void PrintRhoAxes(FILE *fp)
{
	real xBin, yBin;
	int n;
    stream outstr_rhoz;
    char   buf[200];
	real rangeX, rangeY;
#if (NDIM==3)
	real zBin, rangeZ;
#endif

    outstr_rhoz = stropen(saverhoaxestmp, "w!");

	rangeX = gdforce.Box[0];
	rangeY = gdforce.Box[1];
#if (NDIM==3)
	rangeZ = gdforce.Box[2];
#endif

	fprintf(outstr_rhoz,"# %8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s",
		"x", "RhoXTot", "RhoX1", "RhoX2",
		"y", "RhoYTot", "RhoY1", "RhoY2",
		"z", "RhoZTot", "RhoZ1", "RhoZ2");

	for (n=1; n<=cmd.sizeHistRhoAxes; n++) {
		xBin = (n-0.5)*rangeX/cmd.sizeHistRhoAxes-0.5*gdforce.Box[0];
		yBin = (n-0.5)*rangeY/cmd.sizeHistRhoAxes-0.5*gdforce.Box[1];
#if (NDIM==3)
		zBin = (n-0.5)*rangeZ/cmd.sizeHistRhoAxes-0.5*gdforce.Box[2];
		fprintf(outstr_rhoz, RHOFMT,
//			"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
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
	sprintf(buf,"mv %s %s",saverhoaxestmp,saverhoaxes);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);
}

#undef saverhoaxes
#undef saverhoaxestmp
#undef RHOFMT

local void EvalNFrecAxes(void)
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

	rangeX = gdforce.Box[0];
	rangeY = gdforce.Box[1];
#if (NDIM==3)
	rangeZ = gdforce.Box[2];
#endif

	gd.countNFrecAxes = gd.countNFrecAxes + 1;
	if (gd.countNFrecAxes == 1) {
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecX[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecX1[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecX2[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecXD[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecY[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecY1[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecY2[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecYD[j] = 0;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecZ[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecZ1[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecZ2[j] = 0;
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) gd.histNFrecZD[j] = 0;
#endif
	}
	deltaX = rangeX/cmd.sizeHistNFrecAxes;
	deltaY = rangeY/cmd.sizeHistNFrecAxes;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
	deltaZ = rangeZ/cmd.sizeHistNFrecAxes;
#endif

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		i = (int) ( (Pos(p)[0]+0.5*gdforce.Box[0] )/deltaX)+1;
		l = (int) ( (Pos(p)[1]+0.5*gdforce.Box[1] )/deltaY)+1;
#if (NDIM==3)
		j = (int) ( (Pos(p)[2]+0.5*gdforce.Box[2] )/deltaZ)+1;
#endif
		if (i>cmd.sizeHistNFrecAxes) i=cmd.sizeHistNFrecAxes;
			gd.histNFrecX[i] = gd.histNFrecX[i]+1;
		if ( Type(p) == BODY1 )
			gd.histNFrecX1[i] = gd.histNFrecX1[i]+1;
		if ( Type(p) == BODY2 )
			gd.histNFrecX2[i] = gd.histNFrecX2[i]+1;

		if (l>cmd.sizeHistNFrecAxes) l=cmd.sizeHistNFrecAxes;
			gd.histNFrecY[l] = gd.histNFrecY[l]+1;
		if ( Type(p) == BODY1 )
			gd.histNFrecY1[l] = gd.histNFrecY1[l]+1;
		if ( Type(p) == BODY2 )
			gd.histNFrecY2[l] = gd.histNFrecY2[l]+1;

#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		if (j>cmd.sizeHistNFrecAxes) j=cmd.sizeHistNFrecAxes;
			gd.histNFrecZ[j] = gd.histNFrecZ[j]+1;
		if ( Type(p) == BODY1 )
			gd.histNFrecZ1[j] = gd.histNFrecZ1[j]+1;
		if ( Type(p) == BODY2 )
			gd.histNFrecZ2[j] = gd.histNFrecZ2[j]+1;
#endif
	}

	if (gd.countNFrecAxes == cmd.stepAvgNFrecAxes) {
		for (j = 1; j <= cmd.sizeHistNFrecAxes; j++) {
			gd.histNFrecXD[j] = gd.histNFrecX1[j]-gd.histNFrecX2[j];
			gd.histNFrecYD[j] = gd.histNFrecY1[j]-gd.histNFrecY2[j];
#if (NDIM==3)
			gd.histNFrecZD[j] = gd.histNFrecZ1[j]-gd.histNFrecZ2[j];
#endif
		}
		PrintNFrecAxes(stdout);
		gd.countNFrecAxes = 0;
	}

	printf("saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define savenfrecaxestmp "nfrecaxes.tmp"
#define savenfrecx "nfrecx.dat"
#define savenfrecx1 "nfrecx1.dat"
#define savenfrecx2 "nfrecx2.dat"
#define savenfrecxd "nfrecxd.dat"
#define savenfrecy "nfrecy.dat"
#define savenfrecy1 "nfrecy1.dat"
#define savenfrecy2 "nfrecy2.dat"
#define savenfrecyd "nfrecyd.dat"
#if (NDIM==3)
#define savenfrecz "nfrecz.dat"
#define savenfrecz1 "nfrecz1.dat"
#define savenfrecz2 "nfrecz2.dat"
#define savenfreczd "nfreczd.dat"
#endif

local void PrintNFrecAxes(FILE *fp)
{
	PrintNFrec(&gd.histNFrecX[1], cmd.sizeHistNFrecAxes, savenfrecx);
	PrintNFrec(&gd.histNFrecX1[1], cmd.sizeHistNFrecAxes, savenfrecx1);
	PrintNFrec(&gd.histNFrecX2[1], cmd.sizeHistNFrecAxes, savenfrecx2);
	PrintNFrec(&gd.histNFrecXD[1], cmd.sizeHistNFrecAxes, savenfrecxd);
	PrintNFrec(&gd.histNFrecY[1], cmd.sizeHistNFrecAxes, savenfrecy);
	PrintNFrec(&gd.histNFrecY1[1], cmd.sizeHistNFrecAxes, savenfrecy1);
	PrintNFrec(&gd.histNFrecY2[1], cmd.sizeHistNFrecAxes, savenfrecy2);
	PrintNFrec(&gd.histNFrecYD[1], cmd.sizeHistNFrecAxes, savenfrecyd);
#if (NDIM==3)
	PrintNFrec(&gd.histNFrecZ[1], cmd.sizeHistNFrecAxes, savenfrecz);
	PrintNFrec(&gd.histNFrecZ1[1], cmd.sizeHistNFrecAxes, savenfrecz1);
	PrintNFrec(&gd.histNFrecZ2[1], cmd.sizeHistNFrecAxes, savenfrecz2);
	PrintNFrec(&gd.histNFrecZD[1], cmd.sizeHistNFrecAxes, savenfreczd);
#endif
}

local void PrintNFrec(int *rn, int nval, char *filename)
{
    stream outstr_nfrec;
    char   buf[200];
	int *n, *fn, *p, *s, i;
	real vol;

	vol = gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]/cmd.sizeHistNFrecAxes;

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
				((real)n[s[i]])/( ((real)cmd.stepAvgNFrecAxes)*vol),fn[s[i]]);
	}
    fclose(outstr_nfrec);
	sprintf(buf,"mv %s %s",savenfrecaxestmp,filename);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);

}

#undef savenfrecx
#undef savenfrecx1
#undef savenfrecx2
#undef savenfrecxd
#undef savenfrecy
#undef savenfrecy1
#undef savenfrecy2
#undef savenfrecyd
#if (NDIM==3)
#undef savenfrecz
#undef savenfrecz1
#undef savenfrecz2
#undef savenfreczd
#endif
#undef savenfrecaxestmp

local void ComputeFrequencies(int *n, int *fn, int nval)
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

local void EvalVelDist(void)
{
    bodyptr p;
	real deltaV, histSum, vv;
	int j, k;
	double cpustart;

	printf("EvalVelDist: Entrando ... ");
	cpustart = cputime();                       

	gd.countVel = gd.countVel + 1;
	if (gd.countVel == 1) {
		for (j = 1; j <= cmd.sizeHistVel; j++) gd.histVel[j] = 0.0;
	}
	deltaV = cmd.rangeVel/cmd.sizeHistVel;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		vv =0.0;
		DO_COORD(k)
			vv += rsqr(Vel(p)[k]);
		j = (int) (rsqrt(vv)/deltaV)+1;
		if (j>cmd.sizeHistVel) j=cmd.sizeHistVel;
		gd.histVel[j] = gd.histVel[j]+1.0;
	}
	if (gd.countVel == cmd.stepAvgVel) {
		histSum = 0.0;
		for (j = 1; j <= cmd.sizeHistVel; j++)
			histSum= histSum+gd.histVel[j];
		for (j=1; j<=cmd.sizeHistVel; j++)
			gd.histVel[j]= gd.histVel[j]/histSum;
		PrintVelDist(stdout);
		gd.countVel = 0.0;
	}

	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define saveveltmp "vel.tmp"
#define savevel "vel.dat"

local void PrintVelDist(FILE *fp)
{
	real vBin;
	int n;
    stream outstr_vel;
    char   buf[200];

    outstr_vel = stropen(saveveltmp, "w!");
	
	for (n=1; n<=cmd.sizeHistVel; n++) {
		vBin = (n-0.5)*cmd.rangeVel/cmd.sizeHistVel;
		fprintf(outstr_vel,"%8.3f %8.3f\n",vBin,gd.histVel[n]);
	}
    fclose(outstr_vel);
	sprintf(buf,"mv %s %s",saveveltmp,savevel);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);
}

#undef saveveltmp
#undef savevel

local void EvalRdf(void)
{
    bodyptr j1, j2;
	real deltaR, normFac, rr, rrRange, Vol;
	int k, n;
	vector dr;
	double cpustart;
	real normFac11, normFac12, normFac21, normFac22;

	printf("EvalRdf: Entrando ... ");

	cpustart = cputime();                       

	gd.countRdf = gd.countRdf + 1;
	if (gd.countRdf == 1) {
		for (n = 1; n <= cmd.sizeHistRdf; n++) gd.histRdf[n] = 0.0;
		for (n = 1; n <= cmd.sizeHistRdf; n++) gd.histRdf11[n] = 0.0;
		for (n = 1; n <= cmd.sizeHistRdf; n++) gd.histRdf12[n] = 0.0;
		for (n = 1; n <= cmd.sizeHistRdf; n++) gd.histRdf22[n] = 0.0;
	}
	rrRange = rsqr(cmd.rangeRdf);
	deltaR = cmd.rangeRdf/cmd.sizeHistRdf;

	DO_BODY(j1, bodytab, bodytab+gd.nbody-1) {
		DO_BODY(j2, j1+1, bodytab+gd.nbody) {
			if (j1==j2) continue;
			DO_COORD(k) {
				dr[k] = Pos(j1)[k] - Pos(j2)[k];
				dr[k]=dr[k]-((real)(nint(dr[k]/gdforce.Box[k])))*gdforce.Box[k];
			}
			rr=0.0;
			DO_COORD(k)
				rr += rsqr(dr[k]);
			if (rr < rrRange) {
				n = (int) (rsqrt(rr) / deltaR) + 1;
				gd.histRdf[n] = gd.histRdf[n] + 1.;
				if ( Type(j1) == Type(j2) )
					if (Type(j1) == BODY1)
						gd.histRdf11[n] = gd.histRdf11[n] + 1.;
					else
						gd.histRdf22[n] = gd.histRdf22[n] + 1.;
				else
						gd.histRdf12[n] = gd.histRdf12[n] + 1.;
			}
		}
	}

	if (gd.countRdf == cmd.stepAvgRdf) {
		Vol = 1.0;
		DO_COORD(k)
			Vol = Vol*gdforce.Box[k];
#if (NDIM==3)
		normFac = Vol/(2.0*PI*rpow(deltaR,3.0)*gd.nbody*gd.nbody*gd.countRdf);
		if (gd.nbody1 != 0)
			normFac11 = Vol/(2.0*PI*rpow(deltaR,3.0)*gd.nbody1*gd.nbody1*gd.countRdf);
		else
			normFac11=0.;

		if (gd.nbody1 != 0 && gd.nbody2 != 0)
			normFac12 = Vol/(2.0*PI*rpow(deltaR,3.0)*gd.nbody1*gd.nbody2*gd.countRdf);
		else
			normFac12=0.;

		if (gd.nbody2 != 0)
			normFac22 = Vol/(2.0*PI*rpow(deltaR,3.0)*gd.nbody2*gd.nbody2*gd.countRdf);
		else
			normFac22=0.;
#else
		normFac = Vol/(PI*rpow(deltaR,2.0)*gd.nbody*gd.nbody*gd.countRdf);
		if (gd.nbody1 != 0)
			normFac11 = Vol/(PI*rpow(deltaR,2.0)*gd.nbody1*gd.nbody1*gd.countRdf);
		else
			normFac11=0.;

		if (gd.nbody2 != 0 && gd.nbody1 != 0)
			normFac12 = Vol/(PI*rpow(deltaR,2.0)*gd.nbody1*gd.nbody2*gd.countRdf);
		else
			normFac12=0.;

		if (gd.nbody2 != 0)
			normFac22 = Vol/(PI*rpow(deltaR,2.0)*gd.nbody2*gd.nbody2*gd.countRdf);
		else
			normFac22=0.;
#endif

		for (n = 1; n <= cmd.sizeHistRdf; n++) {
#if (NDIM==3)
			gd.histRdf[n] = gd.histRdf[n] * normFac / rsqr((int)n-0.5);
			gd.histRdf11[n] = gd.histRdf11[n] * normFac11 / rsqr((int)n-0.5);
			gd.histRdf12[n] = gd.histRdf12[n] * normFac12 / rsqr((int)n-0.5);
			gd.histRdf22[n] = gd.histRdf22[n] * normFac22 / rsqr((int)n-0.5);
//			gd.histRdf22[n] = gd.histRdf22[n] * normFac22 / 
//				( rsqr((real)n-1.0)+((real)n) - 1.0 + 1.0/3.0 );
#else
			gd.histRdf[n] = gd.histRdf[n] * normFac / ((int)n-0.5);
			gd.histRdf11[n] = gd.histRdf11[n] * normFac11 / ((int)n-0.5);
			gd.histRdf12[n] = gd.histRdf12[n] * normFac12 / ((int)n-0.5);
			gd.histRdf22[n] = gd.histRdf22[n] * normFac22 / ((int)n-0.5);
#endif
		}

		PrintRdf(stdout);
		gd.countRdf = 0;
	}
	
	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define saverdftmp "rdf.tmp"
#define saverdf "rdf.dat"
#define saveudrtmp "udr.tmp"
#define saveudr "udr.dat"
#define savessftmp "ssf.tmp"
#define savessf "ssf.dat"

local void PrintRdf(FILE *fp)
{
	real rBin;
	int n;
    stream outstr_rdf;
    char   buf[200];
	real ssq,fac,dr2,rri,rri3,UN,Pressure;
	realptr UNIntegrand, PressureIntegrand;
	int i, sizeHistK;
	real k, rangeK; 
	realptr SSFIntegrand, SSF;

    outstr_rdf = stropen(saverdftmp, "w!");

/*
	AllocMem (UNIntegrand, cmd.sizeHistRdf, real);
	AllocMem (PressureIntegrand, cmd.sizeHistRdf, real);

	ssq = cmd.sigma11*cmd.sigma11;						// Una sola especie...

	for (n=1; n<=cmd.sizeHistRdf; n++) {
		rBin = ((real)n-0.5)*cmd.rangeRdf/cmd.sizeHistRdf;
		dr2 = rBin*rBin;
		rri=ssq/dr2; rri3=rri*rri*rri;
		UNIntegrand[n] = dr2*(rri3-1.0)*rri3*gd.histRdf[n];
		PressureIntegrand[n] = dr2*rBin*rri3*(rri3-0.5)*rri*gd.histRdf[n];
	}
	fac = 2.0*PI*cmd.density*gdforce.fphi11;
	UN = fac * Integrate(UNIntegrand, cmd.sizeHistRdf);
	fac = (2.0/3.0)*PI*cmd.density*cmd.density*gdforce.fa11;
	Pressure = cmd.density*cmd.temperature 
				- fac * Integrate(PressureIntegrand, cmd.sizeHistRdf);

	fprintf(outstr_rdf,"# U/N, Pressure : %g %g\n",UN,Pressure);
	free(UNIntegrand);
	free(PressureIntegrand);
*/
	fprintf(outstr_rdf,"# r\t Rdf\t Rdf11\t Rdf12\t Rdf22\n");
	fprintf(outstr_rdf,"# <1>\t <2>\t <3>\t <4>\t <5>\n");
	for (n=1; n<=cmd.sizeHistRdf; n++) {
		rBin = ((real)n-0.5)*cmd.rangeRdf/cmd.sizeHistRdf;
		fprintf(outstr_rdf,"%8.4f %8.4f %8.4f %8.4f %8.4f\n", 
		rBin, gd.histRdf[n], gd.histRdf11[n], gd.histRdf12[n], gd.histRdf22[n]);
	}

    fclose(outstr_rdf);
	sprintf(buf,"mv %s %s",saverdftmp,saverdf);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);

/*
    outstr_rdf = stropen(saveudrtmp, "w!");
	for (n=1; n<=cmd.sizeHistRdf; n++) {
		rBin = ((real)n-0.5)*cmd.rangeRdf/cmd.sizeHistRdf;
		if (gd.histRdf[n] != 0)
			fprintf(outstr_rdf,"%8.4f %8.4f\n",
				rBin, -cmd.temperature*rlog(gd.histRdf[n])/LNE );
	}

    fclose(outstr_rdf);
	sprintf(buf,"mv %s %s",saveudrtmp,saveudr);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);
*/

/*
	sizeHistK = 50;
	AllocMem (SSF, sizeHistK, real);
	AllocMem (SSFIntegrand, cmd.sizeHistRdf, real);
	rangeK = 2. * PI * gd.numUcell[0] / gdforce.Box[0];
    outstr_rdf = stropen(savessftmp, "w!");
	for (i=1; i<=sizeHistK; i++) {
		k = ((real)i-0.5)*rangeK/sizeHistK;
		for (n=1; n<=cmd.sizeHistRdf; n++) {
			rBin = ((int)n-0.5)*cmd.rangeRdf/cmd.sizeHistRdf;
			SSFIntegrand[n] = gd.histRdf[n]*rBin*rsin(k*rBin)/k;
		}
		SSF[i] = 1.0 + 4.0*PI*cmd.density
						*Integrate(SSFIntegrand, cmd.sizeHistRdf);
	}

	for (i=1; i<=sizeHistK; i++) {
		k = ((int)i-0.5)*rangeK/sizeHistK;
		fprintf(outstr_rdf,"%8.4f %8.4f\n", k, SSF[i]);
	}

    fclose(outstr_rdf);
	sprintf(buf,"mv %s %s",savessftmp,savessf);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);
	free(SSF);
	free(SSFIntegrand);
*/
}

#undef saverdftmp
#undef saverdf
#undef saveudrtmp
#undef saveudr
#undef savessftmp
#undef savessf

// COMIENZO DE LA EVALUACION DEL PERFIL DE PRESIONES ---------------------------

local void EvalPressAxes(void)
{
    bodyptr p;
	real deltaX, deltaY, histSum;
	int i, l, j;
	double cpustart;
	real rangeX, rangeY;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
	real deltaZ, rangeZ;
#endif
	real vol;

	printf("EvalPressAxes: Entrando ... ");
	cpustart = cputime();                       

	rangeX = gdforce.Box[0];
	rangeY = gdforce.Box[1];
#if (NDIM==3)
	rangeZ = gdforce.Box[2];
#endif

	gd.countRhoAxes = gd.countRhoAxes + 1;
	if (gd.countRhoAxes == 1) {
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoX[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoY[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoX1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoY1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoX2[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoY2[j] = 0.0;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoZ[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoZ1[j] = 0.0;
		for (j = 1; j <= cmd.sizeHistRhoAxes; j++) gd.histRhoZ2[j] = 0.0;
#endif
	}
	deltaX = rangeX/cmd.sizeHistRhoAxes;
	deltaY = rangeY/cmd.sizeHistRhoAxes;
#if (NDIM==3)										// CORRECCION 2-3D CHECAR
	deltaZ = rangeZ/cmd.sizeHistRhoAxes;
#endif

	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		i = (int) ( (Pos(p)[0]+0.5*gdforce.Box[0] )/deltaX)+1;
		l = (int) ( (Pos(p)[1]+0.5*gdforce.Box[1] )/deltaY)+1;
#if (NDIM==3)
		j = (int) ( (Pos(p)[2]+0.5*gdforce.Box[2] )/deltaZ)+1;
#endif
		if (i>cmd.sizeHistRhoAxes) i=cmd.sizeHistRhoAxes;
			gd.histRhoX[i] = gd.histRhoX[i]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoX1[i] = gd.histRhoX1[i]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoX2[i] = gd.histRhoX2[i]+1.0;

		if (l>cmd.sizeHistRhoAxes) l=cmd.sizeHistRhoAxes;
			gd.histRhoY[l] = gd.histRhoY[l]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoY1[l] = gd.histRhoY1[l]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoY2[l] = gd.histRhoY2[l]+1.0;

#if (NDIM==3)										// CORRECCION 2-3D CHECAR
		if (j>cmd.sizeHistRhoAxes) j=cmd.sizeHistRhoAxes;
			gd.histRhoZ[j] = gd.histRhoZ[j]+1.0;
		if ( Type(p) == BODY1 )
			gd.histRhoZ1[j] = gd.histRhoZ1[j]+1.0;
		if ( Type(p) == BODY2 )
			gd.histRhoZ2[j] = gd.histRhoZ2[j]+1.0;
#endif
	}

	if (gd.countRhoAxes == cmd.stepAvgRhoAxes) {
		vol = gdforce.Box[0]*gdforce.Box[1]*deltaX;
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoAxes; j++)
			histSum= histSum+gd.histRhoX[j];
		for (j=1; j<=cmd.sizeHistRhoAxes; j++) {
			gd.histRhoX[j]= gd.histRhoX[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoX1[j]= gd.histRhoX1[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoX2[j]= gd.histRhoX2[j]/((real)cmd.stepAvgRhoAxes*vol);
		}

		vol = gdforce.Box[0]*gdforce.Box[1]*deltaY;
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoAxes; j++)
			histSum= histSum+gd.histRhoY[j];
		for (j=1; j<=cmd.sizeHistRhoAxes; j++) {
			gd.histRhoY[j]= gd.histRhoY[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoY1[j]= gd.histRhoY1[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoY2[j]= gd.histRhoY2[j]/((real)cmd.stepAvgRhoAxes*vol);
		}

#if (NDIM==3)
		vol = gdforce.Box[0]*gdforce.Box[1]*deltaZ;
		histSum = 0.0;
		for (j = 1; j < cmd.sizeHistRhoAxes; j++)
			histSum= histSum+gd.histRhoZ[j];
		for (j=1; j<=cmd.sizeHistRhoAxes; j++) {
			gd.histRhoZ[j]= gd.histRhoZ[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoZ1[j]= gd.histRhoZ1[j]/((real)cmd.stepAvgRhoAxes*vol);
			gd.histRhoZ2[j]= gd.histRhoZ2[j]/((real)cmd.stepAvgRhoAxes*vol);
		}
#endif
		PrintRhoAxes(stdout);
		gd.countRhoAxes = 0.0;
	}

	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define savepressaxestmp "pressaxes.tmp"
#define savepressaxes "pressaxes.dat"
#define PRESSFMT \
"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n"

local void PrintPressAxes(FILE *fp)
{
	real xBin, yBin;
	int n;
    stream outstr_rhoz;
    char   buf[200];
	real rangeX, rangeY;
#if (NDIM==3)
	real zBin, rangeZ;
#endif

    outstr_rhoz = stropen(savepressaxestmp, "w!");

	rangeX = gdforce.Box[0];
	rangeY = gdforce.Box[1];
#if (NDIM==3)
	rangeZ = gdforce.Box[2];
#endif

	for (n=1; n<=cmd.sizeHistRhoAxes; n++) {
		xBin = (n-0.5)*rangeX/cmd.sizeHistRhoAxes-0.5*gdforce.Box[0];
		yBin = (n-0.5)*rangeY/cmd.sizeHistRhoAxes-0.5*gdforce.Box[1];
#if (NDIM==3)
		zBin = (n-0.5)*rangeZ/cmd.sizeHistRhoAxes-0.5*gdforce.Box[2];
		fprintf(outstr_rhoz, PRESSFMT,
//		"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
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
	sprintf(buf,"mv %s %s",savepressaxestmp,savepressaxes);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);
}

#undef savepressaxes
#undef savepressaxestmp
#undef PRESSFMT

// FIN DE LA EVALUACION DEL PERFIL DE PRESIONES --------------------------------


// COMIENZO DE LA EVALUACION DEL COEFICIENTE DE DIFUSION -----------------------

local void EvalDiffusion (void)
{
  vector dr;
  int n, nb, ni;
  bodyptr p;

	for (nb = 0; nb < cmd.nBuffDiffuse; nb ++) {
		if (tBufD[nb].count == 0) {
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				SETV(tBufD[nb].orgR[n], Pos(p));
				SETV(tBufD[nb].rTrue[n], Pos(p));
			}
		}
		if (tBufD[nb].count >= 0) {
			ni = tBufD[nb].count;
			tBufD[nb].rrDiffuse[ni] = 0.;
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				VSub (dr, tBufD[nb].rTrue[n], Pos(p));
				VDiv (dr, dr, gdforce.Box);
				dr[0] = Nint (dr[0]);
				dr[1] = Nint (dr[1]);
#if (NDIM==3)
				dr[2] = Nint (dr[2]);
#endif
				VMul (dr, dr, gdforce.Box);
				VAdd (tBufD[nb].rTrue[n], Pos(p), dr);
				VSub (dr, tBufD[nb].rTrue[n], tBufD[nb].orgR[n]);
				tBufD[nb].rrDiffuse[ni] += VLenSq (dr);
			}
		}
		++ tBufD[nb].count;
	}
	AccumDiffusion ();
}

local void AccumDiffusion (void)
{
	real fac;
	int j, nb;

	for (nb = 0; nb < cmd.nBuffDiffuse; nb ++) {
		if (tBufD[nb].count == cmd.nValDiffuse) {
			for (j = 0; j < cmd.nValDiffuse; j ++)
				gd.rrDiffuseAv[j] += tBufD[nb].rrDiffuse[j];
			tBufD[nb].count = 0;
			++ gd.countDiffuseAv;
			if (gd.countDiffuseAv == cmd.stepAvgDiffuse) {
				fac = 1. / (NDIM * 2 * gd.nbody * cmd.stepDiffuse *
				gd.dtime * cmd.stepAvgDiffuse);
				for (j = 1; j < cmd.nValDiffuse; j ++)
					gd.rrDiffuseAv[j] *= fac / j;
				PrintDiffusion (stdout);
				ZeroDiffusion ();
			}
		}
	}
}

local void InitDiffusion (void)
{
	int nb;

	for (nb = 0; nb < cmd.nBuffDiffuse; nb ++)
		tBufD[nb].count = - nb * cmd.nValDiffuse / cmd.nBuffDiffuse;
	ZeroDiffusion ();
}

local void ZeroDiffusion (void)
{
	int j;

	gd.countDiffuseAv = 0;
	for (j = 0; j < cmd.nValDiffuse; j ++) gd.rrDiffuseAv[j] = 0.;
}

#define savediffusiontmp "diffusion.tmp"
#define savediffusion "diffusion.dat"

local void PrintDiffusion (FILE *fp)
{
	real tVal;
	int j;
    stream outstr;
    char   buf[200];
    stream outstr_transc;

    outstr_transc = stropen(outfile_transc, "a");

    outstr = stropen(savediffusiontmp, "w!");

	fprintf (outstr, "# diffusion\n");
	for (j = 0; j < cmd.nValDiffuse; j ++) {
		tVal = j * cmd.stepDiffuse * gd.dtime;
		fprintf (outstr, "%8.4f %8.4f\n", tVal, gd.rrDiffuseAv[j]);
	}
	fflush (outstr);
    fclose(outstr);
	sprintf(buf,"mv %s %s",savediffusiontmp,savediffusion);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);

	fprintf(outstr_transc,"%d %g %g\n", gd.nstep, gd.tnow,
		gd.rrDiffuseAv[cmd.nValDiffuse-1]);
	fclose(outstr_transc);
}

#undef savediffusiontmp
#undef savediffusion

// FIN DE LA EVALUACION DEL COEFICIENTE DE DIFUSION ----------------------------


// COMIENZO DE LA EVALUACION DE LA FUNCION DE AUTOCORRELACION DE VELOCIDADES ---

local void EvalVelAcf (void)
{
	int n, nb, ni;
	bodyptr p;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBufVAcf[nb].count == 0) {
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				SETV(tBufVAcf[nb].orgVel[n], Vel(p));
			}
		}
		if (tBufVAcf[nb].count >= 0) {
			ni = tBufVAcf[nb].count;
			tBufVAcf[nb].acfVel[ni] = 0.;
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				tBufVAcf[nb].acfVel[ni] += VDot (tBufVAcf[nb].orgVel[n], Vel(p));
			}
		}
		++ tBufVAcf[nb].count;
	}
	AccumVelAcf ();
}

local void AccumVelAcf (void)
{
	real fac;
	int j, nb;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBufVAcf[nb].count == cmd.nValAcf) {
			for (j = 0; j < cmd.nValAcf; j ++) gd.avAcfVel[j] += tBufVAcf[nb].acfVel[j];
			tBufVAcf[nb].count = 0;
			++ gd.countAcfAv;
			if (gd.countAcfAv == cmd.stepAvgAcf) {
				if (!cmd.adjustTemperature &&  cmd.adjustCenterOfMass)
					fac = cmd.stepAcf * gd.dtime / (NDIM * (gd.nbody-1) * cmd.stepAvgAcf);
				else
					fac = cmd.stepAcf * gd.dtime / (NDIM * gd.nbody * cmd.stepAvgAcf);
				gd.intAcfVel = fac * IntegrateVelAcf (gd.avAcfVel, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) gd.avAcfVel[j] /= gd.avAcfVel[0];
				gd.avAcfVel[0] = 1.;
				PrintVelAcf (stdout);
				ZeroVelAcf ();
			}
		}
	}
}

local void InitVelAcf (void)
{
	int nb;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++)
		tBufVAcf[nb].count = - nb * cmd.nValAcf / cmd.nBuffAcf;
	ZeroVelAcf ();
}

local void ZeroVelAcf (void)
{
	int j;

	gd.countAcfAv = 0;
	for (j = 0; j < cmd.nValAcf; j ++) gd.avAcfVel[j] = 0.;
}

#define savevelacftmp "velacf.tmp"
#define savevelacf "velacf.dat"

local void PrintVelAcf (FILE *fp)
{
	real tVal;
	int j;
    stream outstr;
    char   buf[200];
    stream outstr_transc;

    outstr_transc = stropen(outfile_transc, "a");

    outstr = stropen(savevelacftmp, "w!");


	fprintf (outstr, "# vel acf\n");
	for (j = 0; j < cmd.nValAcf; j ++) {
		tVal = j * cmd.stepAcf * gd.dtime;
		fprintf (outstr, "%8.4f %8.4f\n", tVal, gd.avAcfVel[j]);
	}
	fprintf (outstr, "# vel acf integral: %8.3f\n", gd.intAcfVel);
	fflush (outstr);
    fclose(outstr);
	sprintf(buf,"mv %s %s",savevelacftmp,savevelacf);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);

	fprintf(outstr_transc,"%d %g %g\n", gd.nstep, gd.tnow,
		gd.intAcfVel);
	fclose(outstr_transc);
}

#undef savevelacftmp
#undef savevelacf

local real IntegrateVelAcf (real *f, int nf)
{
  real s;
  int i;

  s = 0.5 * (f[0] + f[nf - 1]);
  for (i = 1; i < nf - 1; i ++) s += f[i];
  return (s);
}

// FIN DE LA EVALUACION DE LA FUNCION DE AUTOCORRELACION DE VELOCIDADES --------


// COMIENZO DE LA EVALUACION DE LA FUNCION DE AUTOCORRELACION DE PRESIONES (BULK-VISCOSITY)-----

local void EvaldPAcf (void)
{
	int n, nb, ni;
	bodyptr p;
	real dP;

	dP = gd.pressure-gd.sPressureSave;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBufdPAcf[nb].count == 0) {
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				tBufdPAcf[nb].orgdP[n]=dP;
			}
		}
		if (tBufdPAcf[nb].count >= 0) {
			ni = tBufdPAcf[nb].count;
			tBufdPAcf[nb].acfdP[ni] = 0.;
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				tBufdPAcf[nb].acfdP[ni] += tBufdPAcf[nb].orgdP[n] * dP;
			}
		}
		++ tBufdPAcf[nb].count;
	}
	AccumdPAcf ();
}

local void AccumdPAcf (void)
{
	real fac;
	int j, nb;
	real vol;
	int k;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBufdPAcf[nb].count == cmd.nValAcf) {
			for (j = 0; j < cmd.nValAcf; j ++) gd.avAcfdP[j] += tBufdPAcf[nb].acfdP[j];
			tBufdPAcf[nb].count = 0;
			++ gd.countAcfAv;
			if (gd.countAcfAv == cmd.stepAvgAcf) {
				vol=1.0;
				DO_COORD(k) vol *= gdforce.Box[k];
//				fac = vol * cmd.stepAcf * gd.dtime 
//						/ (cmd.temperature * gd.nbody * cmd.stepAvgAcf);
				fac = vol * cmd.stepAcf * gd.dtime 
						/ (cmd.temperature * cmd.stepAvgAcf);
				gd.intAcfdP = fac * IntegratedPAcf (gd.avAcfdP, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) gd.avAcfdP[j] /= gd.avAcfdP[0];
				gd.avAcfdP[0] = 1.;
				PrintdPAcf (stdout);
				ZerodPAcf ();
			}
		}
	}
}

local void InitdPAcf (void)
{
	int nb;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++)
		tBufdPAcf[nb].count = - nb * cmd.nValAcf / cmd.nBuffAcf;
	ZerodPAcf ();
}

local void ZerodPAcf (void)
{
	int j;

	gd.countAcfAv = 0;
	for (j = 0; j < cmd.nValAcf; j ++) gd.avAcfdP[j] = 0.;
}

#define savedpacftmp "dpacf.tmp"
#define savedpacf "dpacf.dat"

local void PrintdPAcf (FILE *fp)
{
	real tVal;
	int j;
    stream outstr;
    char   buf[200];
    stream outstr_transc;

    outstr_transc = stropen(outfile_transc, "a");

    outstr = stropen(savedpacftmp, "w!");

	fprintf (outstr, "# dp acf\n");
	for (j = 0; j < cmd.nValAcf; j ++) {
		tVal = j * cmd.stepAcf * gd.dtime;
		fprintf (outstr, "%8.4f %8.4f\n", tVal, gd.avAcfdP[j]);
	}
	fprintf (outstr, "# dp acf integral: %8.3f\n", gd.intAcfdP);
	fflush (outstr);
    fclose(outstr);
	sprintf(buf,"mv %s %s",savedpacftmp,savedpacf);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);

	fprintf(outstr_transc,"%d %g %g\n", gd.nstep, gd.tnow,
		gd.intAcfdP);
	fclose(outstr_transc);
}

#undef savedpacftmp
#undef savedpacf

local real IntegratedPAcf (real *f, int nf)
{
  real s;
  int i;

  s = 0.5 * (f[0] + f[nf - 1]);
  for (i = 1; i < nf - 1; i ++) s += f[i];
  return (s);
}

// FIN DE LA EVALUACION DE LA FUNCION DE AUTOCORRELACION DE PRESIONES (BULK-VISCOSITY)----------


// COMIENZO DE LA EVALUACION DE LAS PROPIEDADES DE TRANSPORTE ------------------

local void EvalVacf (void)
{
	vector vecTherm, vecVisc;
	int n, nb, ni;
	bodyptr p;
	double cpustart;
	int a, b;
// - Bulk Viscosity
	vector JK, JP;
// - Shear Viscosity
	vector vecVisc_K, vecVisc_P;
//

	fprintf(stdout,"EvalVacf: Entrando ... %d %g ...", gd.nstep, gd.tnow);

	cpustart = cputime();
// - Bulk Viscosity
	VZero (JK);
	VZero (JP);

// - Shear Viscosity
	VZero (vecVisc_K);
	VZero (vecVisc_P);
//
	VZero (vecVisc);
	VZero (vecTherm);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {	// HACER LA VERSION 2D ...
// - Bulk Viscosity
		JK[0] += Mass(p) * Vel(p)[0] * Vel(p)[0];
		JK[1] += Mass(p) * Vel(p)[1] * Vel(p)[1];
		JK[2] += Mass(p) * Vel(p)[2] * Vel(p)[2];
		JP[0] += 0.5*rf(p)[0][0];
		JP[1] += 0.5*rf(p)[1][1];
		JP[2] += 0.5*rf(p)[2][2];

// - Shear Viscosity
		vecVisc_K[0] += Mass(p) * Vel(p)[1] * Vel(p)[2];
		vecVisc_K[1] += Mass(p) * Vel(p)[2] * Vel(p)[0];
		vecVisc_K[2] += Mass(p) * Vel(p)[0] * Vel(p)[1];

		vecVisc_P[0] += 0.5 * rf(p)[1][2];
		vecVisc_P[1] += 0.5 * rf(p)[2][0];
		vecVisc_P[2] += 0.5 * rf(p)[0][1];
//
		vecVisc[0] += Mass(p) * Vel(p)[1] * Vel(p)[2] + 0.5 * rf(p)[1][2];
		vecVisc[1] += Mass(p) * Vel(p)[2] * Vel(p)[0] + 0.5 * rf(p)[2][0];
		vecVisc[2] += Mass(p) * Vel(p)[0] * Vel(p)[1] + 0.5 * rf(p)[0][1];

		en(p) += VLenSq (Vel(p));
		VVSAdd (vecTherm, 0.5 * en(p), Vel(p));

		vecTherm[0] += 0.5 * VDot(Vel(p), rf(p)[0]);
		vecTherm[1] += 0.5 * VDot(Vel(p), rf(p)[1]);
		vecTherm[2] += 0.5 * VDot(Vel(p), rf(p)[2]);
	}
// - Bulk Viscosity
	JK[0] -= gd.PV_KAvg;
	JK[1] -= gd.PV_KAvg;
	JK[2] -= gd.PV_KAvg;
	JP[0] -= (1.0/3.0)*gd.PV_PAvg;
	JP[1] -= (1.0/3.0)*gd.PV_PAvg;
	JP[2] -= (1.0/3.0)*gd.PV_PAvg;

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBuf[nb].count == 0) {
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				SETV(tBuf[nb].orgVel[n], Vel(p));
			}
		}

		if (tBuf[nb].count >= 0) {
			ni = tBuf[nb].count;
			tBuf[nb].acfVel[ni] = 0.;
			DO_BODY(p, bodytab, bodytab+gd.nbody) {
				n = p-bodytab;
				tBuf[nb].acfVel[ni] += VDot(tBuf[nb].orgVel[n], Vel(p));
			}
		}

		if (tBuf[nb].count == 0) {
// - Bulk Viscosity
			SETV(tBuf[nb].orgBVisc_K, JK);
			SETV(tBuf[nb].orgBVisc_P, JP);
//
			SETV(tBuf[nb].orgVisc_K, vecVisc_K);
			SETV(tBuf[nb].orgVisc_P, vecVisc_P);
//
			SETV(tBuf[nb].orgVisc, vecVisc);
			SETV(tBuf[nb].orgTherm, vecTherm);
		}
// - Bulk Viscosity
//		tBuf[nb].acfBVisc_KK[ni] = VDot(tBuf[nb].orgBVisc_K, JK);
//		tBuf[nb].acfBVisc_KP[ni] = VDot(tBuf[nb].orgBVisc_K, JP)
//									+ VDot(tBuf[nb].orgBVisc_P, JK);
//		tBuf[nb].acfBVisc_PP[ni] = VDot(tBuf[nb].orgBVisc_P, JP);
		tBuf[nb].acfBVisc_KK[ni] = 0.0;
		tBuf[nb].acfBVisc_KP[ni] = 0.0;
		tBuf[nb].acfBVisc_PP[ni] = 0.0;
		DO_COORD(a) {
			DO_COORD(b) {
				tBuf[nb].acfBVisc_KK[ni] += tBuf[nb].orgBVisc_K[a] * JK[b];
				tBuf[nb].acfBVisc_KP[ni] += tBuf[nb].orgBVisc_K[a] * JP[b]
										 +  tBuf[nb].orgBVisc_P[a] * JK[b];
				tBuf[nb].acfBVisc_PP[ni] += tBuf[nb].orgBVisc_P[a] * JP[b];
			}
		}
		tBuf[nb].acfBVisc_KK[ni] /= 9.0;
		tBuf[nb].acfBVisc_KP[ni] /= 9.0;
		tBuf[nb].acfBVisc_PP[ni] /= 9.0;
//
		tBuf[nb].acfVisc_KK[ni] = VDot(tBuf[nb].orgVisc_K, vecVisc_K);
		tBuf[nb].acfVisc_KP[ni] = VDot(tBuf[nb].orgVisc_K, vecVisc_P)
									+ VDot(tBuf[nb].orgVisc_P, vecVisc_K);
		tBuf[nb].acfVisc_PP[ni] = VDot(tBuf[nb].orgVisc_P, vecVisc_P);
//
		tBuf[nb].acfVisc[ni] = VDot(tBuf[nb].orgVisc, vecVisc);
		tBuf[nb].acfTherm[ni] = VDot(tBuf[nb].orgTherm, vecTherm);
		++ tBuf[nb].count;
	}

	AccumVacf();

	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

local void AccumVacf(void)
{
	real fac;
	int j, nb;

// BORRAR
    fprintf(stdout,"AccumVacf: Entrando ... %d %g ...", gd.nstep, gd.tnow);

	for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
		if (tBuf[nb].count == cmd.nValAcf) {
			for (j = 0; j < cmd.nValAcf; j ++) 
				gd.avAcfVel[j] += tBuf[nb].acfVel[j];
			for (j = 0; j < cmd.nValAcf; j ++) {
// - Bulk Viscosity
				gd.avAcfBVisc_KK[j] += tBuf[nb].acfBVisc_KK[j];
				gd.avAcfBVisc_KP[j] += tBuf[nb].acfBVisc_KP[j];
				gd.avAcfBVisc_PP[j] += tBuf[nb].acfBVisc_PP[j];
// - Shear Viscosity
				gd.avAcfVisc_KK[j] += tBuf[nb].acfVisc_KK[j];
				gd.avAcfVisc_KP[j] += tBuf[nb].acfVisc_KP[j];
				gd.avAcfVisc_PP[j] += tBuf[nb].acfVisc_PP[j];
//
				gd.avAcfVisc[j] += tBuf[nb].acfVisc[j];
				gd.avAcfTherm[j] += tBuf[nb].acfTherm[j];
			}
			tBuf[nb].count = 0;
			++ gd.countAcfAv;
			if (gd.countAcfAv == cmd.stepAvgAcf) {
				fac = cmd.stepAcf * gd.dtime / (NDIM * gd.nbody*cmd.stepAvgAcf);
				gd.intAcfVel = fac * IntegrateVacf (gd.avAcfVel, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) 
					gd.avAcfVel[j] /= gd.avAcfVel[0];
				gd.avAcfVel[0] = 1.;

// - Bulk Viscosity
				fac = cmd.density * cmd.stepAcf * gd.dtime /
						(3. * cmd.temperature * gd.nbody * cmd.stepAvgAcf);
				gd.intAcfBVisc_KK = fac * IntegrateVacf (gd.avAcfBVisc_KK, cmd.nValAcf);
				gd.intAcfBVisc_KP = fac * IntegrateVacf (gd.avAcfBVisc_KP, cmd.nValAcf);
				gd.intAcfBVisc_PP = fac * IntegrateVacf (gd.avAcfBVisc_PP, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) {
					gd.avAcfBVisc_KK[j] /= gd.avAcfBVisc_KK[0];
					gd.avAcfBVisc_KP[j] /= gd.avAcfBVisc_KP[0];
					gd.avAcfBVisc_PP[j] /= gd.avAcfBVisc_PP[0];
				}
				gd.avAcfBVisc_KK[0] = 1.;
				gd.avAcfBVisc_KP[0] = 1.;
				gd.avAcfBVisc_PP[0] = 1.;
//

// - Shear Viscosity
				fac = cmd.density * cmd.stepAcf * gd.dtime /
						(3. * cmd.temperature * gd.nbody * cmd.stepAvgAcf);
//
				gd.intAcfVisc_KK = fac * IntegrateVacf (gd.avAcfVisc_KK, cmd.nValAcf);
				gd.intAcfVisc_KP = fac * IntegrateVacf (gd.avAcfVisc_KP, cmd.nValAcf);
				gd.intAcfVisc_PP = fac * IntegrateVacf (gd.avAcfVisc_PP, cmd.nValAcf);
//
				gd.intAcfVisc = fac * IntegrateVacf (gd.avAcfVisc, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) {
//
					gd.avAcfVisc_KK[j] /= gd.avAcfVisc_KK[0];
					gd.avAcfVisc_KP[j] /= gd.avAcfVisc_KP[0];
					gd.avAcfVisc_PP[j] /= gd.avAcfVisc_PP[0];
//
					gd.avAcfVisc[j] /= gd.avAcfVisc[0];
				}
//
				gd.avAcfVisc_KK[0] = 1.;
				gd.avAcfVisc_KP[0] = 1.;
				gd.avAcfVisc_PP[0] = 1.;
//
				gd.avAcfVisc[0] = 1.;


				fac = cmd.density * cmd.stepAcf * gd.dtime /
						(3.0*rsqr(cmd.temperature) * gd.nbody * cmd.stepAvgAcf);
				gd.intAcfTherm = fac * IntegrateVacf (gd.avAcfTherm, cmd.nValAcf);
				for (j = 1; j < cmd.nValAcf; j ++) 
					gd.avAcfTherm[j] /= gd.avAcfTherm[0];
				gd.avAcfTherm[0] = 1.;
				PrintVacf (stdout);
				ZeroVacf ();
			}
		}
	}
}

local void InitVacf(void)
{
	int nb;

// BORRAR
    fprintf(stdout,"InitVacf: Entrando ... %d %g ...", gd.nstep, gd.tnow);

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
// - Bulk Viscosity
		gd.avAcfBVisc_KK[j] = 0.;
		gd.avAcfBVisc_KP[j] = 0.;
		gd.avAcfBVisc_PP[j] = 0.;
		gd.PV_KAvg = gd.PV_KAvgSave;
		gd.PV_PAvg = gd.PV_PAvgSave;
// - Shear Viscosity
		gd.avAcfVisc_KK[j] = 0.;
		gd.avAcfVisc_KP[j] = 0.;
		gd.avAcfVisc_PP[j] = 0.;

		gd.avAcfVisc[j] = 0.;
		gd.avAcfTherm[j] = 0.;
	}
}

#define saveacftmp "acf.tmp"
#define saveacf "acf.dat"

local void PrintVacf(FILE *fp)
{
	real tVal;
	int j;
    stream outstr;
    char   buf[200];
    stream outstr_transc;
	real avAcfBVisc, avAcfBVisc0, intAcfBVisc;

// BORRAR
//    fprintf(stdout,"PrintVacf: Entrando ... %d %g ...", gd.nstep, gd.tnow);

    outstr_transc = stropen(outfile_transc, "a");

    outstr = stropen(saveacftmp, "w!");

	avAcfBVisc0 = gd.avAcfBVisc_KK[0]+gd.avAcfBVisc_KP[0]+gd.avAcfBVisc_PP[0];

	fprintf (outstr, "# acf : %d %g\n",gd.nstep,gd.tnow);
	for (j = 0; j < cmd.nValAcf; j ++) {
		tVal = j * cmd.stepAcf * gd.dtime;
		avAcfBVisc = (gd.avAcfBVisc_KK[j]+gd.avAcfBVisc_KP[j]
					+gd.avAcfBVisc_PP[j])/avAcfBVisc0;
		fprintf (outstr, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", tVal,
		gd.avAcfVel[j],
		avAcfBVisc, gd.avAcfBVisc_KK[j], gd.avAcfBVisc_KP[j], gd.avAcfBVisc_PP[j], 
		gd.avAcfVisc[j], gd.avAcfVisc_KK[j], gd.avAcfVisc_KP[j], gd.avAcfVisc_PP[j], 
		gd.avAcfTherm[j]);
	}

	intAcfBVisc = gd.intAcfBVisc_KK + gd.intAcfBVisc_KP + gd.intAcfBVisc_PP;

	fprintf(outstr,"# acf integrals: Diff, BVisc, BVisc_KK, BVisc_KP, BVisc_PP, Visc, Visc_KK, Visc_KP, Visc_PP, Therm : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
		gd.intAcfVel, 
		intAcfBVisc, gd.intAcfBVisc_KK, gd.intAcfBVisc_KP, gd.intAcfBVisc_PP, 
		gd.intAcfVisc, gd.intAcfVisc_KK, gd.intAcfVisc_KP, gd.intAcfVisc_PP, 
		gd.intAcfTherm);
	fflush(outstr);
    fclose(outstr);
	sprintf(buf,"mv %s %s",saveacftmp,saveacf);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);

	fprintf(outstr_transc,"%d %g %g %g %g %g %g %g %g %g %g %g\n", gd.nstep, gd.tnow,
		gd.intAcfVel, 
		intAcfBVisc, gd.intAcfBVisc_KK, gd.intAcfBVisc_KP, gd.intAcfBVisc_PP, 
		gd.intAcfVisc, gd.intAcfVisc_KK, gd.intAcfVisc_KP, gd.intAcfVisc_PP, 
		gd.intAcfTherm);
	fclose(outstr_transc);
}

#undef saveacftmp
#undef saveacf
#undef outfile_transc

local real IntegrateVacf(real *f, int nf)
{
  real s;
  int i;

  s = 0.5 * (f[0] + f[nf - 1]);
  for (i = 1; i < nf - 1; i ++) s += f[i];
  return (s);
}


// FIN DE LA EVALUACION DE LAS PROPIEDADES DE TRANSPORTE -----------------------

// COMIENZO DE LA EVALUACION DE LAS CORRELACIONES ESPACIO-TEMPORALES -----------

local void EvalSpacetimeCorr ()
{
  real b, c, c0, c1, c2, kVal, s, s1, s2, w;
  int j, k, m, n, nb, nc, nv;
  bodyptr p;

  for (j = 0; j < 24 * cmd.nFunCorr; j ++) gd.valST[j] = 0.;
  kVal = 2. * M_PI / gdforce.Box[0];		// Se supone que la caja es cubica
  DO_BODY(p, bodytab, bodytab+gd.nbody) {
    j = 0;
    for (k = 0; k < 3; k ++) {
      for (m = 0; m < cmd.nFunCorr; m ++) {
        if (m == 0) {
          b = kVal * VComp (Pos(p), k);
          c = rcos (b);
          s = rsin (b);
          c0 = c;
        } else if (m == 1) {
          c1 = c;
          s1 = s;
          c = 2. * c0 * c1 - 1.;
          s = 2. * c0 * s1;
        } else {
          c2 = c1;
          s2 = s1;
          c1 = c;
          s1 = s;
          c = 2. * c0 * c1 - c2;
          s = 2. * c0 * s1 - s2;
        }
        gd.valST[j ++] += Vel(p)[0] * c;
        gd.valST[j ++] += Vel(p)[0] * s;
        gd.valST[j ++] += Vel(p)[1] * c;
        gd.valST[j ++] += Vel(p)[1] * s;
#if (NDIM==3)
        gd.valST[j ++] += Vel(p)[2] * c;
        gd.valST[j ++] += Vel(p)[2] * s;
#endif
        gd.valST[j ++] += c;
        gd.valST[j ++] += s;
      }
    }
  }
  for (nb = 0; nb < cmd.nBuffCorr; nb ++) {
    if (tBufC[nb].count == 0) {
      for (j = 0; j < 24 * cmd.nFunCorr; j ++)
         tBufC[nb].orgST[j] = gd.valST[j];
    }
    if (tBufC[nb].count >= 0) {
      for (j = 0; j < 3 * cmd.nFunCorr; j ++)
         tBufC[nb].acfST[j][tBufC[nb].count] = 0.;
      j = 0;
      for (k = 0; k < 3; k ++) {
        for (m = 0; m < cmd.nFunCorr; m ++) {
          for (nc = 0; nc < 4; nc ++) {
            nv = 3 * m + 2;
            if (nc < 3) {
              w = rsqr(kVal * (m + 1));
              -- nv;
              if (nc == k) -- nv;
              else w *= 0.5;
            } else w = 1.;
            tBufC[nb].acfST[nv][tBufC[nb].count] +=
               w * (gd.valST[j] * tBufC[nb].orgST[j] +
               gd.valST[j + 1] * tBufC[nb].orgST[j + 1]);
            j += 2;
          }
        }
      }
    }
    ++ tBufC[nb].count;
  }
  AccumSpacetimeCorr ();
}

local void AccumSpacetimeCorr ()
{
  int j, n, nb;

  for (nb = 0; nb < cmd.nBuffCorr; nb ++) {
    if (tBufC[nb].count == cmd.nValCorr) {
      for (j = 0; j < 3 * cmd.nFunCorr; j ++) {
        for (n = 0; n < cmd.nValCorr; n ++)
           gd.avAcfST[j][n] += tBufC[nb].acfST[j][n];
      }
      tBufC[nb].count = 0;
      ++ gd.countCorrAv;
      if (gd.countCorrAv == cmd.stepAvgCorr) {
        for (j = 0; j < 3 * cmd.nFunCorr; j ++) {
          for (n = 0; n < cmd.nValCorr; n ++)
             gd.avAcfST[j][n] /= 3. * gd.nbody * cmd.stepAvgCorr;
        }
        PrintSpacetimeCorr (stdout);
        ZeroSpacetimeCorr ();
      }
    }
  }
}

local void InitSpacetimeCorr ()
{
  int nb;

  for (nb = 0; nb < cmd.nBuffCorr; nb ++)
     tBufC[nb].count = - nb * cmd.nValCorr / cmd.nBuffCorr;
  ZeroSpacetimeCorr ();
}

local void ZeroSpacetimeCorr ()
{
  int j, n;

  gd.countCorrAv = 0;
  for (j = 0; j < 3 * cmd.nFunCorr; j ++) {
    for (n = 0; n < cmd.nValCorr; n ++) gd.avAcfST[j][n] = 0.;
  }
}

#define savecorrtmp "corr.tmp"
#define savecorr "corr.dat"

local void PrintSpacetimeCorr (FILE *fp)
{
	real tVal;
	int j, k, n;
	char *header[] = {"cur-long", "cur-trans", "density"};
    stream outstr;
    char   buf[200];

    outstr = stropen(savecorrtmp, "w!");

	fprintf (outstr, "space-time corr\n");
	for (k = 0; k < 3; k ++) {
		fprintf (outstr, "%s\n", header[k]);
		for (n = 0; n < cmd.nValCorr; n ++) {
			tVal = n * cmd.stepCorr * gd.dtime;
			fprintf (outstr, "%7.3f", tVal);
			for (j = 0; j < cmd.nFunCorr; j ++)
				fprintf (outstr, " %8.4f", gd.avAcfST[3 * j + k][n]);
			fprintf (outstr, "\n");
		}
	}
    fclose(outstr);
	sprintf(buf,"mv %s %s",savecorrtmp,savecorr);
	fprintf(gd.outlog,"[nstep: %d] : system: %s ...\n", gd.nstep, buf);
	system(buf);
}

#undef savecorrtmp
#undef savecorr

// FIN DE LA EVALUACION DE LAS CORRELACIONES ESPACIO-TEMPORALES ----------------

local void savestate(string pattern)
{
    char namebuf[256];
    stream str;
    int nchars, ndim, nb;

    sprintf(namebuf, pattern, gd.nstep & 1);      
    str = stropen(namebuf, "w!");
    nchars = strlen(getargv0()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getargv0(), nchars * sizeof(char), str);
    nchars = strlen(getversion()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getversion(), nchars * sizeof(char), str);

    nchars = strlen(cmd.forcecalc_method) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.forcecalc_method, nchars * sizeof(char), str);

    safewrite(&cmd.potType, sizeof(int), str);

    safewrite(&cmd.intMethod, sizeof(int), str);

    nchars = strlen(cmd.dtimestr) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.dtimestr, nchars * sizeof(char), str);

    safewrite(&gd.dtime, sizeof(real), str);
    safewrite(&cmd.temperature, sizeof(real), str);
    safewrite(&cmd.adjustTemperature, sizeof(bool), str);
    safewrite(&cmd.adjustCenterOfMass, sizeof(bool), str);
    safewrite(&cmd.density, sizeof(real), str);

    safewrite(&gd.mass1, sizeof(real), str);
    safewrite(&gd.mass2, sizeof(real), str);
    safewrite(&gd.Lx, sizeof(real), str);
    safewrite(&gd.Ly, sizeof(real), str);
#if (NDIM==3)
    safewrite(&gd.Lz, sizeof(real), str);
#endif
    safewrite(&cmd.eps11, sizeof(real), str);
    safewrite(&cmd.eps12, sizeof(real), str);
    safewrite(&cmd.eps22, sizeof(real), str);
    safewrite(&cmd.sigma11, sizeof(real), str);
    safewrite(&cmd.sigma12, sizeof(real), str);
    safewrite(&cmd.sigma22, sizeof(real), str);
    safewrite(&cmd.Rcut11, sizeof(real), str);
    safewrite(&cmd.Rcut12, sizeof(real), str);
    safewrite(&cmd.Rcut22, sizeof(real), str);

    safewrite(&gd.nbody, sizeof(int), str);
    safewrite(&gd.nbody1, sizeof(int), str);
    safewrite(&gd.nbody2, sizeof(int), str);
    safewrite(&gd.nbc1, sizeof(int), str);
    safewrite(&gd.nbc2, sizeof(int), str);

    safewrite(&cmd.stepEquil, sizeof(int), str);
    safewrite(&cmd.stepAvg, sizeof(int), str);

    safewrite(&cmd.computeRhoAxes, sizeof(bool), str);
    safewrite(&cmd.stepAvgRhoAxes, sizeof(int), str);
    safewrite(&cmd.stepRhoAxes, sizeof(int), str);
    safewrite(&cmd.sizeHistRhoAxes, sizeof(int), str);

	if (cmd.computeRhoAxes) {
		safewrite(&gd.countRhoAxes, sizeof(int), str);
		safewrite(gd.histRhoX, cmd.sizeHistRhoAxes * sizeof(real), str);
		safewrite(gd.histRhoY, cmd.sizeHistRhoAxes * sizeof(real), str);
		safewrite(gd.histRhoX1, cmd.sizeHistRhoAxes * sizeof(real), str);
		safewrite(gd.histRhoY1, cmd.sizeHistRhoAxes * sizeof(real), str);
		safewrite(gd.histRhoX2, cmd.sizeHistRhoAxes * sizeof(real), str);
		safewrite(gd.histRhoY2, cmd.sizeHistRhoAxes * sizeof(real), str);
#if (NDIM==3)
		safewrite(gd.histRhoZ, cmd.sizeHistRhoAxes * sizeof(real), str);
		safewrite(gd.histRhoZ1, cmd.sizeHistRhoAxes * sizeof(real), str);
		safewrite(gd.histRhoZ2, cmd.sizeHistRhoAxes * sizeof(real), str);
#endif
	}

    safewrite(&cmd.computeNFrecAxes, sizeof(bool), str);
    safewrite(&cmd.stepAvgNFrecAxes, sizeof(int), str);
    safewrite(&cmd.stepNFrecAxes, sizeof(int), str);
    safewrite(&cmd.sizeHistNFrecAxes, sizeof(int), str);

	if (cmd.computeNFrecAxes) {
		safewrite(&gd.countNFrecAxes, sizeof(int), str);
		safewrite(gd.histNFrecX, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecY, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecX1, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecY1, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecX2, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecY2, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecXD, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecYD, cmd.sizeHistNFrecAxes * sizeof(real), str);
#if (NDIM==3)
		safewrite(gd.histNFrecZ, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecZ1, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecZ2, cmd.sizeHistNFrecAxes * sizeof(real), str);
		safewrite(gd.histNFrecZD, cmd.sizeHistNFrecAxes * sizeof(real), str);
#endif
	}

	if (cmd.computeVelDist)
        gd.stateVel = TRUE;
    else
        gd.stateVel = FALSE;
    
    safewrite(&gd.stateVel, sizeof(bool), str);
    safewrite(&cmd.computeVelDist, sizeof(bool), str);
    safewrite(&cmd.stepAvgVel, sizeof(int), str);
    safewrite(&cmd.stepVel, sizeof(int), str);
    safewrite(&cmd.sizeHistVel, sizeof(int), str);
    safewrite(&cmd.rangeVel, sizeof(real), str);

	if (cmd.computeVelDist) {
		safewrite(&gd.countVel, sizeof(int), str);
		safewrite(gd.histVel, cmd.sizeHistVel * sizeof(real), str);
	}

    // Una vez fijos rangeRdf y sizeHistRdf no se podran modificar ...
    // Falta tener algo parecido para rangeRdf y sizeHistRdf
    // rangeRdf si funciona ... pero hay que dejar pasar la escritura
    // del archivo rdf.dat varias veces para que se estabilize el nuevo
    // valor de rangeRdf ...
	if (cmd.computeRdf)
        gd.stateRdf = TRUE;
    else
        gd.stateRdf = FALSE;
        
    safewrite(&gd.stateRdf, sizeof(bool), str);
    safewrite(&cmd.computeRdf, sizeof(bool), str);
    safewrite(&cmd.stepAvgRdf, sizeof(int), str);
    safewrite(&cmd.stepRdf, sizeof(int), str);
    safewrite(&cmd.sizeHistRdf, sizeof(int), str);
    safewrite(&cmd.rangeRdf, sizeof(real), str);

	if (cmd.computeRdf) {
		safewrite(&gd.countRdf, sizeof(int), str);
		safewrite(gd.histRdf, cmd.sizeHistRdf * sizeof(real), str);
		safewrite(gd.histRdf11, cmd.sizeHistRdf * sizeof(real), str);
		safewrite(gd.histRdf12, cmd.sizeHistRdf * sizeof(real), str);
		safewrite(gd.histRdf22, cmd.sizeHistRdf * sizeof(real), str);
	}

    safewrite(&cmd.computePressAxes, sizeof(bool), str);
    safewrite(&cmd.stepAvgPressAxes, sizeof(int), str);
    safewrite(&cmd.stepPressAxes, sizeof(int), str);
    safewrite(&cmd.sizeHistPressAxes, sizeof(int), str);

	if (cmd.computePressAxes) {
		safewrite(&gd.countPressAxes, sizeof(int), str);
		safewrite(gd.histPressX, cmd.sizeHistPressAxes * sizeof(real), str);
		safewrite(gd.histPressY, cmd.sizeHistPressAxes * sizeof(real), str);
#if (NDIM==3)
		safewrite(gd.histPressZ, cmd.sizeHistPressAxes * sizeof(real), str);
#endif
	}

    safewrite(&cmd.computeChemPot, sizeof(bool), str);
    safewrite(&cmd.stepAvgChemPot, sizeof(int), str);
    safewrite(&cmd.stepChemPot, sizeof(int), str);
    safewrite(&cmd.sizeHistChemPot, sizeof(int), str);
    safewrite(&cmd.numTestBodies, sizeof(int), str);

	if (cmd.computeChemPot) {
		safewrite(&gd.countChemPot, sizeof(int), str);
		safewrite(gd.histChemPotX, cmd.sizeHistChemPot * sizeof(real), str);
		safewrite(gd.histChemPotY, cmd.sizeHistChemPot * sizeof(real), str);
#if (NDIM==3)
		safewrite(gd.histChemPotZ, cmd.sizeHistChemPot * sizeof(real), str);
#endif
	}

    safewrite(&cmd.computeTransport, sizeof(bool), str);
    safewrite(&cmd.stepAvgAcf, sizeof(int), str);
    safewrite(&cmd.stepAcf, sizeof(int), str);
    safewrite(&cmd.nBuffAcf, sizeof(int), str);
    safewrite(&cmd.nValAcf, sizeof(int), str);
/*
	if (cmd.computeTransport) {
		safewrite(&gd.countAcfAv, sizeof(int), str);
		safewrite(gd.avAcfVel, cmd.nValAcf * sizeof(real), str);
		safewrite(gd.avAcfTherm, cmd.nValAcf * sizeof(real), str);
		safewrite(gd.avAcfVisc, cmd.nValAcf * sizeof(real), str);

		safewrite(tBuf, cmd.nBuffAcf * sizeof(TBuf), str);
		for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
			safewrite(tBuf[nb].acfVel, cmd.nValAcf * sizeof(real), str);
			safewrite(tBuf[nb].orgVel, gd.nbody * sizeof(vector), str);
			safewrite(tBuf[nb].acfTherm, cmd.nValAcf * sizeof(real), str);
			safewrite(tBuf[nb].acfVisc, cmd.nValAcf * sizeof(real), str);
		}

	}
*/

    safewrite(&cmd.computeSTCorr, sizeof(bool), str);
    safewrite(&cmd.stepAvgCorr, sizeof(int), str);
    safewrite(&cmd.stepCorr, sizeof(int), str);
    safewrite(&cmd.nBuffCorr, sizeof(int), str);
    safewrite(&cmd.nValCorr, sizeof(int), str);
    safewrite(&cmd.nFunCorr, sizeof(int), str);

    safewrite(&cmd.lattCorr_kx, sizeof(real), str);
    safewrite(&cmd.lattCorr_ky, sizeof(real), str);
#ifdef THREEDIM
    safewrite(&cmd.lattCorr_kz, sizeof(real), str);
#endif

    nchars = strlen(cmd.options) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.options, nchars * sizeof(char), str);
    safewrite(&cmd.tstop, sizeof(real), str);

    nchars = strlen(cmd.dtoutstr) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.dtoutstr, nchars * sizeof(char), str);

    safewrite(&gd.dtout, sizeof(real), str);

    nchars = strlen(cmd.dtoutinfostr) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.dtoutinfostr, nchars * sizeof(char), str);

    safewrite(&gd.dtoutinfo, sizeof(real), str);

    safewrite(&gd.tnow, sizeof(real), str);
    safewrite(&gd.tout, sizeof(real), str);
    safewrite(&gd.toutinfo, sizeof(real), str);
    safewrite(&gd.nstep, sizeof(int), str);

	if (gd.stepAvgStatus)
		gd.nstepNew = gd.nstepNew-gd.nstepOld;	// In order to stepAvg work ...
	gd.nstepOld = 0;
    safewrite(&gd.nstepNew, sizeof(int), str);		// In order to stepAvg work ...
    safewrite(&gd.nstepOld, sizeof(int), str);		// In order to stepAvg work ...

    safewrite(&gdtree.rsize, sizeof(real), str);

	ndim = NDIM;
    safewrite(&ndim, sizeof(int), str);
    safewrite(&cmd.seed, sizeof(int), str);
    safewrite(&gd.cputotforce, sizeof(real), str);
    safewrite(&gd.cputotout, sizeof(real), str);
    safewrite(&gd.cputotal, sizeof(real), str);
    safewrite(&gd.cpuchempot, sizeof(real), str);

    safewrite(&gd.sTotEnergy, sizeof(real), str);
    safewrite(&gd.ssTotEnergy, sizeof(real), str);
    safewrite(&gd.sPotEnergy, sizeof(real), str);
    safewrite(&gd.ssPotEnergy, sizeof(real), str);
    safewrite(&gd.sKinEnergy, sizeof(real), str);
    safewrite(&gd.ssKinEnergy, sizeof(real), str);
    safewrite(&gd.sPressure, sizeof(real), str);
    safewrite(&gd.ssPressure, sizeof(real), str);
    safewrite(&gd.sChemPot, sizeof(real), str);
    safewrite(&gd.ssChemPot, sizeof(real), str);

	safewrite(gd.numUcell, NDIM * sizeof(int), str);

    nchars = strlen(cmd.icfile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.icfile, nchars * sizeof(char), str);
    nchars = strlen(cmd.icfilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.icfilefmt, nchars * sizeof(char), str);
    nchars = strlen(cmd.snapoutfile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.snapoutfile, nchars * sizeof(char), str);
    nchars = strlen(cmd.snapoutfilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.snapoutfilefmt, nchars * sizeof(char), str);
    nchars = strlen(cmd.statefile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.statefile, nchars * sizeof(char), str);

    safewrite(&cmd.stepState, sizeof(int), str);
    
    safewrite(bodytab, gd.nbody * sizeof(body), str);

    fclose(str);
}

void restorestate(string file)
{
    stream str;
    int nchars, ndim, nb;
    string program, version;

	strcpy(gd.model_comment, "start from a restore file");

    str = stropen(file, "r");
    saferead(&nchars, sizeof(int), str);
    program = (string) allocate(nchars * sizeof(char));
    saferead(program, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    version = (string) allocate(nchars * sizeof(char));
    saferead(version, nchars * sizeof(char), str);
    if (! streq(program, getargv0()) ||        
          ! streq(version, getversion()))
        printf("warning: state file may be outdated\n\n");

    saferead(&nchars, sizeof(int), str);
    cmd.forcecalc_method = (string) allocate(nchars * sizeof(char));
    saferead(cmd.forcecalc_method, nchars * sizeof(char), str);

    saferead(&cmd.potType, sizeof(int), str);

    saferead(&cmd.intMethod, sizeof(int), str);

    saferead(&nchars, sizeof(int), str);
    cmd.dtimestr = (string) allocate(nchars * sizeof(char));
    saferead(cmd.dtimestr, nchars * sizeof(char), str);

    saferead(&gd.dtime, sizeof(real), str);
    saferead(&cmd.temperature, sizeof(real), str);
    saferead(&cmd.adjustTemperature, sizeof(bool), str);
    saferead(&cmd.adjustCenterOfMass, sizeof(bool), str);
    saferead(&cmd.density, sizeof(real), str);

    saferead(&gd.mass1, sizeof(real), str);
    saferead(&gd.mass2, sizeof(real), str);
    saferead(&gd.Lx, sizeof(real), str);
    saferead(&gd.Ly, sizeof(real), str);
#if (NDIM==3)
    saferead(&gd.Lz, sizeof(real), str);
#endif
    saferead(&cmd.eps11, sizeof(real), str);
    saferead(&cmd.eps12, sizeof(real), str);
    saferead(&cmd.eps22, sizeof(real), str);
    saferead(&cmd.sigma11, sizeof(real), str);
    saferead(&cmd.sigma12, sizeof(real), str);
    saferead(&cmd.sigma22, sizeof(real), str);
    saferead(&cmd.Rcut11, sizeof(real), str);
    saferead(&cmd.Rcut12, sizeof(real), str);
    saferead(&cmd.Rcut22, sizeof(real), str);

    saferead(&gd.nbody, sizeof(int), str);
    saferead(&gd.nbody1, sizeof(int), str);
    saferead(&gd.nbody2, sizeof(int), str);
    saferead(&gd.nbc1, sizeof(int), str);
    saferead(&gd.nbc2, sizeof(int), str);

    saferead(&cmd.stepEquil, sizeof(int), str);
    saferead(&cmd.stepAvg, sizeof(int), str);

    saferead(&cmd.computeRhoAxes, sizeof(bool), str);
    saferead(&cmd.stepAvgRhoAxes, sizeof(int), str);
    saferead(&cmd.stepRhoAxes, sizeof(int), str);
    saferead(&cmd.sizeHistRhoAxes, sizeof(int), str);

	if (cmd.computeRhoAxes) {
		saferead(&gd.countRhoAxes, sizeof(int), str);
		gd.histRhoX = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoY = AllocVecR(cmd.sizeHistRhoAxes);

		gd.histRhoX1 = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoY1 = AllocVecR(cmd.sizeHistRhoAxes);

		gd.histRhoX2 = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoY2 = AllocVecR(cmd.sizeHistRhoAxes);

		fprintf(gd.outlog,"stepAvgRhoAxes=%d\n",cmd.stepAvgRhoAxes);
		fprintf(gd.outlog,"stepRhoAxes=%d\n",cmd.stepRhoAxes);
		fprintf(gd.outlog,"sizeHistRhoAxes=%d\n",cmd.sizeHistRhoAxes);

		saferead(gd.histRhoX, cmd.sizeHistRhoAxes * sizeof(real), str);
		saferead(gd.histRhoY, cmd.sizeHistRhoAxes * sizeof(real), str);
		saferead(gd.histRhoX1, cmd.sizeHistRhoAxes * sizeof(real), str);
		saferead(gd.histRhoY1, cmd.sizeHistRhoAxes * sizeof(real), str);
		saferead(gd.histRhoX2, cmd.sizeHistRhoAxes * sizeof(real), str);
		saferead(gd.histRhoY2, cmd.sizeHistRhoAxes * sizeof(real), str);
#if (NDIM==3)
		gd.histRhoZ = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoZ1 = AllocVecR(cmd.sizeHistRhoAxes);
		gd.histRhoZ2 = AllocVecR(cmd.sizeHistRhoAxes);

		saferead(gd.histRhoZ, cmd.sizeHistRhoAxes * sizeof(real), str);
		saferead(gd.histRhoZ1, cmd.sizeHistRhoAxes * sizeof(real), str);
		saferead(gd.histRhoZ2, cmd.sizeHistRhoAxes * sizeof(real), str);
#endif
	}

    saferead(&cmd.computeNFrecAxes, sizeof(bool), str);
    saferead(&cmd.stepAvgNFrecAxes, sizeof(int), str);
    saferead(&cmd.stepNFrecAxes, sizeof(int), str);
    saferead(&cmd.sizeHistNFrecAxes, sizeof(int), str);

	if (cmd.computeNFrecAxes) {
		saferead(&gd.countNFrecAxes, sizeof(int), str);
		gd.histNFrecX = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecY = AllocVecI(cmd.sizeHistNFrecAxes);

		gd.histNFrecX1 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecY1 = AllocVecI(cmd.sizeHistNFrecAxes);

		gd.histNFrecX2 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecY2 = AllocVecI(cmd.sizeHistNFrecAxes);

		gd.histNFrecXD = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecYD = AllocVecI(cmd.sizeHistNFrecAxes);

		fprintf(gd.outlog,"stepAvgNFrecAxes=%d\n",cmd.stepAvgNFrecAxes);
		fprintf(gd.outlog,"stepNFrecAxes=%d\n",cmd.stepNFrecAxes);
		fprintf(gd.outlog,"sizeHistNFrecAxes=%d\n",cmd.sizeHistNFrecAxes);

		saferead(gd.histNFrecX, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecY, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecX1, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecY1, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecX2, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecY2, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecXD, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecYD, cmd.sizeHistNFrecAxes * sizeof(real), str);
#if (NDIM==3)
		gd.histNFrecZ = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecZ1 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecZ2 = AllocVecI(cmd.sizeHistNFrecAxes);
		gd.histNFrecZD = AllocVecI(cmd.sizeHistNFrecAxes);

		saferead(gd.histNFrecZ, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecZ1, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecZ2, cmd.sizeHistNFrecAxes * sizeof(real), str);
		saferead(gd.histNFrecZD, cmd.sizeHistNFrecAxes * sizeof(real), str);
#endif
	}

    saferead(&gd.stateVel, sizeof(bool), str);
    saferead(&cmd.computeVelDist, sizeof(bool), str);
    saferead(&cmd.stepAvgVel, sizeof(int), str);
    saferead(&cmd.stepVel, sizeof(int), str);
    saferead(&cmd.sizeHistVel, sizeof(int), str);
    saferead(&cmd.rangeVel, sizeof(real), str);

	if (cmd.computeVelDist && gd.stateVel) {
		saferead(&gd.countVel, sizeof(int), str);
		gd.histVel = AllocVecR(cmd.sizeHistVel);

		fprintf(gd.outlog,"stepAvgVel=%d\n",cmd.stepAvgVel);
		fprintf(gd.outlog,"stepVel=%d\n",cmd.stepVel);
		fprintf(gd.outlog,"sizeHistVel=%d\n",cmd.sizeHistVel);
		fprintf(gd.outlog,"rangeVel=%d\n",cmd.rangeVel);

		saferead(gd.histVel, cmd.sizeHistVel * sizeof(real), str);
	}


    // Una vez fijos rangeRdf y sizeHistRdf no se podran modificar ...
    // Falta tener algo parecido para rangeRdf y sizeHistRdf
    // rangeRdf si funciona ... pero hay que dejar pasar la escritura
    // del archivo rdf.dat varias veces para que se estabilize el nuevo
    // valor de rangeRdf ...
    saferead(&gd.stateRdf, sizeof(bool), str);
    saferead(&cmd.computeRdf, sizeof(bool), str);
    saferead(&cmd.stepAvgRdf, sizeof(int), str);
    saferead(&cmd.stepRdf, sizeof(int), str);
    saferead(&cmd.sizeHistRdf, sizeof(int), str);
    saferead(&cmd.rangeRdf, sizeof(real), str);

	if (cmd.computeRdf && gd.stateRdf) {
		saferead(&gd.countRdf, sizeof(int), str);
		gd.histRdf = AllocVecR(cmd.sizeHistRdf);
		gd.histRdf11 = AllocVecR(cmd.sizeHistRdf);
		gd.histRdf12 = AllocVecR(cmd.sizeHistRdf);
		gd.histRdf22 = AllocVecR(cmd.sizeHistRdf);

		fprintf(gd.outlog,"stepAvgRdf=%d\n",cmd.stepAvgRdf);
		fprintf(gd.outlog,"stepRdf=%d\n",cmd.stepRdf);
		fprintf(gd.outlog,"sizeHistRdf=%d\n",cmd.sizeHistRdf);
		fprintf(gd.outlog,"rangeRdf=%d\n",cmd.rangeRdf);

		saferead(gd.histRdf, cmd.sizeHistRdf * sizeof(real), str);
		saferead(gd.histRdf11, cmd.sizeHistRdf * sizeof(real), str);
		saferead(gd.histRdf12, cmd.sizeHistRdf * sizeof(real), str);
		saferead(gd.histRdf22, cmd.sizeHistRdf * sizeof(real), str);
	}

    saferead(&cmd.computePressAxes, sizeof(bool), str);
    saferead(&cmd.stepAvgPressAxes, sizeof(int), str);
    saferead(&cmd.stepPressAxes, sizeof(int), str);
    saferead(&cmd.sizeHistPressAxes, sizeof(int), str);

	if (cmd.computePressAxes) {
		saferead(&gd.countPressAxes, sizeof(int), str);
		gd.histPressX = AllocVecR(cmd.sizeHistPressAxes);
		gd.histPressY = AllocVecR(cmd.sizeHistPressAxes);

		fprintf(gd.outlog,"stepAvgPressAxes=%d\n",cmd.stepAvgPressAxes);
		fprintf(gd.outlog,"stepPressAxes=%d\n",cmd.stepPressAxes);
		fprintf(gd.outlog,"sizeHistPressAxes=%d\n",cmd.sizeHistPressAxes);

		saferead(gd.histPressX, cmd.sizeHistPressAxes * sizeof(real), str);
		saferead(gd.histPressY, cmd.sizeHistPressAxes * sizeof(real), str);
#if (NDIM==3)
		gd.histPressZ = AllocVecR(cmd.sizeHistPressAxes);

		saferead(gd.histPressZ, cmd.sizeHistPressAxes * sizeof(real), str);
#endif
	}

    saferead(&cmd.computeChemPot, sizeof(bool), str);
    saferead(&cmd.stepAvgChemPot, sizeof(int), str);
    saferead(&cmd.stepChemPot, sizeof(int), str);
    saferead(&cmd.sizeHistChemPot, sizeof(int), str);
    saferead(&cmd.numTestBodies, sizeof(int), str);

	if (cmd.computeChemPot) {
		saferead(&gd.countChemPot, sizeof(int), str);
		gd.histChemPotX = AllocVecR(cmd.sizeHistChemPot);
		gd.histChemPotY = AllocVecR(cmd.sizeHistChemPot);

		fprintf(gd.outlog,"stepAvgChemPot=%d\n",cmd.stepAvgChemPot);
		fprintf(gd.outlog,"stepChemPot=%d\n",cmd.stepChemPot);
		fprintf(gd.outlog,"sizeHistChemPot=%d\n",cmd.sizeHistChemPot);

		saferead(gd.histChemPotX, cmd.sizeHistChemPot * sizeof(real), str);
		saferead(gd.histChemPotY, cmd.sizeHistChemPot * sizeof(real), str);
#if (NDIM==3)
		gd.histChemPotZ = AllocVecR(cmd.sizeHistChemPot);

		saferead(gd.histChemPotZ, cmd.sizeHistChemPot * sizeof(real), str);
#endif
	}

    saferead(&cmd.computeTransport, sizeof(bool), str);
    saferead(&cmd.stepAvgAcf, sizeof(int), str);
    saferead(&cmd.stepAcf, sizeof(int), str);
    saferead(&cmd.nBuffAcf, sizeof(int), str);
    saferead(&cmd.nValAcf, sizeof(int), str);
/*
	if (cmd.computeTransport) {
		saferead(&gd.countAcfAv, sizeof(int), str);

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

		saferead(gd.avAcfVel, cmd.nValAcf * sizeof(real), str);
		saferead(gd.avAcfTherm, cmd.nValAcf * sizeof(real), str);
		saferead(gd.avAcfVisc, cmd.nValAcf * sizeof(real), str);

		saferead(tBuf, cmd.nBuffAcf * sizeof(TBuf), str);
		for (nb = 0; nb < cmd.nBuffAcf; nb ++) {
			saferead(tBuf[nb].acfVel, cmd.nValAcf * sizeof(real), str);
			saferead(tBuf[nb].orgVel, gd.nbody * sizeof(vector), str);
			saferead(tBuf[nb].acfTherm, cmd.nValAcf * sizeof(real), str);
			saferead(tBuf[nb].acfVisc, cmd.nValAcf * sizeof(real), str);
		}

	}
*/

    saferead(&cmd.computeSTCorr, sizeof(bool), str);
    saferead(&cmd.stepAvgCorr, sizeof(int), str);
    saferead(&cmd.stepCorr, sizeof(int), str);
    saferead(&cmd.nBuffCorr, sizeof(int), str);
    saferead(&cmd.nValCorr, sizeof(int), str);
    saferead(&cmd.nFunCorr, sizeof(int), str);

    saferead(&cmd.lattCorr_kx, sizeof(real), str);
    saferead(&cmd.lattCorr_ky, sizeof(real), str);
#ifdef THREEDIM
    saferead(&cmd.lattCorr_kz, sizeof(real), str);
#endif

    saferead(&nchars, sizeof(int), str);
    cmd.options = (string) allocate(nchars * sizeof(char));
    saferead(cmd.options, nchars * sizeof(char), str);
    saferead(&cmd.tstop, sizeof(real), str);

    saferead(&nchars, sizeof(int), str);
    cmd.dtoutstr = (string) allocate(nchars * sizeof(char));
    saferead(cmd.dtoutstr, nchars * sizeof(char), str);

    saferead(&gd.dtout, sizeof(real), str);

    saferead(&nchars, sizeof(int), str);
    cmd.dtoutinfostr = (string) allocate(nchars * sizeof(char));
    saferead(cmd.dtoutinfostr, nchars * sizeof(char), str);
    saferead(&gd.dtoutinfo, sizeof(real), str);
    saferead(&gd.tnow, sizeof(real), str);
    saferead(&gd.tout, sizeof(real), str);
    saferead(&gd.toutinfo, sizeof(real), str);
    saferead(&gd.nstep, sizeof(int), str);

    saferead(&gd.nstepNew, sizeof(int), str);		// In order to stepAvg work ...
    saferead(&gd.nstepOld, sizeof(int), str);		// In order to stepAvg work ...

    saferead(&gdtree.rsize, sizeof(real), str);

    saferead(&ndim, sizeof(int), str);
	if (ndim != NDIM)
		error("\nrestorestate : ndim = %d; expected %d\n",ndim,NDIM);
    saferead(&cmd.seed, sizeof(int), str);
    saferead(&gd.cputotforce, sizeof(real), str);
    saferead(&gd.cputotout, sizeof(real), str);
    saferead(&gd.cputotal, sizeof(real), str);
    saferead(&gd.cpuchempot, sizeof(real), str);

    saferead(&gd.sTotEnergy, sizeof(real), str);
    saferead(&gd.ssTotEnergy, sizeof(real), str);
    saferead(&gd.sPotEnergy, sizeof(real), str);
    saferead(&gd.ssPotEnergy, sizeof(real), str);
    saferead(&gd.sKinEnergy, sizeof(real), str);
    saferead(&gd.ssKinEnergy, sizeof(real), str);
    saferead(&gd.sPressure, sizeof(real), str);
    saferead(&gd.ssPressure, sizeof(real), str);
    saferead(&gd.sChemPot, sizeof(real), str);
    saferead(&gd.ssChemPot, sizeof(real), str);

	saferead(gd.numUcell, NDIM * sizeof(int), str);

    saferead(&nchars, sizeof(int), str);
    cmd.icfile = (string) allocate(nchars * sizeof(char));
    saferead(cmd.icfile, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    cmd.icfilefmt = (string) allocate(nchars * sizeof(char));
    saferead(cmd.icfilefmt, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    cmd.snapoutfile = (string) allocate(nchars * sizeof(char));
    saferead(cmd.snapoutfile, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    cmd.snapoutfilefmt = (string) allocate(nchars * sizeof(char));
    saferead(cmd.snapoutfilefmt, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    cmd.statefile = (string) allocate(nchars * sizeof(char));
    saferead(cmd.statefile, nchars * sizeof(char), str);

    saferead(&cmd.stepState, sizeof(int), str);
    
    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));
    saferead(bodytab, gd.nbody * sizeof(body), str);

    fclose(str);
}

local void Global_to_Header(void)
{
	hdr.nbody		= gd.nbody;
	hdr.nbody1		= gd.nbody1;
	hdr.nbody2		= gd.nbody2;
	hdr.tnow			= gd.tnow;
	hdr.temperature	= cmd.temperature;
	hdr.density		= cmd.density;
	hdr.mass1		= gd.mass1;
	hdr.mass2		= gd.mass2;
	hdr.Lx			= gd.Lx;
	hdr.Ly			= gd.Ly;
#ifdef THREEDIM
	hdr.Lz			= gd.Lz;
#endif
	hdr.eps11		= cmd.eps11;
	hdr.eps12		= cmd.eps12;
	hdr.eps22		= cmd.eps22;
	hdr.sigma11		= cmd.sigma11;
	hdr.sigma12		= cmd.sigma12;
	hdr.sigma22		= cmd.sigma22;
	hdr.Rcut11		= cmd.Rcut11;
	hdr.Rcut12		= cmd.Rcut12;
	hdr.Rcut22		= cmd.Rcut22;
}

#define stopfilename	"stop"
#define stopsavestate	"stop-state"

void checkstop(void)
{
    char   buf[200];
	FILE *fd;
	double cpudt;

	if ((fd=fopen(stopfilename,"r"))) {
		fclose(fd);
		gd.stopflag = 1;
		sprintf(buf,"rm -f %s", stopfilename);
		system(buf);
        printf("\nsaving a stop run state...\n\n");
		cpudt = cputime()-gd.cpuinit;
		gd.cputotal += cpudt;
		Global_to_Header();
        savestate(stopsavestate);
		gd.cputotal -= cpudt;
	}
}

#undef stopfilename
#undef stopsavestate

void EndRun(void)
{
    char   buf[200];
	FILE *fd;
	double cpudt;

    if (! strnull(cmd.statefile)) {
        printf("\nsaving final run state...\n\n");
		cpudt = cputime()-gd.cpuinit;
		gd.cputotal += cpudt;
		Global_to_Header();
        savestate(savestatetmp);
		gd.cputotal -= cpudt;
        sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
        system(buf);
    }

	fclose(gd.outlog);
	gd.cputotal += cputime()-gd.cpuinit;
	printf("\nFinal CPU time : %lf\n\n",gd.cputotal);
}

#undef savestatetmp
