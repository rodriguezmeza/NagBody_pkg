/* ==============================================================================
!	MODULE: md_lj_tree_io.c														!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date:	January, 2005.												!
!	Purpose: Routines to drive input and output data							!
!	Language: C																	!
!	Use:																		!
!	Routines and functions:														!
!	External modules, routines and headers:	stdinc.h, mathfns.h, vectmath		!
!						vectmath.h, getparam.h									!
!						types.h, stat.h, inout.h								!
!	Comments and notes:															!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.							!
!		e-mail: marioalberto.rodriguez@inin.gob.mx
!		http://www.astro.inin.mx/mar											!
!																				!
!	Major revisions:															!
!	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.							!
!===============================================================================
!	Legal matters:																!
!	The author does not warrant that the program and routines it contains		!
!	listed below are free from error or suitable for particular applications,	!
!	and he disclaims all liability from any consequences arising from their		!
!	use.																		!
!==============================================================================*/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "inout.h"
#include "constant.h"
#include "global_defs.h"
#include "proto_defs.h"

#include <sys/types.h>
#include <sys/stat.h>

/* #include <strings.h> */							// For unix
#include "strings.h"		// For Visual C

local void AdjustTemp(void);
local void EvalProps(void);
local void AccumProps(int);
local void PrintSummary(FILE *);
local void EvalVelDist(void);
local void PrintVelDist(FILE *);
local void EvalRdf(void);
local void PrintRdf(FILE *);
local void PrintSnap(FILE *);

local void outputdata(void);
local void diagnostics(void);

local void outputpvdata(void);
local void outputbindata(void);
local void inputdata_ascii(void);
local void inputdata_bin(void);
local void inputdata_pv(void);

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

local void outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;

#define IC_SNAP_FMT	0
#define IC_NULL_FMT	1
#define IC_PV_FMT		2
#define IC_SNAP_FMT_BIN	3

void inputdata(void)
{
	switch(icfilefmt_int) {
		case IC_SNAP_FMT: 
			inputdata_ascii(); break;
		case IC_NULL_FMT: 
			inputdata_ascii(); break;
		case IC_PV_FMT:
			inputdata_pv(); break;
		case IC_SNAP_FMT_BIN: 
			inputdata_bin(); break;
		default: 
			printf("\n\tinput: Unknown input format...");
			printf("\n\tinput in default snap (ascii) format...\n"); 
			inputdata_ascii(); break;
	}
}

#undef IC_SNAP_FMT
#undef IC_NULL_FMT
#undef IC_PV_FMT
#undef IC_SNAP_FMT_BIN

local void inputdata_ascii(void)
{
    stream instr;
    int ndim;
    bodyptr p;

    model_comment = "ASCII input file";

    instr = stropen(icfile, "r");               
    in_int(instr, &nbody);                      
    if (nbody < 1)
        error("inputdata: nbody = %d is absurd\n", nbody);
    in_int(instr, &ndim);                       
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real(instr, &tnow);                      
    bodytab = (bodyptr) allocate(nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+nbody)
        in_real(instr, &Mass(p));               
	DO_BODY(p, bodytab, bodytab+nbody)
        in_vector(instr, Pos(p));               
	DO_BODY(p, bodytab, bodytab+nbody)
        in_vector(instr, Vel(p));               
    fclose(instr);
    if (scanopt(options, "reset-time"))         
        tnow = 0.0;                             
	DO_BODY(p, bodytab, bodytab+nbody) {
        Type(p) = BODY;
		Id(p) = p-bodytab+1;
	}
}

local void inputdata_pv(void)
{
    stream instr;
    int ndim;
    bodyptr p;
	char gato[1], firstline[20];

    model_comment = "Column form input file";

    instr = stropen(icfile, "r");
	fgets(firstline,200,instr);
	fscanf(instr,"%1s",gato);
    in_int(instr, &nbody);                      
    if (nbody < 1)
        error("inputdata: nbody = %d is absurd\n", nbody);
    in_int(instr, &ndim);                       
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real(instr, &tnow);
    bodytab = (bodyptr) allocate(nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+nbody) {
		in_int(instr, &Id(p));
        in_real(instr, &Mass(p));               
        in_vector(instr, Pos(p));               
        in_vector(instr, Vel(p));
	}

    fclose(instr);

    if (scanopt(options, "reset-time"))
        tnow = 0.0;
	DO_BODY(p, bodytab, bodytab+nbody)
        Type(p) = BODY;
}

local void inputdata_bin(void)
{
    stream instr;
    int ndim;
    bodyptr p;

    model_comment = "Binary input file";

    instr = stropen(icfile, "r");               
    in_int_bin(instr, &nbody);                      
    if (nbody < 1)
        error("inputdata: nbody = %d is absurd\n", nbody);
    in_int_bin(instr, &ndim);                       
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real_bin(instr, &tnow);                      
    bodytab = (bodyptr) allocate(nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+nbody)
        in_real_bin(instr, &Mass(p));               
	DO_BODY(p, bodytab, bodytab+nbody)
        in_vector_bin(instr, Pos(p));               
	DO_BODY(p, bodytab, bodytab+nbody)
        in_vector_bin(instr, Vel(p));               
    fclose(instr);                              
    if (scanopt(options, "reset-time"))         
        tnow = 0.0;                             
	DO_BODY(p, bodytab, bodytab+nbody) {
        Type(p) = BODY;                         
		Id(p) = p-bodytab+1;
	}
}

#define outfile_thermo		"thermo.dat"
#define savestatetmp "savestate-tmp"

void startoutput(void)
{
    stream outstr_thermo;
	char buf[200];

	if (strnull(restorefile)) {
		sprintf(buf,"%s",outfile_thermo);
		if(!(outstr_thermo=fopen(buf,mode)))
		{
			error("can't open file `%s`\n",buf);
		}
		fprintf(outstr_thermo,"%1s%4s%8s%8s%8s%8s%7s%8s%7s%9s%7s%9s\n",
          "#","nstep","time","vSum","TotE","ssTotE","kE","sskE","p","ssp","T","rho");
		fprintf(outstr_thermo,"%1s%4s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n",
          "#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>","<10>","<11>");
		fclose(outstr_thermo);
	}

    outfilefmt_string_to_int(snapoutfilefmt, &outfilefmt_int);

    printf("\n%s\n%s: %s\n\t %s\n", headline0, headline1, headline2, headline3);
    printf("  \t -- %s --\n", model_comment);
	printf("\n  \t running ... %s force calculation method ...\n",forcecalc_method);

    printf("\n%12s%9s%8s%7s", "temperature", "density","nbody", "dtime");
    printf("%9s%10s%8s\n", "dtout", "dtoutinfo", "tstop");
    printf("%8.2f     %8.4f%8d%8.5f", temperature, density, nbody, dtime);
    printf("%8.2f%8.2f %9.4f\n", dtout, dtoutinfo, tstop);

    if (! strnull(options))                     
        printf("\n\toptions: %s\n", options);

	if (! strnull(statefile)) {
		savestate(savestatetmp);
		sprintf(buf,"cp %s %s",savestatetmp,statefile);
		printf("system: %s",buf);
		system(buf);
	}
	
	if (strnull(restorefile))
		AccumProps(0);

	countVel = 0;
    histVel = AllocVecR(sizeHistVel);
	
	fprintf(outlog,"\n\nstepEquil=%d\n",stepEquil);
	fprintf(outlog,"stepAvgVel=%d\n",stepAvgVel);
	fprintf(outlog,"stepVel=%d\n",stepVel);
	fprintf(outlog,"sizeHistVel=%d\n",sizeHistVel);
	fprintf(outlog,"rangeVel=%g\n",rangeVel);

	countRdf = 0;
	histRdf = AllocVecR(sizeHistRdf);
	
	fprintf(outlog,"stepAvgRdf=%d\n",stepAvgRdf);
	fprintf(outlog,"stepRdf=%d\n",stepRdf);
	fprintf(outlog,"sizeHistRdf=%d\n",sizeHistRdf);
	fprintf(outlog,"rangeRdf=%g\n",rangeRdf);
	fprintf(outlog,"stepSnapInit=%d\n",stepSnapInit);
	fprintf(outlog,"stepSnap=%d\n",stepSnap);
	fflush(outlog);
}

#define SNAP_FMT	0
#define NULL_FMT	1
#define PV_FMT		2
#define SNAP_FMT_BIN	3

void output(void)
{
    real cmabs, amabs, teff;
    char   buf[200];
	real drmin, vmax, tmp, amax, rri, rri3;
	bodyptr p, q;
    double cpustart, cpudt;
	int flaginfo=0;

    cpustart = cputime();                       

    teff = tnow + dtime/8.0;  

	EvalProps(); AccumProps(1);

    if (teff >= toutinfo) { 

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
	cputotforce += cpuforce;
    printf("%8d %8.3f %8.5f  %9.5f %9d %8.3f %9.3f\n", nstep,
           tnow, cmabs, amabs, nbbcalc, cpuforce, cputotforce);

    if (! scanopt(forcecalc_method, "direct")  &&
		! scanopt(forcecalc_method, "direct2") &&
		! scanopt(forcecalc_method, "cells")	) {
    printf("\n\t%8s%8s%8s%8s%10s%10s%8s\n",
           "rsize", "tdepth", "ftree",
           "actmax", "nbbtot", "nbctot", "CPUtree");
    printf("\t%8.1f%8d%8.3f%8d%10d%10d%8.3f\n",
           rsize, tdepth, (nbody + ncell - 1) / ((real) ncell),
           actmax, nbbcalc, nbccalc, cputree);
	}

    printf("\n%9s %6s %9s %9s %10s %8s %8s\n", "Averages:",
           "Etot/N", "kinE/N", "PotE/N", "pressure", "vSum", "vvSum");
    printf("       %9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
           totEnergy,kinEnergy,potEnergy,pressure,vSum,vvSum);

#if defined(DIAGNOSTICS)
    printf("\n %8s %8s %8s %10s\n",
           "Etot", "kinE", "PotE", "virSum");
    printf("%9.4f %9.4f %9.4f %9.4f\n",
           etot[0],etot[1],etot[2],virSum);

	drmin = Box[0]; vmax = 0.; 
	DO_BODY(p, bodytab, bodytab+nbody) {
		ABSV(tmp,Vel(p)); vmax = MAX(vmax, tmp);
	}
	DO_BODY(p, bodytab, bodytab+nbody-1)
		DO_BODY(q, p+1, bodytab+nbody) {
			DISTV(tmp, Pos(p), Pos(q));
			drmin = MIN(drmin, tmp);
		}
	rri=ssq/drmin; rri3=rri*rri*rri;
	amax = fa*rri3*(rri3-0.5)*rri*drmin;
    printf("\n %8s %8s %8s %10s\n",
           "drmin", "vmax", "drmin/vmax", "amax");
    printf("%9.4f %9.4f %9.4f %9.4f\n",
           drmin,vmax,drmin/vmax,amax);
#endif
	toutinfo += dtoutinfo;
	}

    teff = tnow + dtime/8.0;  
//	EvalProps(); AccumProps(1);

	if (nstep >= stepEquil && (nstep - stepEquil)%stepVel == 0)
		EvalVelDist();
	if (nstep >= stepEquil && (nstep - stepEquil)%stepRdf == 0)
		EvalRdf();
	if (nstep >= stepSnapInit && (nstep - stepSnapInit)%stepSnap == 0)
		PrintSnap(stdout);
	if (nstep%stepAvg == 0) {
		AccumProps(2);
//		AdjustTemp();						// Re-escaling is in "timestep"	
		PrintSummary(stdout); AccumProps(0);
	}

    if (! strnull(snapoutfile) && teff >= tout) { 
        switch(outfilefmt_int) {
            case SNAP_FMT: 
                printf("\n\tsnap format output"); outputdata(); break;
            case NULL_FMT: 
                printf("\n\tsnap format output"); outputdata(); break;
            case PV_FMT: 
                printf("\n\tpv format output"); outputpvdata(); break;
            case SNAP_FMT_BIN: 
                printf("\n\tsnap-bin format output"); outputbindata(); break;
            default: 
                printf("\n\toutput: Unknown output format...");
                printf("\n\tprinting in default format..."); outputdata(); break;
        }
		tout += dtout;
    }

	if (nstep%stepSnap == 0) {				// stepSnap = dtout/dtime
        if (! strnull(statefile)) {
			cpudt = cputime()-cpuinit;
			cputotal += cpudt;
            savestate(savestatetmp);
			cputotal -= cpudt;
            sprintf(buf,"cp %s %s",savestatetmp,statefile);
            printf("system: %s",buf);
            system(buf);
        }
	}
	cputotout += cputime()-cpustart;
	if (flaginfo)
		printf("Accum. output cpu time = %g\n",cputotout);
}

#undef SNAP_FMT
#undef NULL_FMT
#undef PV_FMT
#undef SNAP_FMT_BIN


local void outputdata(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf, snapoutfile, nstep);
    outstr = stropen(namebuf, "w!");         
    out_int(outstr, nbody);
    out_int(outstr, NDIM);                      
    out_real(outstr, tnow);                     
	DO_BODY(p, bodytab, bodytab+nbody)
        out_real(outstr, Mass(p));              
	DO_BODY(p, bodytab, bodytab+nbody)
        out_vector(outstr, Pos(p));             
	DO_BODY(p, bodytab, bodytab+nbody)
        out_vector(outstr, Vel(p));             
    if (scanopt(options, "out-phi"))           
		DO_BODY(p, bodytab, bodytab+nbody)
            out_real(outstr, Phi(p));           
    if (scanopt(options, "out-acc"))            
		DO_BODY(p, bodytab, bodytab+nbody)
            out_vector(outstr, Acc(p));         
    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
}

local void outputbindata(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf, snapoutfile, nstep);
    outstr = stropen(namebuf, "w!");         
    out_int_bin(outstr, nbody);
    out_int_bin(outstr, NDIM);                      
    out_real_bin(outstr, tnow);
	DO_BODY(p, bodytab, bodytab+nbody)
        out_real_bin(outstr, Mass(p));              
	DO_BODY(p, bodytab, bodytab+nbody)
        out_vector_bin(outstr, Pos(p));             
	DO_BODY(p, bodytab, bodytab+nbody)
        out_vector_bin(outstr, Vel(p));             
    if (scanopt(options, "out-phi"))           
		DO_BODY(p, bodytab, bodytab+nbody)
            out_real_bin(outstr, Phi(p));           
    if (scanopt(options, "out-acc"))            
		DO_BODY(p, bodytab, bodytab+nbody)
            out_vector_bin(outstr, Acc(p));         
    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
}

local void outputpvdata(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf, snapoutfile, nstep);
    outstr = stropen(namebuf, "w!");
	fprintf(outstr,"#   nbody NDIM time\n# %d %d %lf\n",nbody,NDIM,tnow);
	DO_BODY(p, bodytab, bodytab+nbody) {
		out_int_mar(outstr, Id(p));
		out_real_mar(outstr, Mass(p));
        out_vector_mar(outstr, Pos(p));             
        out_vector_mar(outstr, Vel(p));
        if (scanopt(options, "out-phi"))
            out_real_mar(outstr, Phi(p));
        if (scanopt(options, "out-acc"))
            out_vector_mar(outstr, Acc(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, tnow);
}

local void diagnostics(void)			// se mantiene para chequeo cruzado
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
	DO_BODY(p, bodytab, bodytab+nbody) {
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

local void AdjustTemp(void)
{
	real vFac;
    bodyptr p;

	vFac = vMag/rsqrt( 2.0*sKinEnergy/mass );
	DO_BODY(p, bodytab, bodytab+nbody)
        MULVS(Vel(p), Vel(p), vFac);			// Re-scaling vels.
}

local void EvalProps(void) 
{
    bodyptr p;
	real v, vv, dt;
	int k;
	
    dt = dtime;                            
	vSum = vvSum = 0.0;
	DO_BODY(p, bodytab, bodytab+nbody) {
		vv = 0.0;
		DO_COORD(k) {							// No se sincronizan vels
			v = Vel(p)[k];						// - 0.5*Acc(p)[k]*dt;
			vSum += v; vv += rsqr(v);
		}
		vvSum += Mass(p)*vv;
	}
	kinEnergy = 0.5*vvSum/((real) (nbody));
	potEnergy = uSum/((real) (nbody));
	totEnergy = kinEnergy + potEnergy;
	pressure = density * (vvSum + virSum) / ((real) (nbody*NDIM));
/*
    printf("\n%9s %6s %9s %9s %10s %8s %8s\n", "Averages:",
           "Etot/N", "kinE/N", "PotE/N", "pressure", "vSum", "vvSum");
    printf("       %9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
           totEnergy,kinEnergy,potEnergy,pressure,vSum,vvSum);
*/
}

local void AccumProps(int icode) 
{
	if (icode == 0) {
		sTotEnergy = ssTotEnergy = 0.0;
		sKinEnergy = ssKinEnergy = 0.0;
		sPressure = ssPressure = 0.0;
	} if (icode == 1) {
		sTotEnergy += totEnergy;
		ssTotEnergy += rsqr(totEnergy);
		sKinEnergy += kinEnergy;
		ssKinEnergy += rsqr(kinEnergy);
		sPressure += pressure;
		ssPressure += rsqr(pressure);
	} if (icode == 2) {
		sTotEnergy = sTotEnergy/stepAvg;
		ssTotEnergy = rsqrt( (ssTotEnergy/stepAvg > rsqr(sTotEnergy) ? 
			ssTotEnergy/stepAvg - rsqr(sTotEnergy) : 0) );
		sKinEnergy = sKinEnergy/stepAvg;
		ssKinEnergy = rsqrt( (ssKinEnergy/stepAvg > rsqr(sKinEnergy) ?
			ssKinEnergy/stepAvg - rsqr(sKinEnergy) : 0) );
		sPressure = sPressure/stepAvg;
		ssPressure = rsqrt( (ssPressure/stepAvg > rsqr(sPressure) ?
			ssPressure/stepAvg - rsqr(sPressure) : 0) );
	}
}

local void PrintSummary(FILE *fp)
{
    stream outstr_thermo;

    outstr_thermo = stropen(outfile_thermo, "a");
	fprintf(outstr_thermo,
		"%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f",
		nstep, tnow, vSum, sTotEnergy, ssTotEnergy,
		sKinEnergy, ssKinEnergy, sPressure, ssPressure, temperature, density);
	fprintf(outstr_thermo,"\n");
    fclose(outstr_thermo);
}

#undef outfile_thermo

local void EvalVelDist(void)
{
    bodyptr p;
	real deltaV, histSum, vv;
	int j, k;
	double cpustart;

	printf("\nEvalVelDist: Entrando ... ");
	cpustart = cputime();                       

	countVel = countVel + 1;
	if (countVel == 1) {
		for (j = 1; j < sizeHistVel; j++) histVel[j] = 0.0;
	}
	deltaV = rangeVel/sizeHistVel;
	DO_BODY(p, bodytab, bodytab+nbody) {
		vv =0.0;
		DO_COORD(k)
			vv += rsqr(Vel(p)[k]);
		j = (int) (rsqrt(vv)/deltaV)+1;
		if (j>sizeHistVel) j=sizeHistVel;
		histVel[j] = histVel[j]+1.0;
	}
	if (countVel == stepAvgVel) {
		histSum = 0.0;
		for (j = 1; j < sizeHistVel; j++)
			histSum= histSum+histVel[j];
		for (j=1; j<=sizeHistVel; j++)
			histVel[j]= histVel[j]/histSum;
		PrintVelDist(stdout);
		countVel = 0.0;
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
	
	for (n=1; n<=sizeHistVel; n++) {
		vBin = (n-0.5)*rangeVel/sizeHistVel;
		fprintf(outstr_vel,"%8.3f %8.3f\n",vBin,histVel[n]);
	}
    fclose(outstr_vel);
	sprintf(buf,"mv %s %s",saveveltmp,savevel);
	printf("\nsystem: %s ...\n",buf);
	system(buf);
}

#undef savevel


local void EvalRdf(void)
{
    bodyptr j1, j2;
	real deltaR, normFac, rr, rrRange, Vol;
	int k, n;
	vector dr;
	double cpustart;

	printf("\nEvalRdf: Entrando ... ");

	cpustart = cputime();                       

	countRdf = countRdf + 1;
	if (countRdf == 1) {
		for (n = 1; n <= sizeHistRdf; n++) histRdf[n] = 0.0;
	}
	rrRange = rsqr(rangeRdf);
	deltaR = rangeRdf/sizeHistRdf;
	DO_BODY(j1, bodytab, bodytab+nbody-1) {
		DO_BODY(j2, j1+1, bodytab+nbody) {
			DO_COORD(k) {
				dr[k] = Pos(j1)[k] - Pos(j2)[k];
				dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
			}
			rr=0.0;
			DO_COORD(k)
				rr += rsqr(dr[k]);
			if (rr < rrRange) {
				n = (int) (rsqrt(rr) / deltaR) + 1;
				histRdf[n] = histRdf[n] + 1.;
			}
		}
	}
	if (countRdf == stepAvgRdf) {
		Vol = 1.0;
		DO_COORD(k)
			Vol = Vol*Box[k];
		if (NDIM == 3)
			normFac = Vol /	(2.0 * PI * rpow(deltaR, 3.0) * nbody * nbody * countRdf);
		else if (NDIM == 2)
			normFac = Vol /	(PI * rpow(deltaR, 2.0) * nbody * nbody * countRdf);
		else error("\n\nWrong NDIM!\n\n");
		
		for (n = 1; n <= sizeHistRdf; n++)
			if (NDIM == 3)
				histRdf[n] = histRdf[n] * normFac / rsqr(n-0.5);
			else if (NDIM == 2)
				histRdf[n] = histRdf[n] * normFac / ((int)n-0.5);
			else error("\n\nWrong NDIM!\n\n");

		PrintRdf(stdout);
		countRdf = 0;
	}
	
	printf("Saliendo : CPU time = %lf\n",cputime()-cpustart); 
}

#define saverdftmp "rdf.tmp"
#define saverdf "rdf.dat"

local void PrintRdf(FILE *fp)
{
	real rBin;
	int n;
    stream outstr_rdf;
    char   buf[200];

    outstr_rdf = stropen(saverdftmp, "w!");
	
	for (n=1; n<=sizeHistRdf; n++) {
		rBin = ((int)n-0.5)*rangeRdf/sizeHistRdf;
		fprintf(outstr_rdf,"%8.4f %8.4f\n",rBin,histRdf[n]);
	}
    fclose(outstr_rdf);
	sprintf(buf,"mv %s %s",saverdftmp,saverdf);
	printf("\nsystem: %s ...\n",buf);
	system(buf);
}

#undef saverdf

#define savesnaptmp "snap.tmp"
#define savesnap "snap.dat"

local void PrintSnap(FILE *fp)
{
    stream outstr_snap;
    bodyptr p;
    char   buf[200];

    outstr_snap = stropen(savesnaptmp, "w!");

	DO_BODY(p, bodytab, bodytab+nbody) {
        out_vector_mar(outstr_snap, Pos(p));             
        out_vector_mar(outstr_snap, Vel(p));
        if (scanopt(options, "out-phi"))
            out_real_mar(outstr_snap, Phi(p));
        if (scanopt(options, "out-acc"))
            out_vector_mar(outstr_snap, Acc(p));
        fprintf(outstr_snap,"\n");
    }
    fclose(outstr_snap);
	sprintf(buf,"mv %s %s",savesnaptmp,savesnap);
	printf("\nsystem: %s\n\n",buf);
	system(buf);
}

#undef savesnap


local void savestate(string pattern)
{
    char namebuf[256];
    stream str;
    int nchars, ndim;

    sprintf(namebuf, pattern, nstep & 1);      
    str = stropen(namebuf, "w!");
    nchars = strlen(getargv0()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getargv0(), nchars * sizeof(char), str);
    nchars = strlen(getversion()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getversion(), nchars * sizeof(char), str);

    nchars = strlen(forcecalc_method) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(forcecalc_method, nchars * sizeof(char), str);

    safewrite(&dtime, sizeof(real), str);
    safewrite(&temperature, sizeof(real), str);
    safewrite(&density, sizeof(real), str);
    safewrite(&stepEquil, sizeof(int), str);

    safewrite(&stepAvgVel, sizeof(int), str);
    safewrite(&stepVel, sizeof(int), str);
    safewrite(&sizeHistVel, sizeof(int), str);
    safewrite(&rangeVel, sizeof(real), str);

    safewrite(&stepAvgRdf, sizeof(int), str);
    safewrite(&stepRdf, sizeof(int), str);
    safewrite(&sizeHistRdf, sizeof(int), str);
    safewrite(&rangeRdf, sizeof(real), str);

    safewrite(&theta, sizeof(real), str);
    safewrite(&usequad, sizeof(bool), str);
    nchars = strlen(options) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(options, nchars * sizeof(char), str);
    safewrite(&tstop, sizeof(real), str);
    safewrite(&dtout, sizeof(real), str);
    safewrite(&dtoutinfo, sizeof(real), str);
    safewrite(&tnow, sizeof(real), str);
    safewrite(&tout, sizeof(real), str);
    safewrite(&nstep, sizeof(int), str);
    safewrite(&rsize, sizeof(real), str);
    safewrite(&nbody, sizeof(int), str);

	ndim = NDIM;
    safewrite(&ndim, sizeof(int), str);
    safewrite(&seed, sizeof(int), str);
    safewrite(&cputotforce, sizeof(real), str);
    safewrite(&cputotout, sizeof(real), str);
    safewrite(&cputotal, sizeof(real), str);

    safewrite(&sTotEnergy, sizeof(real), str);
    safewrite(&ssTotEnergy, sizeof(real), str);
    safewrite(&sKinEnergy, sizeof(real), str);
    safewrite(&ssKinEnergy, sizeof(real), str);
    safewrite(&sPressure, sizeof(real), str);
    safewrite(&ssPressure, sizeof(real), str);

    nchars = strlen(icfile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(icfile, nchars * sizeof(char), str);
    nchars = strlen(icfilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(icfilefmt, nchars * sizeof(char), str);
    nchars = strlen(snapoutfile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(snapoutfile, nchars * sizeof(char), str);
    nchars = strlen(snapoutfilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(snapoutfilefmt, nchars * sizeof(char), str);
    nchars = strlen(statefile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(statefile, nchars * sizeof(char), str);
//    nchars = strlen(restorefile) + 1;
//    safewrite(&nchars, sizeof(int), str);
//    safewrite(restorefile, nchars * sizeof(char), str);

    safewrite(bodytab, nbody * sizeof(body), str);

    fclose(str);
}


void restorestate(string file)
{
    stream str;
    int nchars, ndim;
    string program, version;

    model_comment = "start from a restore file";

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
    forcecalc_method = (string) allocate(nchars * sizeof(char));
    saferead(forcecalc_method, nchars * sizeof(char), str);

    saferead(&dtime, sizeof(real), str);
    saferead(&temperature, sizeof(real), str);
    saferead(&density, sizeof(real), str);
    saferead(&stepEquil, sizeof(int), str);

    saferead(&stepAvgVel, sizeof(int), str);
    saferead(&stepVel, sizeof(int), str);
    saferead(&sizeHistVel, sizeof(int), str);
    saferead(&rangeVel, sizeof(real), str);

    saferead(&stepAvgRdf, sizeof(int), str);
    saferead(&stepRdf, sizeof(int), str);
    saferead(&sizeHistRdf, sizeof(int), str);
    saferead(&rangeRdf, sizeof(real), str);

    saferead(&theta, sizeof(real), str);
    saferead(&usequad, sizeof(bool), str);
    saferead(&nchars, sizeof(int), str);
    options = (string) allocate(nchars * sizeof(char));
    saferead(options, nchars * sizeof(char), str);
    saferead(&tstop, sizeof(real), str);
    saferead(&dtout, sizeof(real), str);
    saferead(&dtoutinfo, sizeof(real), str);
    saferead(&tnow, sizeof(real), str);
    saferead(&tout, sizeof(real), str);
    saferead(&nstep, sizeof(int), str);
    saferead(&rsize, sizeof(real), str);
    saferead(&nbody, sizeof(int), str);

    saferead(&ndim, sizeof(int), str);
	if (ndim != NDIM)
		error("\nrestorestate : ndim = %d; expected %d\n",ndim,NDIM);
    saferead(&seed, sizeof(int), str);
    saferead(&cputotforce, sizeof(real), str);
    saferead(&cputotout, sizeof(real), str);
    saferead(&cputotal, sizeof(real), str);

    saferead(&sTotEnergy, sizeof(real), str);
    saferead(&ssTotEnergy, sizeof(real), str);
    saferead(&sKinEnergy, sizeof(real), str);
    saferead(&ssKinEnergy, sizeof(real), str);
    saferead(&sPressure, sizeof(real), str);
    saferead(&ssPressure, sizeof(real), str);

    saferead(&nchars, sizeof(int), str);
    icfile = (string) allocate(nchars * sizeof(char));
    saferead(icfile, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    icfilefmt = (string) allocate(nchars * sizeof(char));
    saferead(icfilefmt, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    snapoutfile = (string) allocate(nchars * sizeof(char));
    saferead(snapoutfile, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    snapoutfilefmt = (string) allocate(nchars * sizeof(char));
    saferead(snapoutfilefmt, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    statefile = (string) allocate(nchars * sizeof(char));
    saferead(statefile, nchars * sizeof(char), str);
//    saferead(&nchars, sizeof(int), str);
//    restorefile = (string) allocate(nchars * sizeof(char));
//    saferead(restorefile, nchars * sizeof(char), str);

    bodytab = (bodyptr) allocate(nbody * sizeof(body));
    saferead(bodytab, nbody * sizeof(body), str);

    fclose(str);
}

local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"snap-ascii") == 0) *outfmt_int = 0;
    if (strnull(outfmt_str)) *outfmt_int = 1;
    if (strcmp(outfmt_str,"snap-pv") == 0) *outfmt_int = 2;
    if (strcmp(outfmt_str,"snap-bin") == 0) *outfmt_int = 3;
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
		stopflag = 1;
		sprintf(buf,"rm -f %s", stopfilename);
		system(buf);
        printf("\nsaving a stop run state...\n\n");
		cpudt = cputime()-cpuinit;
		cputotal += cpudt;
        savestate(stopsavestate);
		cputotal -= cpudt;
	}
}

#undef stopfilename
#undef stopsavestate

void code_endrun(void)
{
    char   buf[200];
	FILE *fd;
	double cpudt;

    if (! strnull(statefile)) {
        printf("\nsaving final run state...\n\n");
		cpudt = cputime()-cpuinit;
		cputotal += cpudt;
        savestate(savestatetmp);
		cputotal -= cpudt;
        sprintf(buf,"cp %s %s",savestatetmp,statefile);
        system(buf);
    }

        printf("\nremoving temporal files (%s %s)...\n\n",
			saveveltmp, saverdftmp);
		if ((fd=fopen(saveveltmp,"r"))) {
			fclose(fd);
			sprintf(buf,"rm -f %s", saveveltmp);
			system(buf);
		}
		if ((fd=fopen(saverdftmp,"r"))) {
			fclose(fd);
			sprintf(buf,"rm -f %s", saverdftmp);
			system(buf);
		}

	fclose(outlog);
//	printf("\nFinal CPU time : %lf\n",cputime()-cpuinit);
	cputotal += cputime()-cpuinit;
	printf("\nFinal CPU time : %lf\n",cputotal);			// Checar el funcionamiento con restorefile
}

#undef saveveltmp
#undef saverdftmp
#undef savestatetmp
#undef savesnaptmp
