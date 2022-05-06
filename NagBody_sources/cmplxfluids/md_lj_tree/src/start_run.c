/* ==============================================================================
!	MODULE: start_run.c															!
!	Written by: M.A. Rodr'iguez-Meza.											!
!	Fecha de creaci'on: Enero 2005.												!
!	Purpose: Initialize md_lj_tree												!
!	Language: C																	!
!	Use: 'startrun();'															!
!	Routines and functions: testdata, ReadParameterFile, PrintParameterFile		!
!	Modules, routines and external headers:										!
!	Coments and notes:															!
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
#include "machines.h"
#include "lic.h"
#include "global_defs.h"
#include "proto_defs.h"
#include <string.h>		// Incluido para quitar el warning: 
				// "Implicit declaration of built-in function 'strcpy' y 'strchr'"

local void Compute_nbody(void);
local void Compute_Box_Units(void);
local void Compute_Parameters(void);
local void Print_Parameters(void);
local void testdata_sc_random(void);
local void testdata_sc(void);
local void testdata(void);
local void ReadParameterFile(char *);
local void PrintParameterFile(char *);
local void forcecalc_method_string_to_int(string,int *);
local void icfilefmt_string_to_int(string,int *);

//local int seed;

local vector N;							// Needed in testdata_sc

#define parameter_null	"parameters_null"
#define logfile			"md_lj_tree.log"

void startrun(void)
{
	real dt1, dt2;

//    lic();

	cputotforce = 0;					// Initializing cpu time counters
	cputotout = 0.;						// If !strnull(restorefile)
	cputotal = 0.;						// these values are modified below

    paramfile = GetParam("paramfile");
    if (!strnull(paramfile)) {			// Parameter file is read
        ReadParameterFile(paramfile);
		if (strnull(restorefile))
			strcpy(mode,"w");
		else
			strcpy(mode,"a");
		if(!(outlog=fopen(logfile,mode))) {
			fprintf(stdout,"error opening file '%s' \n",logfile);
			exit(0);
		}

		if (!strnull(restorefile) && !strnull(icfile))
			error("\nstartrun : Use only one of the options restorefile or icfile\n\n");

        if (strnull(restorefile)) {
            if (! strnull(icfile)) {
				icfilefmt_string_to_int(icfilefmt, &icfilefmt_int);
                inputdata();
				Compute_Box_Units();
				Compute_Parameters();
            } else {                                 
                if (nbody < 1)                      
                    error("startrun: absurd value for nbody\n");
#ifdef VISUALC
                srand(seed);			// on digital C
#else
                srandom(seed);			// on unix
#endif
				Compute_nbody();
				Compute_Box_Units();
                testdata();                       
				Compute_Parameters();
                tnow = 0.0;                         
            }
			rsize = 1.0;                            
            nstep = 0;
            tout = tnow;
            toutinfo = tnow;
        } else {										// Reading restorefile
			fprintf(outlog,"\n\nAdded after restart from restart file\n");
            restorestate(restorefile);
			fprintf(outlog,"\nnbody=%d\n",nbody);
            if (GetParamStat("temperature") & ARGPARAM)
                temperature = GetdParam("temperature");
            if (GetParamStat("density") & ARGPARAM) 
                density = GetdParam("density");
			if (GetParamStat("stepEquil") & ARGPARAM)
				stepEquil = GetiParam("stepEquil");         
			if (GetParamStat("stepAvgVel") & ARGPARAM)
				stepAvgVel = GetiParam("stepAvgVel");         
			if (GetParamStat("stepVel") & ARGPARAM)
				stepVel = GetiParam("stepVel");         
			if (GetParamStat("sizeHistVel") & ARGPARAM)
				sizeHistVel = GetiParam("sizeHistVel");
			if (GetParamStat("rangeVel") & ARGPARAM)
				rangeVel = GetdParam("rangeVel");
			if (GetParamStat("stepAvgRdf") & ARGPARAM)
				stepAvgRdf = GetiParam("stepAvgRdf");         
			if (GetParamStat("stepRdf") & ARGPARAM)
				stepRdf = GetiParam("stepRdf");         
			if (GetParamStat("sizeHistRdf") & ARGPARAM)
				sizeHistRdf = GetiParam("sizeHistRdf");
			if (GetParamStat("rangeRdf") & ARGPARAM)
				rangeRdf = GetdParam("rangeRdf");
            if (GetParamStat("options") & ARGPARAM)
                options = GetParam("options");
			if (GetParamStat("theta") & ARGPARAM) 
				theta = GetdParam("theta");
			if (GetParamStat("usequad") & ARGPARAM) 
				usequad = GetbParam("usequad");
            if (GetParamStat("tstop") & ARGPARAM)
                tstop = GetdParam("tstop");
            if (GetParamStat("dtout") & ARGPARAM)
				dtout = (sscanf(GetParam("dtout"), "%lf/%lf", &dt1, &dt2) == 2 ?
					dt1/dt2 : GetdParam("dtout"));
            if (scanopt(options, "new-tout"))
				tout = tnow + dtout;
            if (GetParamStat("dtoutinfo") & ARGPARAM)
				dtoutinfo = (sscanf(GetParam("dtoutinfo"), "%lf/%lf", &dt1, &dt2) == 2 ?
					dt1/dt2 : GetdParam("dtoutinfo"));
            if (scanopt(options, "new-toutinfo"))
				toutinfo = tnow + dtoutinfo;
			Compute_Box_Units();
			Compute_Parameters();
        }
		forcecalc_method_string_to_int(forcecalc_method, &forcemethod_int);
		Print_Parameters();
        PrintParameterFile(paramfile);
        return;
    }								// Finish reading parameterfile

    icfile = GetParam("icfile");
    icfilefmt = GetParam("icfilefmt");
    snapoutfile = GetParam("snapout");
    snapoutfilefmt = GetParam("snapoutfmt");
    statefile = GetParam("statefile");
    restorefile = GetParam("restorefile");

	if (!strnull(restorefile) && !strnull(icfile))
		error("\nstartrun : Use only one of the options restorefile or icfile\n\n");

	if (strnull(restorefile))
		strcpy(mode,"w");
	else
		strcpy(mode,"a");
	if(!(outlog=fopen(logfile,mode))) {
		fprintf(stdout,"error opening file '%s' \n",logfile);
		exit(0);
	}
    if (strnull(restorefile)) {
        temperature = GetdParam("temperature");
        density = GetdParam("density");
		stepEquil = GetiParam("stepEquil");         
		stepAvgVel = GetiParam("stepAvgVel");
		stepVel = GetiParam("stepVel");
		sizeHistVel = GetiParam("sizeHistVel");
		rangeVel = GetdParam("rangeVel");
		stepAvgRdf = GetiParam("stepAvgRdf");
		stepRdf = GetiParam("stepRdf");
		sizeHistRdf = GetiParam("sizeHistRdf");
		rangeRdf = GetdParam("rangeRdf");
        theta = GetdParam("theta");
        usequad = GetbParam("usequad");
		dtime = (sscanf(GetParam("dtime"), "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : GetdParam("dtime"));
        tstop = GetdParam("tstop");
		dtout = (sscanf(GetParam("dtout"), "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : GetdParam("dtout"));
		dtoutinfo = (sscanf(GetParam("dtoutinfo"), "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : GetdParam("dtoutinfo"));
        forcecalc_method = GetParam("forcecalc_method");
        options = GetParam("options");
        if (! strnull(icfile)) {
			icfilefmt_string_to_int(icfilefmt, &icfilefmt_int);
            inputdata();
			Compute_Box_Units();
			Compute_Parameters();
            seed=GetiParam("seed");	// to always have defaults
        } else {                                 
            nbody = GetiParam("nbody");         
            if (nbody < 1)                      
                error("startrun: absurd value for nbody\n");
            seed = GetiParam("seed");         
#ifdef VISUALC
            srand(seed);			// on digital C
#else
            srandom(seed);			// on unix
#endif
			Compute_nbody();
			Compute_Box_Units();
            testdata();                       
			Compute_Parameters();
            tnow = 0.0;
        }
        rsize = 1.0;                            
        nstep = 0;                              
        tout = tnow;
        toutinfo = tnow;
    } else {										// Reading restorefile
		fprintf(outlog,"\n\nAdded after restart from restart file\n");
        restorestate(restorefile);
		fprintf(outlog,"\nnbody=%d\n",nbody);
        if (GetParamStat("temperature") & ARGPARAM) 
            temperature = GetdParam("temperature");
        if (GetParamStat("density") & ARGPARAM) 
            density = GetdParam("density");
		if (GetParamStat("stepEquil") & ARGPARAM)
			stepEquil = GetiParam("stepEquil");         
		if (GetParamStat("stepAvgVel") & ARGPARAM)
			stepAvgVel = GetiParam("stepAvgVel");         
		if (GetParamStat("stepVel") & ARGPARAM)
			stepVel = GetiParam("stepVel");         
		if (GetParamStat("sizeHistVel") & ARGPARAM)
			sizeHistVel = GetiParam("sizeHistVel");
		if (GetParamStat("rangeVel") & ARGPARAM)
			rangeVel = GetdParam("rangeVel");
		if (GetParamStat("stepAvgRdf") & ARGPARAM)
			stepAvgRdf = GetiParam("stepAvgRdf");         
		if (GetParamStat("stepRdf") & ARGPARAM)
			stepRdf = GetiParam("stepRdf");         
		if (GetParamStat("sizeHistRdf") & ARGPARAM)
			sizeHistRdf = GetiParam("sizeHistRdf");
		if (GetParamStat("rangeRdf") & ARGPARAM)
			rangeRdf = GetdParam("rangeRdf");
        if (GetParamStat("options") & ARGPARAM)
            options = GetParam("options");
		if (GetParamStat("theta") & ARGPARAM) 
            theta = GetdParam("theta");
		if (GetParamStat("usequad") & ARGPARAM) 
            usequad = GetbParam("usequad");
        if (GetParamStat("tstop") & ARGPARAM)
            tstop = GetdParam("tstop");
        if (GetParamStat("dtout") & ARGPARAM)
			dtout = (sscanf(GetParam("dtout"), "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : GetdParam("dtout"));
        if (scanopt(options, "new-tout"))
            tout = tnow + dtout;
        if (GetParamStat("dtoutinfo") & ARGPARAM)
			dtoutinfo = (sscanf(GetParam("dtoutinfo"), "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : GetdParam("dtoutinfo"));
        if (scanopt(options, "new-toutinfo"))
            toutinfo = tnow + dtoutinfo;
		Compute_Box_Units();
		Compute_Parameters();
    }
    forcecalc_method_string_to_int(forcecalc_method, &forcemethod_int);
	Print_Parameters();
    if (strnull(paramfile))
        PrintParameterFile(parameter_null);
}

#undef logfile

local void Compute_Parameters(void)
{
	real rri, rri3;
	int k;

	Rcut = rpow(2.0, 1.0/6.0)*sigma;	// Cut radius
	RcutSq = Rcut*Rcut;

	stepAvg = (int) (dtout/dtime);
	stepSnap = stepAvg;
	stepSnapInit = 0;
	vMag = rsqrt( NDIM*(1.0 -			// Thermal velocity
			1.0/nbody)*kB*temperature/mass );
	fphi = 4.0*eps;						// Needed in ljforcecalc
	ssq = sigma*sigma;					// Needed in ljforcecalc
	fa = 48.0*eps/ssq;					// Needed in ljforcecalc

	MULVS(cells, Box, 1.0/Rcut);		// Only needed for cellmethod
	AllocMem(cellList, VProd (cells)	// Only needed for cellmethod
				+ nbody, int);
	
	lAvg = 0.; DO_COORD(k) lAvg += Box[k];
	lAvg = (lAvg/NDIM) / rpow(nbody, 1./3.);
	dtcrit = lAvg/vMag;
	rri=ssq/rsqr(lAvg); rri3=rri*rri*rri;
	fAvg = fa*rri3*(rri3-0.5)*rri*lAvg;
}

local void Compute_nbody(void)
{
	fprintf(outlog,"\n\nInitial nbody=%d",nbody);
	if (NDIM == 3) {
		N[0] = (int) rpow(nbody,1./3.);		// nbody = N^3
		N[1] = N[0]; N[2] = N[0];
		nbody = N[0]*N[1]*N[2];
	} else if (NDIM == 2) {
		N[0] = (int) rsqrt(nbody);			// nbody = N^2
		N[1] = N[0];
		nbody = N[0]*N[1];
	} else error("\n\nWrong NDIM!\n\n");
}

local void Compute_Box_Units(void)
{
	eps = 1.0; sigma = 1.0;	mass = 1.0;		// Units system
	kB = 1.0;

	if (NDIM == 3) {						// Size of box side x
		Box[0] = rpow(((real) nbody * mass)/density,1.0/3.0); 
		Box[1] = Box[0]; Box[2] = Box[0];
	} else if (NDIM == 2) {					// Size of box side x
		Box[0] = rsqrt(((real) nbody * mass)/density); 
		Box[1] = Box[0];
	} else error("\n\nWrong NDIM!\n\n");
}

local void Print_Parameters(void)
{
	fprintf(outlog,"\nParameters:\n");
#if (NDIM==3)
	fprintf(outlog,"\nBox: %g %g %g",Box[0],Box[1],Box[2]);
	fprintf(outlog,"\ncells: %d %d %d",cells[0],cells[1],cells[2]);
#else
#if (NDIM==2)
	fprintf(outlog,"\nBox: %g %g",Box[0],Box[1]);
	fprintf(outlog,"\ncells: %d %d",cells[0],cells[1]);
#endif
#endif
	fprintf(outlog,"\nsize of cellList : %d\n",VProd (cells)+nbody);
	fprintf(outlog,"\neps= %g",eps);
	fprintf(outlog,"\nsigma= %g",sigma);
	fprintf(outlog,"\ntemperature= %g",temperature);
	fprintf(outlog,"\nmass= %g",mass);
	fprintf(outlog,"\ndensity= %g",density);
	fprintf(outlog,"\nRcut= %g",Rcut);
	fprintf(outlog,"\nRcutSq= %g",RcutSq);
	fprintf(outlog,"\nstepAvg= %d\n\n",stepAvg);
	fprintf(outlog,"\nlAvg= %g",lAvg);
	fprintf(outlog,"\ndtcrit= %g",dtcrit);
	fprintf(outlog,"\nfAvg= %g",fAvg);
	fflush(outlog);
}

#undef parameter_null

local void testdata(void)
{
	testdata_sc();
//	testdata_sc_random();
}


// Checar que la siguiente rutina cumple con Compute_nbody, Compute_Box_Units, y
// que la caja esta entre -Box/2 y Box/2 ...
local void testdata_sc_random(void)
{
    vector vcm;
    bodyptr p;
	real Ekin, ScKin, velsq, tmass;
	int k;

    model_comment = "Random Cubic Box Model";

    bodytab = (bodyptr) allocate(nbody * sizeof(body));

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	DO_BODY(p, bodytab, bodytab+nbody) {
		Id(p) = p-bodytab+1;
        Type(p) = BODY;
        Mass(p) = mass;
		DO_COORD(k)
			Pos(p)[k]	= xrandom(-0.5*Box[k], 0.5*Box[k]);	// Box center is (0,0,0)
		DO_COORD(k)
			Vel(p)[k] = grandom(0.0,1.0);
		tmass += Mass(p);
        ADDMULVS(vcm, Vel(p), Mass(p));
    }
	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }

	Ekin *= 0.5;
	fprintf(outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	if (NDIM == 3)
		ScKin = rsqrt(temperature* (real)(3*nbody-3) /Ekin );
	else if (NDIM == 2)
		ScKin = rsqrt(temperature* (real)(2*nbody-2) /Ekin );
	else error("\n\nWrong NDIM!\n\n");

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+nbody) {
        MULVS(Vel(p), Vel(p), ScKin);
        DOTVP(velsq, Vel(p), Vel(p));           
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(outlog,"\nEkin after=%g\n",Ekin);
}

local void testdata_sc(void)
{
    vector vcm;
    bodyptr p;
	real Ekin, ScKin, velsq, tmass;
	int k,n, i,j,l;
	vector delta, c;
	real f, os;

    model_comment = "Simple Cubic Box Model";

    bodytab = (bodyptr) allocate(nbody * sizeof(body));

	os =0.05*Box[0];						// offset
	f = 1.0 - 2.0*os/Box[0];				// scaling delta to 1 - 2 os/Lk

	DO_COORD(k)
		delta[k]=f*Box[k]/((real) (N[k]-1));

	if (NDIM == 3)
		fprintf(outlog,"\ndelta: %g %g %g\n",
				delta[0],delta[1],delta[2]);
	else if (NDIM == 2)
		fprintf(outlog,"\ndelta: %g %g\n",
				delta[0],delta[1]);
	else error("\n\nWrong NDIM!\n\n");

	p = bodytab;
	if (NDIM == 3)
		for (i=1; i<=N[0]; i++) {
			c[0] = ((real) (i-1))*delta[0] + os;
			for (j=1; j<=N[1]; j++) {
				c[1] = ((real) (j-1))*delta[1] + os;
				for (l=1; l<=N[2]; l++) {
					c[2] = ((real) (l-1))*delta[2] + os;
					SETV(Pos(p),c);
					++p;
				}
			}
		}
	else if (NDIM == 2)
		for (i=1; i<=N[0]; i++) {
			c[0] = ((real) (i-1))*delta[0] + os;
			for (j=1; j<=N[1]; j++) {
				c[1] = ((real) (j-1))*delta[1] + os;
				SETV(Pos(p),c);
				++p;
			}
		}
	else error("\n\nWrong NDIM!\n\n");
	fprintf(outlog,"\nCreated bodies = %d",nbody);

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	DO_BODY(p, bodytab, bodytab+nbody) {
		Id(p) = p-bodytab+1;
        Type(p) = BODY;
        Mass(p) = mass;
		DO_COORD(k)
			Vel(p)[k] = grandom(0.0,1.0);
		tmass += Mass(p);
        ADDMULVS(vcm, Vel(p), Mass(p));
    }
	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	if (NDIM == 3)
		ScKin = rsqrt(temperature* (real)(3*nbody-3) /Ekin );
	else if (NDIM == 2)
		ScKin = rsqrt(temperature* (real)(2*nbody-2) /Ekin );
	else error("\n\nWrong NDIM!\n\n");

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+nbody) {
        MULVS(Vel(p), Vel(p), ScKin);
        DOTVP(velsq, Vel(p), Vel(p));           
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(outlog,"\nEkin after=%g\n",Ekin);

	DO_BODY(p, bodytab, bodytab+nbody)
		VVSAdd(Pos(p), -0.5, Box);		// Box center is (0,0,0)

}

local void ReadParameterFile(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define BOOLEAN 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200];
  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

  nt=0;

//	SPName(icfile,1);
  strcpy(tag[nt],"icfile"); 
  icfile = (string) malloc(1);
  addr[nt]=icfile;
  id[nt++]=STRING;

//	SPName(icfilefmt,1);
  strcpy(tag[nt],"icfilefmt");
  icfilefmt = (string) malloc(1);
  addr[nt]=icfilefmt;
  id[nt++]=STRING;

  strcpy(tag[nt],"snapout"); 
  snapoutfile=(string) malloc(100);
  addr[nt]=snapoutfile;
  id[nt++]=STRING;

  strcpy(tag[nt],"snapoutfmt"); 
  snapoutfilefmt=(string) malloc(100);
  addr[nt]=snapoutfilefmt;
  id[nt++]=STRING;

//	RPName(dtime);
  strcpy(tag[nt],"dtime"); 
  addr[nt]=&dtime;
  id[nt++]=DOUBLE;
  
//  RPName(theta);
  strcpy(tag[nt],"theta"); 
  addr[nt]=&theta;
  id[nt++]=DOUBLE;

//	BPName(usequad);
  strcpy(tag[nt],"usequad"); 
  addr[nt]=&usequad;
  id[nt++]=BOOLEAN;


  strcpy(tag[nt],"temperature");
  addr[nt]=&temperature;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"density");
  addr[nt]=&density;
  id[nt++]=DOUBLE;

//	IPName(stepEquil);
  strcpy(tag[nt],"stepEquil"); 
  addr[nt]=&stepEquil;
  id[nt++]=INT;

  strcpy(tag[nt],"stepAvgVel"); 
  addr[nt]=&stepAvgVel;
  id[nt++]=INT;

  strcpy(tag[nt],"stepVel"); 
  addr[nt]=&stepVel;
  id[nt++]=INT;

  strcpy(tag[nt],"sizeHistVel"); 
  addr[nt]=&sizeHistVel;
  id[nt++]=INT;

  strcpy(tag[nt],"rangeVel"); 
  addr[nt]=&rangeVel;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"stepAvgRdf"); 
  addr[nt]=&stepAvgRdf;
  id[nt++]=INT;

  strcpy(tag[nt],"stepRdf"); 
  addr[nt]=&stepRdf;
  id[nt++]=INT;

  strcpy(tag[nt],"sizeHistRdf"); 
  addr[nt]=&sizeHistRdf;
  id[nt++]=INT;

  strcpy(tag[nt],"rangeRdf"); 
  addr[nt]=&rangeRdf;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"options"); 
  options=(string)malloc(100);
  addr[nt]=options;
  id[nt++]=STRING;

  strcpy(tag[nt],"forcecalc_method"); 
  forcecalc_method=(string)malloc(100);
  addr[nt]=forcecalc_method;
  id[nt++]=STRING;

  strcpy(tag[nt],"tstop"); 
  addr[nt]=&tstop;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"dtout"); 
  addr[nt]=&dtout;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"dtoutinfo"); 
  addr[nt]=&dtoutinfo;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"nbody"); 
  addr[nt]=&nbody;
  id[nt++]=INT;

  strcpy(tag[nt],"seed"); 
  addr[nt]=&seed;
  id[nt++]=INT;

  strcpy(tag[nt],"statefile"); 
  statefile=(string) malloc(100);
  addr[nt]=statefile;
  id[nt++]=STRING;

  strcpy(tag[nt],"restorefile"); 
  restorefile=(string) malloc(100);
  addr[nt]=restorefile;
  id[nt++]=STRING;      

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
//                *buf2=NULL;							// Original
                *buf2=(int) NULL;						// With cast
            if(buf1[0]=='%')
		continue;
            for(i=0,j=-1;i<nt;i++)
                if(strcmp(buf1,tag[i])==0) {
                    j=i;
		    tag[i][0]=0;
		    break;
            }
            if(j>=0) {
                switch(id[j]) {
		    case DOUBLE:
		      *((double*)addr[j])=atof(buf2); 
		      break;
		    case STRING:
		      strcpy(addr[j],buf2);
		      break;
		    case INT:
		      *((int*)addr[j])=atoi(buf2);
		      break;
		    case BOOLEAN:
                        if (strchr("tTyY1", *buf2) != NULL) {          
                            *((bool*)addr[j])=TRUE;
                        } else 
                            if (strchr("fFnN0", *buf2) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("GetbParam: %s=%s not bool\n", buf1, buf2);
                            }
		      break;
                }
            } else {
                fprintf(stdout,
                    "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
                    fname,buf1);
                errorFlag=1;
            }
        }					// end of while
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1); 
    }						// end of the if((fd=fopen...))
  
    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            fprintf(stdout,
                "Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                tag[i],fname);
            exit(0);
        }
    }
#undef DOUBLE 
#undef STRING 
#undef INT 
#undef BOOLEAN
#undef MAXTAGS
}

local void PrintParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];
    
    sprintf(buf,"%s%s",fname,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
            "%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Parameter input file for:",headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",headline1,headline2,"%",
            headline3);
        fprintf(fdout,"%s\n%s\n",
            "%-------------------------------------------------------------------",
            "%");
        fprintf(fdout,"%-35s%s\n","forcecalc_method",forcecalc_method);
        fprintf(fdout,"%-35s%s\n","icfile",icfile);
        fprintf(fdout,"%-35s%s\n","icfilefmt",icfilefmt);
        fprintf(fdout,"%-35s%s\n","snapout",snapoutfile);
        fprintf(fdout,"%-35s%s\n","snapoutfmt",snapoutfilefmt);
        fprintf(fdout,"%-35s%s\n","statefile",statefile);
        fprintf(fdout,"%-35s%s\n","restorefile",restorefile);
        fprintf(fdout,"%-35s%g\n","theta",theta);
        fprintf(fdout,"%-35s%s\n","usequad",usequad ? "true" : "false");
        fprintf(fdout,"%-35s%g\n","temperature",temperature);
        fprintf(fdout,"%-35s%g\n","density",density);
        fprintf(fdout,"%-35s%d\n","stepEquil",stepEquil);
        fprintf(fdout,"%-35s%d\n","stepAvgVel",stepAvgVel);
        fprintf(fdout,"%-35s%d\n","stepVel",stepVel);
        fprintf(fdout,"%-35s%d\n","sizeHistVel",sizeHistVel);
        fprintf(fdout,"%-35s%g\n","rangeVel",rangeVel);
        fprintf(fdout,"%-35s%d\n","stepAvgRdf",stepAvgRdf);
        fprintf(fdout,"%-35s%d\n","stepRdf",stepRdf);
        fprintf(fdout,"%-35s%d\n","sizeHistRdf",sizeHistRdf);
        fprintf(fdout,"%-35s%g\n","rangeRdf",rangeRdf);
        fprintf(fdout,"%-35s%g\n","dtime",dtime);
        fprintf(fdout,"%-35s%g\n","tstop",tstop);
        fprintf(fdout,"%-35s%g\n","dtout",dtout);
        fprintf(fdout,"%-35s%g\n","dtoutinfo",dtoutinfo);
        fprintf(fdout,"%-35s%s\n","options",options);
        fprintf(fdout,"%-35s%d\n","nbody",nbody);
        fprintf(fdout,"%-35s%d\n\n","seed",seed);	
							// is important "\n\n" so can be read as
    }						// input parameter file
    fclose(fdout);
}

local void forcecalc_method_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"barnes") == 0) *method_int = 0;
    if (strnull(method_str)) *method_int = 1;
    if (strcmp(method_str,"normal") == 0) *method_int = 2;
    if (strcmp(method_str,"nblist") == 0) *method_int = 3;
    if (strcmp(method_str,"direct") == 0) *method_int = 4;
    if (strcmp(method_str,"cells") == 0) *method_int = 5;
    if (strcmp(method_str,"direct2") == 0) *method_int = 6;
    if (strcmp(method_str,"normal2") == 0) *method_int = 7;
    if (strcmp(method_str,"barnes2") == 0) *method_int = 8;
}

local void icfilefmt_string_to_int(string icfmt_str,int *infmt_int)
{
    *infmt_int=-1;
    if (strcmp(icfmt_str,"snap-ascii") == 0) *infmt_int = 0;
    if (strnull(icfmt_str)) *infmt_int = 1;
    if (strcmp(icfmt_str,"snap-pv") == 0) *infmt_int = 2;
    if (strcmp(icfmt_str,"snap-bin") == 0) *infmt_int = 3;
}

