/* =============================================================================
	MODULE: startrun.c					[gbsph]
	Written by: M.A. Rodriguez-Meza
	Fecha de creacion: Enero 2005
	Purpose: Initialize gbsph
	Language: C
	Use: 'startrun();'
	Routines and functions: ReadParameterFile, PrintParameterFile
	Modules, routines and external headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 24, 2007; October 04, 2007:
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "switches.h"

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/vectmath.h"
//#include "../../../General_libs/math/mathfns.h"
//#include "../../../General_libs/general/getparam.h"
//#include "../../../General_libs/math/numrec.h"
//#include "../../../General_libs/general/lic.h"

#include "globaldefs.h"

#include "protodefs.h"

//#include <string.h>

local void force_models_driver(void);

local void ReadParameterFile(char *fname);
local void PrintParameterFile(char *fname);

local void Compute_Parameters(void);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_restorefile(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

local void testdata(void);

void StartRun(string head0, string head1, string head2, string head3)
{
    double cpustart;
	bodyptr p;

    cpustart = cputime();                       

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);
    printf("\t-- Version %s --\n",getversion());
    printf("%s\n", copyright);

//    LicDriver(gd.headline0);

	gd.stopflag = 0;
//	gd.cputotforce = 0;
	gd.cputotout = 0.;
//	gd.cputotal = 0.;


    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();


	StartOutput();

	fprintf(stdout,"\n\nStartRun CPU time: %g\n",cputime()-cpustart);
}


local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
	if (strnull(cmd.restorefile))
		startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-gbsph"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
	cmd.dtimestr = GetParam("dtime");
	cmd.dtoutstr = GetParam("dtout");
	cmd.stepState = GetiParam("stepState");
	cmd.dtoutinfostr = GetParam("dtoutinfo");

    cmd.restorefile = GetParam("restorefile");

	cmd.icfile = GetParam("icfile");    
	cmd.icfilefmt = GetParam("icfilefmt");    
	cmd.snapoutfile = GetParam("snapout");
	cmd.snapoutfilefmt = GetParam("snapoutfmt");
	cmd.statefile = GetParam("statefile");

	cmd.eps = GetdParam("eps");

#if !defined(QUICKSCAN)
	cmd.theta = GetdParam("theta");
#endif
	cmd.usequad = GetbParam("usequad");

	cmd.dm_lambda = GetdParam("dm_lambda");
	cmd.dm_alpha = GetdParam("dm_alpha");
//	cmd.dm_inv_avgphi = GetdParam("dm_inv_avgphi");
	cmd.G = GetdParam("G");
	cmd.dm_a = GetdParam("dm_a");
	cmd.dm_time = GetdParam("dm_time");

	cmd.eps_pot = GetdParam("eps_pot");
	cmd.sigma_pot = GetdParam("sigma_pot");
	cmd.x_pot = GetdParam("x_pot");
	cmd.y_pot = GetdParam("y_pot");
	cmd.z_pot = GetdParam("z_pot");

	cmd.tstop = GetdParam("tstop");

	cmd.forcecalc_method = GetParam("forcecalc_method");
	cmd.force_models = GetParam("force_models");

	cmd.computeTransport = GetbParam("computeTransport"); // Transport ...

	cmd.options = GetParam("options");

	cmd.nbody = GetiParam("nbody");         

	cmd.seed = GetiParam("seed");
}

#undef parameter_null


#define logfile			"gbsph.log"

local void startrun_Common()
{
	real dt1, dt2;
	bodyptr p;

	if (strnull(cmd.restorefile))
		strcpy(gd.mode,"w");
	else
		strcpy(gd.mode,"a");

	if(!(gd.outlog=fopen(logfile,gd.mode)))
		error("\nstart_Common: error opening file '%s' \n",logfile);

    if (strnull(cmd.restorefile)) {
		gd.dtime = (sscanf(cmd.dtimestr, "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtimestr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtime : dt2 must be finite\n");

		gd.pos_pot[0]=cmd.x_pot; gd.pos_pot[1]=cmd.y_pot; gd.pos_pot[2]=cmd.z_pot;

		gd.dtout = (sscanf(cmd.dtoutstr, "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtoutstr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtout : dt2 must be finite\n");

		gd.dtoutinfo=(sscanf(cmd.dtoutinfostr,"%lf/%lf",&dt1,&dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtoutinfostr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtoutinfo :dt2 must be finite\n");

        force_models_driver();
        if (gd.IncludeGrav_Spline) force_setkernel();
        if (! strnull(cmd.icfile)) {
            inputdata();                        
			Compute_Parameters();
        } else {
	    if (cmd.nbody < 1)
                error("startrun: absurd value for nbody\n");

			idum=cmd.seed;
 			Compute_Parameters();
            testdata();
            gd.tnow = 0.0;
        }
        if (cmd.theta < 0.004)
            error("startrun: too small theta value\n");
        gd.rsize = 1.0;				// In some cases with need values below
									// rsize = 1.000001; rsize = 1.001;
        gd.nstep = 0;
        gd.nstep_grav = 0;                          
        gd.tout = gd.tnow;
		gd.toutinfo = gd.tnow;
        for (p = bodytab; p < bodytab+cmd.nbody; p++)
            dm_Mass(p)=Mass(p);
                            
    } else {
        startrun_restorefile();
    }
}

#undef logfile

local void startrun_ParamStat(void)
{
	real dt1, dt2;

	if (GetParamStat("eps") & ARGPARAM)     
		cmd.eps = GetdParam("eps");
#if !defined(QUICKSCAN)
	if (GetParamStat("theta") & ARGPARAM)
		cmd.theta = GetdParam("theta");
#endif
	if (GetParamStat("usequad") & ARGPARAM)
		cmd.usequad = GetbParam("usequad");

	if (GetParamStat("computeTransport") & ARGPARAM)
		cmd.computeTransport = GetbParam("computeTransport");

	if (GetParamStat("options") & ARGPARAM)
		cmd.options = GetParam("options");
	if (GetParamStat("tstop") & ARGPARAM)
		cmd.tstop = GetdParam("tstop");

	if (GetParamStat("dtout") & ARGPARAM) {
		cmd.dtoutstr = GetParam("dtout");
		gd.dtout = (sscanf(cmd.dtoutstr, "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtoutstr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtout : dt2 must be finite\n");
	}

	if (GetParamStat("stepState") & ARGPARAM) 
		cmd.stepState = GetiParam("stepState");

	if (GetParamStat("dtoutinfo") & ARGPARAM) {
		cmd.dtoutinfostr = GetParam("dtoutinfo");
		gd.dtoutinfo = (sscanf(cmd.dtoutinfostr,"%lf/%lf",&dt1,&dt2) == 2 ?
			dt1/dt2 : atof(cmd.dtoutinfostr));
		if ( dt2 == 0. )
			error("\n\nstartrun_ParamStat: dtoutinfo : dt2 must be finite\n");
	}
}


local void CheckParameters(void)
{
	if (!strnull(cmd.restorefile) && !strnull(cmd.icfile))
		fprintf(stdout,"\nCheckParameters: Warning! : %s\n\n",
			"You are using options restorefile and icfile at the same time");

	if (cmd.nbody < 1)
		error("CheckParameters: absurd value for nbody\n");
	if (cmd.stepState < 1)
		error("CheckParameters: absurd value for stepState\n");
	if (gd.dtime == 0.)
		error("\n\nCheckParameters: dtime must be != 0\n");

    if (!(gd.nstep == 0) && !(cmd.tstop - gd.tnow > gd.dtime) ) {
        error("\nNothing to do!...\nCheck nstep, tstop and tnow values...\n\n");
    }
}

#define parameter_restore		"parameters_restore.in-gbsph"

local void startrun_restorefile(void)
{
	real dt1, dt2;

	gd.model_comment = "Restore from state file";
	fprintf(gd.outlog,"\n\nAdded after restart from restart file\n");
	restorestate(cmd.restorefile);
	fprintf(gd.outlog,"\nnbody=%d\n",cmd.nbody);
	fflush(gd.outlog);
	startrun_ParamStat();
	CheckParameters();

	if (scanopt(cmd.options, "new-tout"))
		gd.tout = gd.tnow + gd.dtout;
	if (scanopt(cmd.options, "new-toutinfo"))
		gd.toutinfo = gd.tnow + gd.dtoutinfo;

	Compute_Parameters();

	gd.pos_pot[0]=cmd.x_pot; gd.pos_pot[1]=cmd.y_pot; gd.pos_pot[2]=cmd.z_pot;
	force_models_driver();
    
	if (gd.IncludeGrav_Spline) force_setkernel();
	
	PrintParameterFile(parameter_restore);
}

#undef parameter_restore

local void Compute_Parameters(void)
{
//	stepAvg = (int) (dtout/dtime);
//	gd.stepSnap = (int) (cmd.dtout/cmd.dtime);
//	gd.stepSnap = (int) (gd.dtout/gd.dtime);
//	stepSnap = (int) (freq/freqout);
	gd.stepSnapInit = 0;
}


void ReadParameterFile(char *fname)
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

	SPName(cmd.icfile,"icfile",100);
	SPName(cmd.icfilefmt,"icfilefmt",100);
	SPName(cmd.snapoutfile,"snapout",100);
	SPName(cmd.snapoutfilefmt,"snapoutfmt",100);
	RPName(cmd.eps,"eps");
	SPName(cmd.options,"options",100);
	SPName(cmd.dtimestr,"dtime",100);
	RPName(cmd.tstop,"tstop");
	SPName(cmd.dtoutstr,"dtout",100);
	SPName(cmd.dtoutinfostr,"dtoutinfo",100);
	IPName(cmd.nbody,"nbody");
	IPName(cmd.seed,"seed");
	SPName(cmd.statefile,"statefile",100);
	IPName(cmd.stepState,"stepState");
	SPName(cmd.restorefile,"restorefile",100);

	BPName(cmd.computeTransport,"computeTransport");	// Transport ...

	RPName(cmd.theta,"theta");
	BPName(cmd.usequad,"usequad");

	SPName(cmd.forcecalc_method,"forcecalc_method",100);
	SPName(cmd.force_models,"force_models",100);
	RPName(cmd.dm_lambda,"dm_lambda");
	RPName(cmd.dm_alpha,"dm_alpha");
//	RPName(cmd.dm_inv_avgphi,"dm_inv_avgphi");
	RPName(cmd.G,"G");
	RPName(cmd.dm_a,"dm_a");
	RPName(cmd.dm_time,"dm_time");

	RPName(cmd.eps_pot,"eps_pot");
	RPName(cmd.sigma_pot,"sigma_pot");

	RPName(cmd.x_pot,"x_pot");
	RPName(cmd.y_pot,"y_pot");
	RPName(cmd.z_pot,"z_pot");

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;

            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                        *buf2=(char)NULL;

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
                            } else
                                error("getbparam: %s=%s not bool\n", buf1, buf2);
                        break;
                }
            } else 
                error("Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
                        fname,buf1);
        }
        fclose(fd);
    } else
      error("Parameter file %s not found.\n\n", fname);
  
    for(i=0;i<nt;i++)
        if(*tag[i])
            error("Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                    tag[i],fname);

#undef DOUBLE 
#undef STRING 
#undef INT 
#undef MAXTAGS
}

#define FMTT	"%-35s%s\n"
#define FMTI	"%-35s%d\n"
#define FMTR	"%-35s%g\n"
#define FMTLI	"%-35s%ld\n"

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
        fprintf(fdout,"%s %s\n","% Parameter input file for:",gd.headline0);
        fprintf(fdout,"%s %s\n","%",gd.headline);
        fprintf(fdout,"%s %s\n","%",gd.model_comment);
        fprintf(fdout,"%s %s\n","%",copyright);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s\n%s\n",
            "%-------------------------------------------------------------------",
            "%");
        fprintf(fdout,FMTT,"forcecalc_method",cmd.forcecalc_method);
        fprintf(fdout,FMTT,"force_models",cmd.force_models);
        fprintf(fdout,FMTT,"icfile",cmd.icfile);
        fprintf(fdout,FMTT,"icfilefmt",cmd.icfilefmt);
        fprintf(fdout,FMTT,"snapout",cmd.snapoutfile);
        fprintf(fdout,FMTT,"snapoutfmt",cmd.snapoutfilefmt);
        fprintf(fdout,FMTT,"statefile",cmd.statefile);
        fprintf(fdout,FMTI,"stepState",cmd.stepState);
        fprintf(fdout,FMTT,"restorefile",cmd.restorefile);

        fprintf(fdout,FMTT,"computeTransport",
				cmd.computeTransport ? "true" : "false");

        fprintf(fdout,FMTR,"theta",cmd.theta);
        fprintf(fdout,FMTT,"usequad",cmd.usequad ? "true" : "false");
        fprintf(fdout,FMTR,"dm_lambda",cmd.dm_lambda);
        fprintf(fdout,FMTR,"dm_alpha",cmd.dm_alpha);
//        fprintf(fdout,FMTR,"dm_inv_avgphi",cmd.dm_inv_avgphi);
        fprintf(fdout,FMTR,"G",cmd.G);
        fprintf(fdout,FMTR,"dm_a",cmd.dm_a);
        fprintf(fdout,FMTR,"dm_time",cmd.dm_time);

        fprintf(fdout,FMTR,"eps_pot",cmd.eps_pot);
        fprintf(fdout,FMTR,"sigma_pot",cmd.sigma_pot);
        fprintf(fdout,FMTR,"x_pot",cmd.x_pot);
        fprintf(fdout,FMTR,"y_pot",cmd.y_pot);
        fprintf(fdout,FMTR,"z_pot",cmd.z_pot);
        fprintf(fdout,FMTT,"dtime",cmd.dtimestr);
        fprintf(fdout,FMTR,"eps",cmd.eps);
        fprintf(fdout,FMTR,"tstop",cmd.tstop);
        fprintf(fdout,FMTT,"dtout",cmd.dtoutstr);
        fprintf(fdout,FMTT,"dtoutinfo",cmd.dtoutinfostr);
        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,FMTI,"nbody",cmd.nbody);
        fprintf(fdout,FMTLI,"seed",cmd.seed);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

local void force_models_driver(void)
{
    string force_models_input;

    gd.IncludeGrav=0;
    gd.IncludeGrav_Plummer=0;
    gd.IncludeGrav_Spline=0;
    gd.IncludeSF_eps=0;
    force_models_input=GetParam("force_models");
    if(scanopt(cmd.force_models,"gravity-plummer")) gd.IncludeGrav_Plummer=1;
    if(scanopt(cmd.force_models,"gravity-spline")) gd.IncludeGrav_Spline=1;
    if (gd.IncludeGrav_Plummer && gd.IncludeGrav_Spline)
        error("Gravity Plummer and Spline not at the same time!...");
    if(scanopt(cmd.force_models,"scalar-field-potential-eps")) gd.IncludeSF_eps=1;
    if(strnull(force_models_input)) gd.IncludeGrav_Plummer=1;

    if (gd.IncludeGrav_Plummer || gd.IncludeGrav_Spline) gd.IncludeGrav=1;
}

#define MFRAC  0.999  

local void testdata(void)
{
    real rsc, vsc, r, v, x, y;
    vector rcm, vcm;
    bodyptr p;

    gd.model_comment = "3D Plummer Model";
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
                                                
    rsc = (3 * PI) / 16;                        
    vsc = rsqrt(1.0 / rsc);                     
    CLRV(rcm);                                  
    CLRV(vcm);                                  
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        Type(p) = BODY;                         
        Mass(p) = 1.0 / cmd.nbody;                  
        x = xrandom(0.0, MFRAC);                
        r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1); 
        pickshell(Pos(p), NDIM, rsc * r);       
        do {                                    
            x = xrandom(0.0, 1.0);              
            y = xrandom(0.0, 0.1);              
        } while (y > x*x * rpow(1 - x*x, 3.5)); 
        v = x * rsqrt(2.0 / rsqrt(1 + r*r));    
        pickshell(Vel(p), NDIM, vsc * v);       
        ADDMULVS(rcm, Pos(p), 1.0 / cmd.nbody);     
        ADDMULVS(vcm, Vel(p), 1.0 / cmd.nbody);     
    }
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        SUBV(Pos(p), Pos(p), rcm);              
        SUBV(Vel(p), Vel(p), vcm);              
    }

}
#undef MFRAC

