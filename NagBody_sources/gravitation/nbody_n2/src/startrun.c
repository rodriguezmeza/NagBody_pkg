/*==============================================================================
	MODULE: startrun.c				[nbody_n2]
	Written by: Mario A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Initialize nbody_n2
	Language: C
	Use: 'startrun();'
	Routines and functions: testdata, ReadParameterFile, PrintParameterFile
	Modules, routines and external headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 2007; November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

local void testdata(void);
local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_restorefile(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);


void StartRun(string head0, string head1, string head2, string head3)
{
    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

	gd.stopflag = 0;

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

	StartOutput();
}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
	if (strnull(cmd.restorefile))
		startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-nbody_n2"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
    cmd.icfile = GetParam("icfile");
    cmd.icfilefmt = GetParam("icfilefmt");
    cmd.snapoutfile = GetParam("snapout");
    cmd.snapoutfilefmt = GetParam("snapoutfmt");
    cmd.statefile = GetParam("statefile");
	cmd.stepState = GetiParam("stepState");
    cmd.restorefile = GetParam("restorefile");
	cmd.eps = GetdParam("eps");
	cmd.dtimestr = GetParam("dtime");
	cmd.tstop = GetdParam("tstop");
	cmd.dtoutstr = GetParam("dtout");
	cmd.dtoutinfostr = GetParam("dtoutinfo");
	cmd.options = GetParam("options");
	cmd.seed = GetiParam("seed");
	cmd.nbody = GetiParam("nbody");
	PrintParameterFile(parameter_null); // Esta arriba ... quitar ...
}

#undef parameter_null

#define logfile			"nbody_n2.log"

local void startrun_Common(void)
{
	real dt1, dt2;

	if (strnull(cmd.restorefile))
		strcpy(gd.mode,"w");
	else
		strcpy(gd.mode,"a");

	if(!(gd.outlog=fopen(logfile, gd.mode)))
		error("\nstart_Common: error opening file '%s' \n",logfile);

	if (strnull(cmd.restorefile)) {
		gd.dtime = (sscanf(cmd.dtimestr, "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtimestr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtime : dt2 must be finite\n");

		gd.dtout = (sscanf(cmd.dtoutstr, "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtoutstr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtout : dt2 must be finite\n");

		gd.dtoutinfo=(sscanf(cmd.dtoutinfostr,"%lf/%lf",&dt1,&dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtoutinfostr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtoutinfo :dt2 must be finite\n");

		if (! strnull(cmd.icfile)) {
			inputdata();
			CheckParameters();
		} else {
			idum=cmd.seed;
			xsrandom(idum);
			CheckParameters();
			testdata();
			gd.tnow = 0.0;
		}
		gd.nstep = 0;
		gd.tout = gd.tnow;
		gd.toutinfo = gd.tnow;
	} else
		startrun_restorefile();
}

#undef logfile

local void startrun_restorefile(void)
{
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
}

local void startrun_ParamStat(void)
{
	real dt1, dt2;

	if (GetParamStat("eps") & ARGPARAM) 
		cmd.eps = GetdParam("eps");
	if (GetParamStat("stepState") & ARGPARAM) 
		cmd.stepState = GetiParam("stepState");
	if (GetParamStat("options") & ARGPARAM)
		cmd.options = GetParam("options");
	if (GetParamStat("tstop") & ARGPARAM)
		cmd.tstop = GetdParam("tstop");
	if (GetParamStat("dtime") & ARGPARAM) {
		cmd.dtimestr = GetParam("dtime");
		gd.dtime = (sscanf(cmd.dtimestr, "%lf/%lf", &dt1, &dt2) == 2 ?
			dt1/dt2 : atof(cmd.dtimestr));
		if ( dt2 == 0. )
			error("\n\nstartrun_ParamStat: dtime : dt2 must be finite\n");
	}
	if (GetParamStat("dtout") & ARGPARAM) {
		cmd.dtoutstr = GetParam("dtout");
		gd.dtout = (sscanf(cmd.dtoutstr, "%lf/%lf", &dt1, &dt2) == 2 ?
			dt1/dt2 : atof(cmd.dtoutstr));
		if ( dt2 == 0. )
			error("\n\nstartrun_ParamStat: dtout : dt2 must be finite\n");
	}
	if (GetParamStat("dtoutinfo") & ARGPARAM) {
		cmd.dtoutstr = GetParam("dtoutinfo");
		gd.dtoutinfo = (sscanf(cmd.dtoutinfostr,"%lf/%lf",&dt1,&dt2) == 2 ?
			dt1/dt2 : atof(cmd.dtoutinfostr));
		if ( dt2 == 0. )
			error("\n\nstartrun_ParamStat: dtoutinfo : dt2 must be finite\n");
	}
}

local void CheckParameters(void)
{
	if (cmd.nbody < 1)
		error("CheckParameters: absurd value for nbody\n");
	if (cmd.stepState < 1)
		error("CheckParameters: absurd value for stepState\n");
	if (gd.dtime == 0.)
		error("\n\nCheckParameters: dtime must be != 0\n");
}


#define MFRAC  0.999                           

local void testdata(void)
{
    real rsc, vsc, r, v, x, y;
    vector rcm, vcm;
    bodyptr p;

    strcpy(gd.model_comment, "Plummer sphere created");

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

    rsc = (3 * PI) / 16;
    vsc = rsqrt(1.0 / rsc);
    CLRV(rcm);
    CLRV(vcm);
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
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
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        SUBV(Pos(p), Pos(p), rcm);
        SUBV(Vel(p), Vel(p), vcm);
    }
}

#undef MFRAC

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

	if((fd=fopen(fname,"r"))) {
		while(!feof(fd)) {
			fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
//            *buf2=(int)NULL;
// Removing the warning:
// warning: cast from pointer to integer of different size [-Wpointer-to-int-cast]
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
                                error("getbparam: %s=%s not bool\n",buf1,buf2);
                            }
						break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
					fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1); 
    }
  
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

#define FMTT	"%-35s%s\n"
#define FMTI	"%-35s%d\n"
#define FMTR	"%-35s%g\n"

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
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
            gd.headline3);
        fprintf(fdout,"%s\n%s\n",
		"%-------------------------------------------------------------------",
		"%");
        fprintf(fdout,FMTT,"icfile",cmd.icfile);
        fprintf(fdout,FMTT,"icfilefmt",cmd.icfilefmt);
        fprintf(fdout,FMTT,"snapout",cmd.snapoutfile);
        fprintf(fdout,FMTT,"snapoutfmt",cmd.snapoutfilefmt);
        fprintf(fdout,FMTT,"statefile",cmd.statefile);
        fprintf(fdout,FMTI,"stepState",cmd.stepState);
        fprintf(fdout,FMTT,"restorefile",cmd.restorefile);
        fprintf(fdout,FMTR,"eps",cmd.eps);
        fprintf(fdout,FMTT,"dtime",cmd.dtimestr);
        fprintf(fdout,FMTR,"tstop",cmd.tstop);
        fprintf(fdout,FMTT,"dtout",cmd.dtoutstr);
        fprintf(fdout,FMTT,"dtoutinfo",cmd.dtoutinfostr);
        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,FMTI,"nbody",cmd.nbody);
        fprintf(fdout,FMTI,"seed",cmd.seed);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

