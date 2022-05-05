/*==============================================================================
	MODULE: startrun.c				[galaxy_starscream]
	Written by: M.A. Rodriguez-Meza
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
	Copyright: (c) 2005-2014 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/general/getparam.h"
#include "globaldefs.h"
#include "protodefs.h"
#include <string.h>

local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);


void StartRun(string head0, string head1, string head2, string head3)
{
    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();
}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
    startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-galaxy_starscream"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
    cmd.ingal1 = GetParam("ingal1");
    cmd.ingal2 = GetParam("ingal2");

    cmd.snapoutfile = GetParam("snapout");
    cmd.snapoutfilefmt = GetParam("snapoutfmt");
	cmd.options = GetParam("options");

    cmd.a = GetdParam("a");
    cmd.b = GetdParam("b");
    cmd.g = GetdParam("g");
//
    cmd.radius = GetdParam("radius");
    cmd.p = GetdParam("p");
}

#undef parameter_null

#define logfile			"galaxy_starscream.log"

local void startrun_Common(void)
{
    strcpy(gd.mode,"w!");

	if(!(gd.outlog=fopen(logfile, gd.mode)))
		error("\nstart_Common: error opening file '%s' \n",logfile);

    CheckParameters();
}

#undef logfile

local void startrun_ParamStat(void)
{

	if (GetParamStat("options") & ARGPARAM)
		cmd.options = GetParam("options");
}

local void CheckParameters(void)
{
    if (strnull(cmd.ingal1))
        error("\nCheckParameters: you should suply an ingal1 file name\n\n");
//    if (strnull(cmd.ingal2))
//        error("\nCheckParameters: you should suply an ingal2 file name\n\n");

    if (strnull(cmd.snapoutfile))
        error("\nCheckParameters: you should suply an output file name\n\n");

    if (cmd.radius < 0)
        error("CheckParameters: absurd value for radius\n");
    if (cmd.p < 0)
        error("CheckParameters: absurd value for p\n");

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

    SPName(cmd.ingal1,"ingal1",1);
    SPName(cmd.ingal2,"ingal2",1);

    SPName(cmd.snapoutfile,"snapout",100);
	SPName(cmd.snapoutfilefmt,"snapoutfmt",100);
	SPName(cmd.options,"options",100);

    RPName(cmd.a,"a");
    RPName(cmd.b,"b");
    RPName(cmd.g,"g");
//
    RPName(cmd.radius,"radius");
    RPName(cmd.p,"p");

	if((fd=fopen(fname,"r"))) {
		while(!feof(fd)) {
			fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2=(int)NULL;
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

        fprintf(fdout,FMTT,"ingal1",cmd.ingal1);
        fprintf(fdout,FMTT,"ingal2",cmd.ingal2);

        fprintf(fdout,FMTT,"snapout",cmd.snapoutfile);
        fprintf(fdout,FMTT,"snapoutfmt",cmd.snapoutfilefmt);
        fprintf(fdout,FMTT,"options",cmd.options);
//
        fprintf(fdout,FMTR,"a",cmd.a);
        fprintf(fdout,FMTR,"b",cmd.b);
        fprintf(fdout,FMTR,"g",cmd.g);
//
        fprintf(fdout,FMTR,"radius",cmd.radius);
        fprintf(fdout,FMTR,"p",cmd.p);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

