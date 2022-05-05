/*==============================================================================
	MODULE: startrun.c				[galaxy_hernquist]
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
//#include <string.h>

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

#define parameter_null	"parameters_null-galaxy_hernquist"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
    cmd.snapoutfile = GetParam("snapout");
    cmd.snapoutfilefmt = GetParam("snapoutfmt");
	cmd.options = GetParam("options");
	cmd.seed = GetiParam("seed");
    cmd.Ngas = GetiParam("ngas");
    cmd.Nhalo = GetiParam("nhalo");
    cmd.Ndisk = GetiParam("ndisk");
    cmd.Nbulge = GetiParam("nbulge");
    cmd.Mgas = GetdParam("mgas");
    cmd.Mhalo = GetdParam("mhalo");
    cmd.Mdisk = GetdParam("mdisk");
    cmd.Mbulge = GetdParam("mbulge");
//
    cmd.masscut = GetdParam("masscut");
//
    cmd.ag = GetdParam("ag");
    cmd.ab = GetdParam("ab");
    cmd.ad = GetdParam("ad");
    cmd.ah = GetdParam("ah");
//
    cmd.zg = GetdParam("zg");
    cmd.zd = GetdParam("zd");
    cmd.gammah = GetdParam("gammah");
    cmd.gammab = GetdParam("gammab");
//
    cmd.rmaxg = GetdParam("rmaxg");
//    cmd.rmaxb = GetdParam("rmaxb");
    cmd.rmaxd = GetdParam("rmaxd");
//    cmd.rmaxh = GetdParam("rmaxh");
//
}

#undef parameter_null

#define logfile			"galaxy_hernquist.log"

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
    if (cmd.Ngas < 0)
        error("CheckParameters: absurd value for Ngas\n");
    if (cmd.Nhalo < 1)
        error("CheckParameters: absurd value for Nhalo\n");
    if (cmd.Ndisk < 1)
        error("CheckParameters: absurd value for Ndisk\n");
    if (cmd.Nbulge < 1)
        error("CheckParameters: absurd value for Nbulge\n");
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

	SPName(cmd.snapoutfile,"snapout",100);
	SPName(cmd.snapoutfilefmt,"snapoutfmt",100);
	SPName(cmd.options,"options",100);
    IPName(cmd.Ngas,"ngas");
    IPName(cmd.Nhalo,"nhalo");
    IPName(cmd.Ndisk,"ndisk");
    IPName(cmd.Nbulge,"nbulge");
    RPName(cmd.Mgas,"mgas");
    RPName(cmd.Mhalo,"mhalo");
    RPName(cmd.Mdisk,"mdisk");
    RPName(cmd.Mbulge,"mbulge");
//
    RPName(cmd.masscut,"masscut");
//
    RPName(cmd.ag,"ag");
    RPName(cmd.ab,"ab");
    RPName(cmd.ad,"ad");
    RPName(cmd.ah,"ah");
//
    RPName(cmd.zg,"zg");
    RPName(cmd.zd,"zd");
    RPName(cmd.gammah,"gammah");
    RPName(cmd.gammab,"gammab");
//
    RPName(cmd.rmaxg,"rmaxg");
//    RPName(cmd.rmaxb,"rmaxb");
    RPName(cmd.rmaxd,"rmaxd");
//    RPName(cmd.rmaxh,"rmaxh");
//
//
	IPName(cmd.seed,"seed");

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
        fprintf(fdout,FMTT,"snapout",cmd.snapoutfile);
        fprintf(fdout,FMTT,"snapoutfmt",cmd.snapoutfilefmt);
        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,FMTI,"ngas",cmd.Ngas);
        fprintf(fdout,FMTI,"nhalo",cmd.Nhalo);
        fprintf(fdout,FMTI,"ndisk",cmd.Ndisk);
        fprintf(fdout,FMTI,"nbulge",cmd.Nbulge);
        fprintf(fdout,FMTR,"mgas",cmd.Mgas);
        fprintf(fdout,FMTR,"mhalo",cmd.Mhalo);
        fprintf(fdout,FMTR,"mdisk",cmd.Mdisk);
        fprintf(fdout,FMTR,"mbulge",cmd.Mbulge);
//
        fprintf(fdout,FMTR,"masscut",cmd.masscut);
//
        fprintf(fdout,FMTR,"ag",cmd.ag);
        fprintf(fdout,FMTR,"ab",cmd.ab);
        fprintf(fdout,FMTR,"ad",cmd.ad);
        fprintf(fdout,FMTR,"ah",cmd.ah);
//
        fprintf(fdout,FMTR,"zg",cmd.zg);
        fprintf(fdout,FMTR,"zd",cmd.zd);
        fprintf(fdout,FMTR,"gammah",cmd.gammah);
        fprintf(fdout,FMTR,"gammab",cmd.gammab);
//
        fprintf(fdout,FMTR,"rmaxg",cmd.rmaxg);
//        fprintf(fdout,FMTR,"rmaxb",cmd.rmaxb);
        fprintf(fdout,FMTR,"rmaxd",cmd.rmaxd);
//        fprintf(fdout,FMTR,"rmaxh",cmd.rmaxh);
//
        fprintf(fdout,FMTI,"seed",cmd.seed);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

