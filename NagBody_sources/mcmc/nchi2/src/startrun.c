/*==============================================================================
    MODULE: startrun.c				[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date: January 2018
    Purpose: routines to initialize the main code
    Language: C
    Use: 'StartRun();'
    Routines and functions:
    Modules, routines and external headers:
    Coments and notes:
    Info: Mario A. Rodriguez-Meza
        Depto. de Fisica, ININ
        Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
        e-mail: marioalberto.rodriguez@inin.gob.mx
        https://github.com/rodriguezmeza

    Mayor revisions:
    Copyright: (c) 2005-2020 Mar.  All Rights Reserved
================================================================================
    Legal matters:
    The author does not warrant that the program and routines it contains
    listed below are free from error or suitable for particular applications,
    and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);
local void startrun_restorefile(void);

local void scanbOption(string, bool *, int *, int, int, string);


void StartRun(string head0, string head1, string head2, string head3)
{
    real aTime;
    aTime = cputime();
    
    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

	StartOutput();
    fprintf(stdout,"\nTime to StartRun: %g\n\n",cputime()-aTime);
    fflush(stdout);
}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
    startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-nchi2"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
// Data options
    cmd.inputObsDataFile = GetParam("inputObsDataFile");
    cmd.usingcolumns = GetParam("usingcolumns");
    cmd.witherrors = GetParam("witherrors");
    cmd.errorstype = GetiParam("errorstype");
// Model options:
    cmd.model = GetParam("modeltype");
    cmd.suffixModel = GetParam("suffixModel");
    cmd.model_paramfile = GetParam("modelParamfile");
// Output options:
    cmd.xoutmin = GetParam("xoutmin");
    cmd.xoutmax = GetParam("xoutmax");
    cmd.Nx = GetParam("Nx");
// Minimization parameters:
    cmd.minMethod = GetParam("minMethod");
    cmd.gridN = GetiParam("gridN");
    cmd.gridNSimplex = GetiParam("gridNSimplex");
    cmd.ftolmin = GetdParam("ftolMinimization");
    cmd.ucparams = GetbParam("useCentralParamValues");
// Other options:
    cmd.epsq = GetdParam("epsQuad");
    cmd.options = GetParam("options");
}

#undef parameter_null

#define logfile			"nchi2.log"

local void startrun_Common(void)
{
    char *pch;
    int i, npcols;
    short flag;
    char *pusingcolumns[30];
    char inputfiletmp[200], usingcolumnstmp[100];

    real xtmp;
    int itmp;

    setFilesDirs_log();
    strcpy(gd.mode,"w");
    if(!(gd.outlog=fopen(gd.logfilePath, gd.mode)))
        error("\nstart_Common: error opening file '%s' \n",gd.logfilePath);

// Input options:
    if (strnull(cmd.inputObsDataFile)) {
        error("\nstartrun_Common: no inputObsDataFile was given...\n");
    } else {
        strcpy(inputfiletmp,cmd.inputObsDataFile);
        fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",inputfiletmp);
        gd.nfiles=0;
        pch = strtok(inputfiletmp," ,");
        while (pch != NULL) {
            gd.filenames[gd.nfiles] = (string) malloc(100);
            strcpy(gd.filenames[gd.nfiles],pch);
            ++gd.nfiles;
            fprintf(stdout,"%s\n",gd.filenames[gd.nfiles-1]);
            pch = strtok (NULL, " ,");
        }
        fprintf(stdout,"num. of files in inputObsDataFile %s =%d\n",cmd.inputObsDataFile,gd.nfiles);
    }

    scanbOption(cmd.witherrors, gd.errorbars, &gd.nwitherrorbars, gd.nfiles, 1,
        "witherrors");

    strcpy(usingcolumnstmp,cmd.usingcolumns);
    fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",usingcolumnstmp);
    npcols=0;
    pch = strtok(usingcolumnstmp," ,");
    while (pch != NULL) {
        pusingcolumns[npcols] = (string) malloc(10);
        strcpy(pusingcolumns[npcols],pch);
        ++npcols;
        fprintf(stdout,"%s\n",pusingcolumns[npcols-1]);
        pch = strtok (NULL, " ,");
    }
    fprintf(stdout,
        "num. of sets of columns in usingcolumns %s =%d\n",cmd.usingcolumns,npcols);

    if (npcols != gd.nfiles)
        error("\nStartRunParameterfile: nusingcolumns must be equal to number of files\n\n");

      if (strnull(cmd.usingcolumns))
          error("startrun_Common: no observed data colummns are selected\n");

      gd.vcol1 = (int *) allocate(npcols*sizeof(int));
      gd.vcol2 = (int *) allocate(npcols*sizeof(int));

      flag = 1;
      for (i=0; i<npcols; i++) {
          if (gd.errorbars[i] && flag==1) {
              flag=0;
              if (cmd.errorstype==2) {
                  gd.vcol3 = (int *) allocate(npcols*sizeof(int));
                  gd.vcol4 = (int *) allocate(npcols*sizeof(int));
              } else {
                  if (cmd.errorstype==1)
                      gd.vcol3 = (int *) allocate(npcols*sizeof(int));
              }
          }
      }

    for (i=0; i<npcols; i++) {
    if (!strnull(cmd.witherrors) && gd.errorbars[i]==1) {
        switch (cmd.errorstype){
            case 2:
                if (!(sscanf(pusingcolumns[i], "%d:%d:%d:%d",
                             &gd.vcol1[i], &gd.vcol2[i], &gd.vcol3[i], &gd.vcol4[i]) == 4))
                    error("\nStartRunCmdline: usingcolumns must be in the form c1:c2:c3:c4\n\n");
                fprintf(stdout,"\nObserved data columns to read: %d %d %d %d\n",
                    gd.vcol1[i],gd.vcol2[i],gd.vcol3[i],gd.vcol4[i]);
                break;
            case 1:
                if (!(sscanf(pusingcolumns[i], "%d:%d:%d",
                             &gd.vcol1[i], &gd.vcol2[i], &gd.vcol3[i]) == 3))
                    error("\nStartRunCmdline: usingcolumns must be in the form c1:c2:c3\n\n");
                fprintf(stdout,"\nObserved data columns to read: %d %d %d\n",
                    gd.vcol1[i],gd.vcol2[i],gd.vcol3[i]);
                break;
            case 0:
                if (!(sscanf(pusingcolumns[i], "%d:%d", &gd.vcol1[i], &gd.vcol2[i]) == 2))
                    error("\nStartRunCmdline: usingcolumns must be in the form c1:c2\n\n");
                fprintf(stdout,"\nObserved data columns to read: %d %d\n",
                    gd.vcol1[i],gd.vcol2[i]);
                break;
            default: error("\nInputObsDataTable: errortype = %d is absurd\n\n",cmd.errorstype);
        }
    }
    }

    pObs = (ObsDataTableptr) allocate(gd.nfiles*sizeof(ObsDataTable));

// Output options:
    if (!strnull(cmd.xoutmin)) {
        gd.xoutmin = (sscanf(cmd.xoutmin, "%lf", &xtmp) == 1 ?
                xtmp : atof(cmd.xoutmin));
    }
    if (!strnull(cmd.xoutmax)) {
        gd.xoutmax = (sscanf(cmd.xoutmax, "%lf", &xtmp) == 1 ?
                xtmp : atof(cmd.xoutmax));
    }
    if (!strnull(cmd.Nx)) {
        gd.Nx = (sscanf(cmd.Nx, "%ld", &itmp) == 1 ?
                itmp : atoi(cmd.Nx));
    }

    CheckParameters();
    minmethod_string_to_int(cmd.minMethod, &gd.minmethod_int);
//
    setFilesDirs();

    if ( (strcmp(cmd.model,"USER") == 0) || (strcmp(cmd.model,"user") == 0)) {
        if (!strnull(cmd.model_paramfile)) {
            fprintf(stdout,"\n\nUser model :: using parameter file: %s\n",
                    cmd.model_paramfile);
            ReadModelParameterFile_user(cmd.model_paramfile);
            PrintModelParameterFile_user(cmd.model_paramfile);
            CheckParameters_model_user();
        }
    } else {
        if (!strnull(cmd.model_paramfile)) {
            fprintf(stdout,"\n\nModel :: using parameter file: %s\n",
                    cmd.model_paramfile);
            ReadModelParameterFile(cmd.model_paramfile);
            PrintModelParameterFile(cmd.model_paramfile);
            CheckParameters_model();
        }
    }

    set_model();
}

#undef logfile

local void startrun_ParamStat(void)
{
    real xtmp;
    int itmp;

// Model:
    if (GetParamStat("suffixModel") & ARGPARAM)
        cmd.suffixModel = GetParam("suffixModel");
// Output options:
    if (GetParamStat("xoutmin") & ARGPARAM) {
        cmd.xoutmin = GetParam("xoutmin");
        if (!strnull(cmd.xoutmin)) {
            gd.xoutmin = (sscanf(cmd.xoutmin, "%lf", &xtmp) == 1 ?
                    xtmp : atof(cmd.xoutmin));
        }
    }

// Minimization parameters:
    if (GetParamStat("minMethod") & ARGPARAM) {
        cmd.minMethod = GetParam("minMethod");
        fprintf(gd.outlog,"\n\nrunning instead %s minimization method ...\n",
                cmd.minMethod);
    }
    if (GetParamStat("gridN") & ARGPARAM)
        cmd.gridN = GetiParam("gridN");
    if (GetParamStat("gridNSimplex") & ARGPARAM)
        cmd.gridNSimplex = GetiParam("gridNSimplex");
    if (GetParamStat("ftolMinimization") & ARGPARAM)
        cmd.ftolmin = GetdParam("ftolMinimization");
    if (GetParamStat("useCentralParamValues") & ARGPARAM)
        cmd.ucparams = GetbParam("useCentralParamValues");

// Other options
    if (GetParamStat("epsQuad") & ARGPARAM)
        cmd.epsq = GetdParam("epsQuad");
    if (GetParamStat("options") & ARGPARAM)
		cmd.options = GetParam("options");
}

local void CheckParameters(void)
{
    if (!strnull(cmd.Nx)) {
        if (gd.Nx < 2)
            error("CheckParameters: absurd value for Nx: %d\n", gd.Nx);
    }
// Minimization parameters:
    if (cmd.gridN < 3)
            error("CheckParameters: absurd value for gridN: %d... must be greater than 3.\n", cmd.gridN);
    if (cmd.gridNSimplex < 10)
            error("CheckParameters: absurd value for gridNSimplex: %d... must be greater than 10.\n", cmd.gridNSimplex);
    if (cmd.ftolmin < 0)
            error("CheckParameters: absurd value for ftolMinimization: %g\n",
                  cmd.ftolmin);
// Other options:
    if (cmd.epsq < 0 || cmd.epsq > 1)
            error("CheckParameters: absurd value for epsQuad: %g\n",
                  cmd.epsq);
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

// Data options
    SPName(cmd.inputObsDataFile,"inputObsDataFile",200);
    SPName(cmd.usingcolumns,"usingcolumns",100);
    SPName(cmd.witherrors,"witherrors",100);
    IPName(cmd.errorstype,"errorstype");
// Model:
    SPName(cmd.model,"modeltype",100);
    SPName(cmd.suffixModel,"suffixModel",100);
    SPName(cmd.model_paramfile,"modelParamfile",100);
// Output options:
    SPName(cmd.xoutmin,"xoutmin",100);
    SPName(cmd.xoutmax,"xoutmax",100);
    SPName(cmd.Nx,"Nx",100);
// Minimization parameters:
    SPName(cmd.minMethod,"minMethod",100);
	SPName(cmd.options,"options",100);
    IPName(cmd.gridN,"gridN");
    IPName(cmd.gridNSimplex,"gridNSimplex");
    RPName(cmd.ftolmin,"ftolMinimization");
    BPName(cmd.ucparams,"useCentralParamValues");
// Other options:
    RPName(cmd.epsq,"epsQuad");

	if((fd=fopen(fname,"r"))) {
		while(!feof(fd)) {
			fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
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
    sprintf(buf,"%s/%s%s",gd.tmpDir,fname,"-usedvalues");
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
// Data options
        fprintf(fdout,FMTT,"inputObsDataFile",cmd.inputObsDataFile);
        fprintf(fdout,FMTT,"usingcolumns",cmd.usingcolumns);
        fprintf(fdout,FMTT,"witherrors",cmd.witherrors);
        fprintf(fdout,FMTI,"errorstype",cmd.errorstype);
// Model:
        fprintf(fdout,FMTT,"modeltype",cmd.model);
        fprintf(fdout,FMTT,"suffixModel",cmd.suffixModel);
        fprintf(fdout,FMTT,"modelParamfile",cmd.model_paramfile);
// Output options:
        fprintf(fdout,FMTT,"xoutmin",cmd.xoutmin);
        fprintf(fdout,FMTT,"xoutmax",cmd.xoutmax);
        fprintf(fdout,FMTT,"Nx",cmd.Nx);
// Minimization parameters:
        fprintf(fdout,FMTT,"minMethod",cmd.minMethod);
        fprintf(fdout,FMTI,"gridN",cmd.gridN);
        fprintf(fdout,FMTI,"gridNSimplex",cmd.gridNSimplex);
        fprintf(fdout,FMTR,"ftolMinimization",cmd.ftolmin);
        fprintf(fdout,FMTT,"useCentralParamValues",cmd.ucparams ? "true" : "false");
// Other options:
        fprintf(fdout,FMTR,"epsQuad",cmd.epsq);
        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR


local void scanbOption(string optionstr, bool *option, int *noption,
    int nfiles, int flag, string message)
{
    char *pch;
    char *poptionstr[30],  optiontmp[100];
    int i;

    fprintf(stdout,"\nProcessing '%s' option:\n", message);

    if (!strnull(optionstr)) {
        strcpy(optiontmp,optionstr);
        fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",optiontmp);
        *noption=0;
        pch = strtok(optiontmp," ,");
        while (pch != NULL) {
            poptionstr[*noption] = (string) malloc(10);
            strcpy(poptionstr[*noption],pch);
            ++(*noption);
            fprintf(stdout,"%s\n",poptionstr[*noption-1]);
            pch = strtok (NULL, " ,");
        }
        fprintf(stdout,"num. of tokens in option %s =%d\n",
            optionstr,*noption);

        if (flag == 0)
            if (*noption != nfiles)
                error("\nscanOption: noption = %d must be equal to number of files\n\n",*noption);
        if (*noption > MAXLINES)
            error("\nscanOption: noption = %d must be less than the maximum num. of lines\n\n",*noption);

        for (i=0; i<*noption; i++) {
            if (strchr("tTyY1", *poptionstr[i]) != NULL) {
                option[i]=TRUE;
            } else {
                if (strchr("fFnN0", *poptionstr[i]) != NULL)
                    option[i]=FALSE;
                else
                    error("\nscanOption: not bool in %s",poptionstr[i]);
            }
            fprintf(stdout,"option: %s\n", option[i] ? "true" : "false");
        }

        fprintf(stdout,"\nnoptions, nfiles: %d %d\n",*noption,nfiles);
        if (flag == 1) {
            if (*noption > nfiles)
                error("\nscanOption: noption = %d must be less or equal to number of files\n\n",*noption);
            else {
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[0];
                    fprintf(stdout,"option: %s\n", option[i] ? "true" : "false");
                }
                for (i=*noption; i<nfiles; i++) {
                    option[i]=option[0];
                    fprintf(stdout,"option: %d\n",option[i]);
                }
            }
        }
    } else {
        for (i=0; i<nfiles; i++) {
            option[i]=FALSE;
            fprintf(stdout,"option: %s\n", option[i] ? "true" : "false");
        }
    }
}
