/*==============================================================================
	MODULE: startrun.c			[model]
	Written by: M.A. Rodr'guez-Meza
	Starting date: January 2005
	Purpose: Initialize model
	Language: C
	Use: 'startrun();'
	Routines and functions: testdata, ReadParameterFile, PrintParameterFile
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: July 23, 2007
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
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

local void SnapDataTransform(void);
local void SnapDataTransform_long(void);

local void CheckInOutFmt(void);
local int in_long_fmt;
local int out_long_fmt;

local long saveidum;

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

	StartOutput();
}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
	startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-model"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

#undef parameter_null

local void ReadParametersCmdline(void)
{
	cmd.in = GetParam("in");    
	cmd.infmt = GetParam("infmt");    
	cmd.out = GetParam("out");
	cmd.outfmt = GetParam("outfmt");

	cmd.options = GetParam("options");

	cmd.model = GetParam("model-type");
	cmd.nbody = GetiParam("nbody");

	cmd.Mtotal=GetdParam("Mtotal");
	cmd.Rmax=GetdParam("Rmax");
	cmd.vcmx=GetdParam("vcmx");
	cmd.vcmy=GetdParam("vcmy");
	cmd.vcmz=GetdParam("vcmz");
	cmd.cmx=GetdParam("cmx");
	cmd.cmy=GetdParam("cmy");
	cmd.cmz=GetdParam("cmz");
	cmd.absvel=GetdParam("absvel");
	cmd.omega0=GetdParam("omega0");
	cmd.a_p=GetdParam("a_p");
	cmd.m_p=GetdParam("m_p");
	cmd.SoundSpeed=GetdParam("SoundSpeed");
	cmd.a_spheroid=GetdParam("a_spheroid");
	cmd.b_spheroid=GetdParam("b_spheroid");
	cmd.c_spheroid=GetdParam("c_spheroid");
	cmd.factor=GetdParam("factor");

	idum=GetiParam("seed");
}

#define logfile			"models.log"

local void startrun_Common(void)
{
	char *pch;
	char inputfiletmp[100], inputfilefmttmp[100];
	int i;

    if(!(gd.outlog=fopen(logfile,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",logfile);
        exit(0);
    }

	gd.headerfmt = "snap-blj-ascii";	// headerfmt debe de inicializarse siempre!!

	if (! strnull(cmd.in))	{
//
		strcpy(inputfiletmp,cmd.in);
		fprintf (gd.outlog,"\nSplitting string \"%s\" in tokens:\n",
			inputfiletmp);
		gd.nfiles=0;
		pch = strtok(inputfiletmp," ,");
		while (pch != NULL) {
			gd.filenames[gd.nfiles] = (string) malloc(30);
			strcpy(gd.filenames[gd.nfiles],pch);
			++gd.nfiles;
			fprintf (gd.outlog,"%s\n",gd.filenames[gd.nfiles-1]);
			pch = strtok (NULL, " ,");
		}
		fprintf (gd.outlog,"num. of inputfiles %s =%d\n",cmd.in,gd.nfiles);
//
		if (!strnull(cmd.infmt)) {
			strcpy(inputfilefmttmp,cmd.infmt);
			fprintf (gd.outlog,"\nSplitting string \"%s\" in tokens:\n",
				inputfilefmttmp);
			gd.nfilefmts=0;
			pch = strtok(inputfilefmttmp," ,");
			while (pch != NULL) {
				gd.filenamefmts[gd.nfilefmts] = (string) malloc(30);
				strcpy(gd.filenamefmts[gd.nfilefmts],pch);
				++gd.nfilefmts;
				fprintf (gd.outlog,"%s\n",gd.filenamefmts[gd.nfilefmts-1]);
				pch = strtok (NULL, " ,");
			}
			fprintf (gd.outlog,"num. of inputfile formats %s =%d\n",
				cmd.infmt,gd.nfilefmts);
			if (gd.nfiles != gd.nfilefmts)
				error("\nGet_Parameters: number of file formats must be equal to number of inputfiles\n\n");
		} else {
			for (i=0; i<gd.nfiles; i++) {
				gd.filenamefmts[i] = (string) malloc(30);
				strcpy(gd.filenamefmts[i],"snap-ascii");
			}
			gd.nfilefmts = gd.nfiles;
		}
//
// Inputdata_gadget driver for SPECIAL Particle data structure to manipulate I/O
// N > 10^6 purpose...
		CheckInOutFmt();
		if (in_long_fmt==1) {
			InputData_long();
			SnapDataTransform_long();
		} else {
			InputData();
			SnapDataTransform();
		}

//		if (strcmp(cmd.outfmt,"snap-ascii-long") == 0) {
//			InputData_long();
//			SnapDataTransform_long();
//		} else {
//			if (strcmp(cmd.outfmt,"powmes-ascii-long") == 0) {
//				InputData_long();
//				SnapDataTransform_long();
//			} else {
//				InputData();
//				SnapDataTransform();
//			}
//		}

		CheckParameters();
	} else {
		CheckParameters();
		saveidum=idum;
		xsrandom(idum);
		testdata();
		gd.tnow = 0.0;
	}

	gd.tout = gd.tnow;
	return;
}

#undef logfile

local void startrun_ParamStat(void)
{
// INCLUIR LOS PARAMETROS QUE PUEDAN SER UTILES DE AGREGAR 
// DESPUES DE UN ARCHIVO DE PARAMETROS ...
}

local void CheckParameters(void)
{
	if (cmd.nbody < 1)                      
		error("\nCheckParameters: absurd value for nbody\n\n");
	if (cmd.absvel < 0)
		error("\nCheckParameters: absurd value for absvel\n\n");
	if (cmd.a_spheroid < 0)
		error("\nCheckParameters: absurd value for a_spheroid\n\n");
	if (cmd.b_spheroid < 0)
		error("\nCheckParameters: absurd value for b_spheroid\n\n");
	if (cmd.c_spheroid < 0)
		error("\nCheckParameters: absurd value for c_spheroid\n\n");

    if (strnull(cmd.out))
		error("\nCheckParameters: you should suply an output file name\n\n");

	if (gd.nfiles != 1 && strcmp(cmd.outfmt,"powmes-ascii-long") == 0)
		error("\nCheckParameters: with an output long-file nfiles must be 1\n\n");

	if (gd.nfiles != 1 && strcmp(cmd.infmt,"snap-ascii-long") == 0)
		error("\nCheckParameters: with an output long-file nfiles must be 1\n\n");

//	if (strcmp(cmd.outfmt,"powmes-ascii-long") == 0 && 
//		strcmp(cmd.infmt,"gadget11-ascii-long") != 0)
//		error("\nCheckParameters: you should suply an input long-file\n\n");

//	if (strcmp(cmd.outfmt,"snap-ascii-long") == 0 && 
//		strcmp(cmd.infmt,"gadget11-ascii-long") != 0)
//		error("\nCheckParameters: you should suply an input long-file\n\n");

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

	SPName(cmd.in,"in",1);
	SPName(cmd.infmt,"infmt",1);
	SPName(cmd.out,"out",1);
	SPName(cmd.outfmt,"outfmt",100);
 	IPName(idum,"seed");
	SPName(cmd.options,"options",100);
 	IPName(cmd.nbody,"nbody");
	SPName(cmd.model,"model",100);
 	RPName(cmd.Mtotal,"Mtotal");
	RPName(cmd.Rmax,"Rmax");
	RPName(cmd.vcmx,"vcmx");
	RPName(cmd.vcmy,"vcmy");
	RPName(cmd.vcmz,"vcmz");
	RPName(cmd.cmx,"cmx");
	RPName(cmd.cmy,"cmy");
	RPName(cmd.cmz,"cmz");
	RPName(cmd.absvel,"absvel");
	RPName(cmd.omega0,"omega0");
	RPName(cmd.a_p,"a_p");
	RPName(cmd.m_p,"m_p");
	RPName(cmd.SoundSpeed,"SoundSpeed");
	RPName(cmd.a_spheroid,"a_spheroid");
	RPName(cmd.b_spheroid,"b_spheroid");
	RPName(cmd.c_spheroid,"c_spheroid");
	RPName(cmd.factor,"factor");
    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
			fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
				continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2=(int) NULL;
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
                                error("getbparam: %s=%s not bool\n", buf1, buf2);
                            }
						break;
				}
            } else {
                fprintf(stdout,
                    "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
                    fname,buf1);
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

    sprintf(buf,"%s%s",fname,"-models");
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
        fprintf(fdout,FMTT,"in",cmd.in);
        fprintf(fdout,FMTT,"infmt",cmd.infmt);
        fprintf(fdout,FMTT,"out",cmd.out);
        fprintf(fdout,FMTT,"outfmt",cmd.outfmt);

        fprintf(fdout,FMTI,"seed",saveidum);

        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,FMTI,"nbody",cmd.nbody);
        fprintf(fdout,FMTT,"model-type",cmd.model);

        fprintf(fdout,FMTR,"Mtotal",cmd.Mtotal);
        fprintf(fdout,FMTR,"Rmax",cmd.Rmax);
        fprintf(fdout,FMTR,"vcmx",cmd.vcmx);
        fprintf(fdout,FMTR,"vcmy",cmd.vcmy);
        fprintf(fdout,FMTR,"vcmz",cmd.vcmz);
        fprintf(fdout,FMTR,"cmx",cmd.cmx);
        fprintf(fdout,FMTR,"cmy",cmd.cmy);
        fprintf(fdout,FMTR,"cmz",cmd.cmz);
        fprintf(fdout,FMTR,"absvel",cmd.absvel);
        fprintf(fdout,FMTR,"omega0",cmd.omega0);
        fprintf(fdout,FMTR,"a_p",cmd.a_p);
        fprintf(fdout,FMTR,"m_p",cmd.m_p);
        fprintf(fdout,FMTR,"SoundSpeed",cmd.SoundSpeed);
        fprintf(fdout,FMTR,"a_spheroid",cmd.a_spheroid);
        fprintf(fdout,FMTR,"b_spheroid",cmd.b_spheroid);
        fprintf(fdout,FMTR,"c_spheroid",cmd.c_spheroid);
        fprintf(fdout,FMTR,"factor",cmd.factor);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

local void SnapDataTransform(void)
{
	bodyptr p;
	int k;

	if (scanopt(cmd.options, "mass-transform")) {
		DO_BODY(p, bodytab, bodytab+cmd.nbody) {
			Mass(p) *= cmd.factor;
		}
	}
	if (scanopt(cmd.options, "pos-transform")) {
		DO_BODY(p, bodytab, bodytab+cmd.nbody) {
			for (k=0; k<NDIM; k++)
				Pos(p)[k] *= cmd.factor;
		}
	}
	if (scanopt(cmd.options, "vel-transform")) {
		DO_BODY(p, bodytab, bodytab+cmd.nbody) {
			for (k=0; k<NDIM; k++)
				Vel(p)[k] *= cmd.factor;
		}
	}
}


local void SnapDataTransform_long(void)
{
	bodyptr_long p;
	int k;
	
	if (scanopt(cmd.options, "mass-transform")) {
		DO_BODY(p, bodytab_long, bodytab_long+cmd.nbody) {
			Mass_long(p) *= cmd.factor;
		}
	}
	if (scanopt(cmd.options, "pos-transform")) {
		DO_BODY(p, bodytab_long, bodytab_long+cmd.nbody) {
			for (k=0; k<NDIM; k++)
				Pos_long(p)[k] *= cmd.factor;
		}
	}
	if (scanopt(cmd.options, "vel-transform")) {
		DO_BODY(p, bodytab_long, bodytab_long+cmd.nbody) {
			for (k=0; k<NDIM; k++)
				Vel_long(p)[k] *= cmd.factor;
		}
	}
}


local void CheckInOutFmt(void)
{
	in_long_fmt = 0;
	out_long_fmt = 0;

	if (strcmp(cmd.infmt,"snap-ascii-long") == 0)
		in_long_fmt = 1;
	if (strcmp(cmd.outfmt,"snap-ascii-long") == 0)
		out_long_fmt = 1;
	if (strcmp(cmd.infmt,"gadget11-ascii-long") == 0)
		in_long_fmt = 1;
	if (strcmp(cmd.outfmt,"gadget11-normal-body-ascii-long") == 0)
		out_long_fmt = 1;
	if (strcmp(cmd.infmt,"heitmann-ascii-long") == 0)
		in_long_fmt = 1;
	if (strcmp(cmd.outfmt,"powmes-ascii-long") == 0)
		out_long_fmt = 1;

//	if (strcmp(cmd.outfmt,"powmes-ascii-long") == 0 && 
//		strcmp(cmd.infmt,"gadget11-ascii-long") != 0)
//		error("\nCheckParameters: you should suply an input long-file\n\n");
	
//	if (strcmp(cmd.outfmt,"snap-ascii-long") == 0 && 
//		strcmp(cmd.infmt,"gadget11-ascii-long") != 0)
//		error("\nCheckParameters: you should suply an input long-file\n\n");

	if (in_long_fmt != out_long_fmt)
		error("\nCheckInOutFmt: you have to give infmt and ofmt of the same type\n\n");
}

