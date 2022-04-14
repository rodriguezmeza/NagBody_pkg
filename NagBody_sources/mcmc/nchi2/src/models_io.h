/*==============================================================================
	HEADER: models_io.h			[nchi2]
	Written by: M.A. Rodriguez-Meza
	Starting date: January 2018
	Purpose: proto definitios of some model routines
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main)
	External headers: None
	Comments and notes:
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		https://github.com/rodriguezmeza

	Major revisions:
	Copyright: (c) 2005-2020 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _models_io_h
#define _models_io_h

//=============================================================
// Begin: Model reading and writing parameters

global void ReadModelParameterFile(char *fname)
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
    
    READPARAMS();
//
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

#define FMTT    "%-35s%s\n"
#define FMTI    "%-35s%d\n"
#define FMTR    "%-35s%g\n"
global void PrintModelParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];

    sprintf(buf,"%s/%s%s%s",gd.tmpDir,fname,cmd.suffixModel,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
                "%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Model parameters input file for:",gd.headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
                gd.headline3);
        fprintf(fdout,"%s\n%s\n",
                "%-------------------------------------------------------------------",
                "%");
//
        WRITEPARAMS();
//
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}
#undef FMTT
#undef FMTI
#undef FMTR

global void CheckParameters_model(void)
{
        if (md.nparams < 1)
            error("CheckParameters_model: absurd value for nparams: %d\n",
                  md.nparams);

    if (md.p1max <= md.p1min)
            error("CheckParameters_model: absurd values for p1min and p1max: %e %e... must be greater than 3.\n",
                  md.p1min, md.p1max);
}

// End: USER model reading and writing parameters
//=============================================================

#endif // !_models_io_h
