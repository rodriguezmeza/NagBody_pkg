/*==============================================================================
    MODULE: nchi2_io.c		[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date:	January 2018
    Purpose: Routines to drive input and output data
    Language: C
    Use:
    Routines and functions:
    External modules, routines and headers:
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

#include "globaldefs.h"

local void outputdata(void);

void InputObsDataTable(char *fname, int nfile)
{
    char fpfnametmp[100];
    stream outstr;
    int i, nObsData;
//
    fprintf(gd.outlog,"\n\nReading observed data from file %s...\n",
            cmd.inputObsDataFile);

    if (!strnull(cmd.witherrors) && gd.errorbars[nfile]==1) {
        switch (cmd.errorstype){

            case 2:
                fprintf(stdout,"Reading four columns (errt=%d)...\n",cmd.errorstype);
                inout_InputData_4c(fname, gd.vcol1[nfile],
                                   gd.vcol2[nfile],
                                   gd.vcol3[nfile], gd.vcol4[nfile], &nObsData);
                break;
            case 1:
                fprintf(stdout,"Reading three columns (errt=%d)...\n",
                        cmd.errorstype);
                inout_InputData_3c(fname, gd.vcol1[nfile],
                                   gd.vcol2[nfile],
                                   gd.vcol3[nfile], &nObsData);
                break;
            case 0:
                    fprintf(stdout,"Reading two columns (errt=%d)...\n",cmd.errorstype);
                    inout_InputData(fname, gd.vcol1[nfile],
                                    gd.vcol2[nfile], &nObsData);
                break;

            default: error("\nUnknown errors type %s\n\n",cmd.errorstype);
        }
    }

    printf("\nnObsData = %d in file: %s\n",nObsData,fname);
//

    if (nObsData < 1)
        error("\n\nInputObsDataTable: nObsData = %d is absurd\n\n", nObsData);

    pObs[nfile] = (ObsDataTableptr) allocate(sizeof(ObsDataTable));
    nObs(pObs[nfile]) = nObsData;

    xObs(pObs[nfile])=dvector(1,nObsData);
    yObs(pObs[nfile])=dvector(1,nObsData);
    sigmaObs(pObs[nfile])=dvector(1,nObsData);
    sigmamObs(pObs[nfile])=dvector(1,nObsData);

    fprintf(gd.outlog,"nObsData for file %d : %d\n", nfile, nObs(pObs));

        for (i=1; i<=nObs(pObs[nfile]); i++) {
        xObs(pObs[nfile])[i] = inout_xval[i-1];
        yObs(pObs[nfile])[i] = inout_yval[i-1];
        if (cmd.errorstype==0) {
            sigmaObs(pObs[nfile])[i] = 1.0;
        } else {
            if (cmd.errorstype==1) {
                sigmaObs(pObs[nfile])[i] = inout_zval[i-1];
            } else {
                sigmaObs(pObs[nfile])[i] = inout_zval[i-1];
                sigmamObs(pObs[nfile])[i] = inout_wval[i-1];
            }
        }
    }

    sprintf(fpfnametmp,"%s/obsdatatabletmp_%d.dat",gd.tmpDir,nfile);
    outstr = stropen(fpfnametmp,"w!");
    if (cmd.errorstype==1 || cmd.errorstype==0) {
        for (i=1; i<=nObs(pObs[nfile]); i++) {
            fprintf(outstr,"%g %g %g\n",
                xObs(pObs[nfile])[i],yObs(pObs[nfile])[i],sigmaObs(pObs[nfile])[i]);
        }
    } else
        for (i=1; i<=nObs(pObs[nfile]); i++)
            fprintf(outstr,"%g %g %g %g\n",
                xObs(pObs[nfile])[i],yObs(pObs[nfile])[i],
                    sigmaObs(pObs[nfile])[i],sigmamObs(pObs[nfile])[i]);
        
    fclose(outstr);
}

void StartOutput(void)
{
    
    fprintf(stdout,"\n  \t -- %s --\n", gd.model_comment);
    //
    fprintf(stdout,"  \t -- %s --\n\n", gd.minmethod_comment);
    if (! strnull(cmd.options))
        fprintf(stdout,"\n\toptions: %s\n", cmd.options);
    
}

void output(void)
{
// Dummy
}


#define SOLFMT \
"%g %g\n"

local void outputdata(void)
{
    // Dummy
    fprintf(stdout,"\nIn outputdata");
}


global void plotting_model_table(real (*ymodel)(real, real *),char *fname)
{
    stream outstr;
    int i, N;
    real x1, x2, dx, x, y;

    outstr = stropen(fname,"w!");

    if (!strnull(cmd.xoutmin)) {
        x1 = gd.xoutmin;
    } else {
        x1 = xObs(pObs[gd.filecount])[1];
    }

    if (!strnull(cmd.xoutmax)) {
        x2 = gd.xoutmax;
    } else {
        x2 = xObs(pObs[gd.filecount])[nObs(pObs[gd.filecount])];
    }

    if (!strnull(cmd.Nx)) {
        N = gd.Nx;
    } else {
        N = nObs(pObs[gd.filecount]);
    }

    dx = (x2-x1)/((real)(N-1));
    for (i=1; i<=N; i++) {
        x = x1 + dx*((real)(i-1));
        y = ymodel(x,pd.params);
        fprintf(outstr,SOLFMT,x,y);
    }

    fclose(outstr);
}
#undef SOLFMT


// I/O directories:
global void setFilesDirs_log(void)
{
    char buf[200];
    
    sprintf(gd.tmpDir,"tmp");
    
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.tmpDir,gd.tmpDir);
    system(buf);

    sprintf(gd.logfilePath,"%s/nchi2%s.log",gd.tmpDir,cmd.suffixModel);
}

global void setFilesDirs(void)
{
    char buf[200];
    
    sprintf(gd.outDir,"Output");
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.outDir,gd.outDir);
    fprintf(gd.outlog,"system: %s\n",buf);
    system(buf);

    sprintf(gd.fpfnamemodeltest,"%s/%s%s_test.dat",
            gd.tmpDir,cmd.model,cmd.suffixModel);
    sprintf(gd.fpfnamemodeltable,"%s/%s%s.dat",
            gd.outDir,cmd.model,cmd.suffixModel);
}

void EndRun(void)
{
	fclose(gd.outlog);
	printf("\nFinal CPU time : %lf\n\n", cputime() - gd.cpuinit);
}


