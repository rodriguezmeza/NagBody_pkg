/*==============================================================================
    MODULE: nchi2.c				[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date:	January 2018
    Purpose:
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

local void min_Chi2(void (*funcs)(double, double [], double *, double [], int));

void MainLoop(void)
{
    real pv, qv, Chi2val;

// Reading data, setting model and initial parameters values
    set_model();
    if (strcmp(cmd.options,"test") == 0)
        test_model();
    if (strcmp(cmd.options,"testonly") == 0) {
        test_model();
        return;
    }

// Setting initial paramateres values
    grid_Chi2(gd.modelfun);
// Finding minimum value of Chi2
    min_Chi2(gd.modelfunwd);

// Writing useful info:
    Chi2val = gd.modelChi2(gd.modelfun, pd.params);
//    pv = P_value(Chi2val, nObs(pObs[gd.filecount])-pd.nparams);
//    qv = Q_value(Chi2val, nObs(pObs[gd.filecount])-pd.nparams);
    pv = P_value(Chi2val, pObsT.no-pd.nparams);
    qv = Q_value(Chi2val, pObsT.no-pd.nparams);
    fprintf(stdout,"\nP and Q values and their sum: %g %g %g\n",pv, qv, pv+qv);

    plotting_model_table(gd.modelfun,gd.fpfnamemodeltable);
//    gd.modelend();
    Model_end();
}

#define AMOEBA 1
#define NULLMETHOD 0
#define POWELL 2        // Direction set method
#define FRPRMN 3        // Conjugate gradient method
#define DFPMIN 5
#define AMEBSA 6
#define LEVMARQ 7

local void min_Chi2(void (*funcs)(double, double [], double *, double [], int))
{
    switch (gd.minmethod_int) {
        case AMOEBA:
//            min_Chi2_amoeba(gd.modelfunwd);
            min_Chi2_amoeba(funcs);
            break;
//
        case POWELL:
            min_Chi2_powell(gd.modelfunwd);
            break;
//
        case FRPRMN:
            min_Chi2_frprmn(gd.modelfunwd);
            break;
//
        case DFPMIN:
            min_Chi2_dfpmin(gd.modelfunwd);
            break;
//
        case AMEBSA:
            min_Chi2_amebsa(gd.modelfunwd);
            break;
//
        case LEVMARQ:
            min_Chi2_LevMarq(gd.modelfunwd);
            break;
//
        case NULLMETHOD:
            min_Chi2_amoeba(gd.modelfunwd);
            break;
//
        default:
            min_Chi2_amoeba(gd.modelfunwd);
            break;
    }
}

void minmethod_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"simplex") == 0) {
        *method_int = AMOEBA;
        strcpy(gd.minmethod_comment, "Simplex minimization method");
    }
//
    if (strcmp(method_str,"powell") == 0) {
        *method_int = POWELL;
        strcpy(gd.minmethod_comment, "Powell minimization method");
    }
//
    if (strcmp(method_str,"conjugategradient") == 0) {
        *method_int = FRPRMN;
        strcpy(gd.minmethod_comment, "Conjugate gradient minimization method");
    }
//
    if (strcmp(method_str,"quasinewton") == 0) {
        *method_int = DFPMIN;
        strcpy(gd.minmethod_comment, "Quasi-Newton minimization method");
    }
//
    if (strcmp(method_str,"amebsa") == 0) {
        *method_int = AMEBSA;
        strcpy(gd.minmethod_comment, "Simulated annealing + simplex minimization method");
    }
//
    if (strcmp(method_str,"levmarq") == 0) {
        *method_int = LEVMARQ;
        strcpy(gd.minmethod_comment, "Levenberg-Marquardt minimization method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.minmethod_comment,
               "given null minimization method ... running deafult (amoeba)");
        fprintf(stdout,"\n\tintegration: default integration method (amoeba)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = AMOEBA;
        strcpy(gd.minmethod_comment,
               "Unknown minimization method ... running deafult (amoeba)");
        fprintf(stdout,"\n\tMinimization: Unknown method... %s ",cmd.minMethod);
        fprintf(stdout,
                "\n\tRunning default minimization method (amoeba)...\n");
    }
}

#undef LEVMARQ
#undef AMOEBA
#undef POWELL
#undef FRPRMN
#undef DFPMIN
#undef AMEBSA
#undef NULLMETHOD
