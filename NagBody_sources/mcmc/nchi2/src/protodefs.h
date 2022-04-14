/*==============================================================================
    HEADER: protodefs.h				[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date: January 2018
    Purpose: Definitions of global prototypes
    Language: C
    Use: '#include "protodefs.h"
    Use in routines and functions:
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

#ifndef _protodefs_h
#define _protodefs_h

void minmethod_string_to_int(string,int *);

// IO routines
void InputObsDataTable(char *fname, int nfile);
void output(void);
global void plotting_model_table(real (*ymodel)(real, real *),char *fname);
global void setFilesDirs_log(void);
global void setFilesDirs(void);

void set_model(void);
void test_model(void);
global void set_model_Interp(real (*ymodel)(real, real *), real params[]);
global real model_Interp(real x);
global void Model_end(void);

void MainLoop(void);
void StartRun(string, string, string, string);
void StartOutput(void);
void EndRun(void);

// Min Chi2
global void min_Chi2_LevMarq(void (*funcs)(double, double [],
                                            double *, double [], int));
global void grid_Chi2(real (*ymodel)(real, real *));
global real Chi2(real (*ymodel)(real, real *));
global real P_value(real Chi2, real nu);
global real Q_value(real Chi2, real nu);

void amoeba_drv(void);
void powell_drv(void);
void frprmn_drv(void);
void dfpmin_drv(void);
void amebsa_drv(void);

global void min_Chi2_amoeba(void (*funcs)(double, double [], double *, double [], int));
global void min_Chi2_powell(void (*funcs)(double, double [], double *, double [], int));
global void min_Chi2_frprmn(void (*funcs)(double, double [], double *, double [], int));
global void min_Chi2_dfpmin(void (*funcs)(double, double [], double *, double [], int));
global void min_Chi2_amebsa(void (*funcs)(double, double [], double *, double [], int));

//

#endif // ! _protodefs_h
