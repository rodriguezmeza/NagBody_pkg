/*==============================================================================
	HEADER: models.h			[nchi2]
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

#ifndef _models_h
#define _models_h

// Additional auxiliary routines or general model independent

global void model_wdyda(double x, double a[], double *y, double dyda[], int na);
global void model_wd1yda(double x, double a[], double *y, double d1yda[], int na);
global void model_wd2yda(double x, double a[], double **d2yda, int na);
global void model_wd2yda_chi2(double a[], double **d2yda, int na);

// User model :: reading/printing parameters
// Public interfaces but exclusive of user model:
global void ReadModelParameterFile_user(char *);
global void PrintModelParameterFile_user(char *);
global void CheckParameters_model_user(void);
//

// Other models :: reading/printing parameters
// Public interfaces but exclusive of the other models:
global void ReadModelParameterFile(char *);
global void PrintModelParameterFile(char *);
global void CheckParameters_model(void);
//

#endif // !_models_h
