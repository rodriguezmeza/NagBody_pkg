/*==============================================================================
	MODULE: startrun.c			[analysis_grav]
	Written by: Mario A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Initialize datanaly_md
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

	Mayor revisions: January 2007;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
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

local void Cmd_to_Global(void);
local void ComputeParameters(void);

#define parameter_null	"parameters_null"
#define logfile			"analysis_grav.log"

local void CheckInOutFmt(void);

void startrun(void)
{
	RestartFlag=0;					// De gadget ... se debe inicializar en 0...

	gd.x_autoscale = TRUE;
	gd.y_autoscale = TRUE;
	gd.z_autoscale = TRUE;

//    LicDriver(gd.headline0);

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

	CheckInOutFmt();

	Cmd_to_Global();
	ComputeParameters();
}

local void startrun_parameterfile(void)
{
	char *pch;
	char bodiesIDtmp[256];					// 8 x MAXBODIESID < 256
	char bodiesSetstmp[256];					// 8 x MAXBODIESSETS < 256
	char yrangeSetstmp[256];
	real dt1, dt2;
	int xrangeflag=1, yrangeflag=1, zrangeflag=1;
	int i;

    if(!(gd.outlog=fopen(logfile,"w")))
        error("error opening file '%s' \n",logfile);

	ReadParameterFile(cmd.paramfile);

	gd.headerfmt = "snap-tlj-ascii";

	if (!strnull(cmd.bodiesID)) {
		strcpy(bodiesIDtmp,cmd.bodiesID);
		fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",bodiesIDtmp);
		gd.nbodiesID=0;
		pch = strtok(bodiesIDtmp," ,");
		while (pch != NULL) {
			gd.bodyID[gd.nbodiesID] = atoi(pch);
			fprintf(stdout,"%d\n",gd.bodyID[gd.nbodiesID]);
			++gd.nbodiesID;
			pch = strtok (NULL, " ,");
		}
		fprintf(stdout,"num. of bodies ID in bodiesID %s =%d\n",
			cmd.bodiesID,gd.nbodiesID);
	} //else
//		if (scanopt(cmd.data_analysis_type, "bodies-anim"))
//			error("\n\nstartrun_parameterfile: analysis_type=bodies-anim %s\n",
//				"must be run with nbodiesID not null");
	if (gd.nbodiesID > MAXBODIESID)
		error("\n\nstartrun_parameterfile: nbodiesID must be <= MAXBODIESID\n\n");
	for (i=0; i<gd.nbodiesID; i++)
		if (gd.bodyID[i]<0)
			error("\n\nstartrun_parameterfile: bodyID < zero\n\n");
//
	if (!strnull(cmd.bodiesSets)) {
		strcpy(bodiesSetstmp,cmd.bodiesSets);
		fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",bodiesSetstmp);
		gd.nbodiesSets=0;
		pch = strtok(bodiesSetstmp," ,");
		while (pch != NULL) {
			if (!(sscanf(pch, "%d:%d", 
					&gd.bodyIDMin[gd.nbodiesSets], 
					&gd.bodyIDMax[gd.nbodiesSets]) == 2))
				error("\nstartrun_parameterfile: bodies set must be in the form d1:d2\n\n");
			fprintf(stdout,"%d %d\n",
				gd.bodyIDMin[gd.nbodiesSets],gd.bodyIDMax[gd.nbodiesSets]);
			++gd.nbodiesSets;
			pch = strtok (NULL, " ,");
		}
		fprintf(stdout,"num. of bodies sets in bodiesSets %s =%d\n",
			cmd.bodiesSets,gd.nbodiesSets);
	}
	if (gd.nbodiesSets > MAXBODIESSETS)
		error("\n\nstartrun_parameterfile: nbodiesSets must be <= MAXBODIESSETS\n\n");
	for (i=0; i<gd.nbodiesSets; i++) {
		if (gd.bodyIDMin[i]<0)
			error("\n\nstartrun_parameterfile: bodyIDMin < zero\n\n");
		if (gd.bodyIDMax[i]<0)
			error("\n\nstartrun_parameterfile: bodyIDMax < zero\n\n");
		if (gd.bodyIDMin[i]>gd.bodyIDMax[i])
			error("\n\nstartrun_parameterfile: bodyIDMin > bodyIDMax \n\n");
	}
//
	if (scanopt(cmd.data_analysis_type, "bodies-anim")) {
		if (strnull(cmd.bodiesID)) {
			if (strnull(cmd.bodiesSets)) {
				error("\n\nstartrun_parameterfile: analysis_type=bodies-anim %s\n\n", 
					"must be run with nbodiesID or nbodiesSets not null\n\n");
			}
		}
	}
//
	if (scanopt(cmd.data_analysis_type, "rhotheta-anim")) {
		if (strnull(cmd.bodiesSets)) {
			error("\n\nstartrun_parameterfile: analysis_type=rhotheta-anim %s\n\n", 
				"must be run with bodiesSets not null\n\n");
		}
		if (gd.nbodiesSets!=3) {
			error("\n\nstartrun_parameterfile: analysis_type=rhotheta-anim %s\n\n", 
				"must be run with bodiesSets equal to 3\n\n");
		}
	}
//
	if (scanopt(cmd.data_analysis_type, "locate-bodiesid")) {
		if (strnull(cmd.bodiesSets)) {
			error("\n\nstartrun_parameterfile: analysis_type=locate-bodiesid %s\n\n", 
				"must be run with bodiesSets not null\n\n");
		}
		if (gd.nbodiesSets!=3) {
			error("\n\nstartrun_parameterfile: analysis_type=locate-bodiesid %s\n\n", 
				"must be run with bodiesSets equal to 3\n\n");
		}
	}
//
		if (scanopt(cmd.data_analysis_type,"groups_catalog") && strnull(cmd.foffile))
			error("\n\ngroups_catalog datanaly_type option must be accompained with a fof file\n");

//		if (GetParamStat("stepEquil") & ARGPARAM)
//			stepEquil = GetiParam("stepEquil");

		if (GetParamStat("stepAvgRhoTheta") & ARGPARAM)
			cmd.stepAvgRhoTheta = GetiParam("stepAvgRhoTheta");         
		if (GetParamStat("sizeHistRhoTheta") & ARGPARAM)
			cmd.sizeHistRhoTheta = GetiParam("sizeHistRhoTheta");
		if (GetParamStat("RhoDeltaZ") & ARGPARAM)
			cmd.RhoDeltaZ = GetdParam("RhoDeltaZ");
		if (GetParamStat("RhoR") & ARGPARAM)
			cmd.RhoR = GetdParam("RhoR");
		if (GetParamStat("RhoDeltaR") & ARGPARAM)
			cmd.RhoDeltaR = GetdParam("RhoDeltaR");

		if (GetParamStat("stepAvgRho") & ARGPARAM)
			cmd.stepAvgRho = GetiParam("stepAvgRho");         
		if (GetParamStat("sizeHistRho") & ARGPARAM)
			cmd.sizeHistRho = GetiParam("sizeHistRho");

		if (GetParamStat("stepAvgVcR") & ARGPARAM)
			cmd.stepAvgVcR = GetiParam("stepAvgVcR");         
		if (GetParamStat("sizeHistVcR") & ARGPARAM)
			cmd.sizeHistVcR = GetiParam("sizeHistVcR");
		if (GetParamStat("rangeR") & ARGPARAM)
			cmd.rangeR = GetdParam("rangeR");

		if (cmd.reductionFac <= 0.0)
			error("\nstartrun_parameterfile: reductionFac must be positive : \n\n",
				cmd.reductionFac);

		if (GetParamStat("options") & ARGPARAM)
			cmd.options = GetParam("options");
		if (strnull(cmd.data_analysis_type))
			error("\nstartrun_parameterfile: run from command line as:\n%s",
				"\t'analysis_md analysis_type=analysis_type_name'\n\n");

		if (cmd.isnap<0 || cmd.fsnap<0)
			error("\n\nstartrun_parameterfile: isnap and fsnap must be >= 0\n");
		if (cmd.isnap > cmd.fsnap)
			error("\n\nstartrun_parameterfile: isnap must be < fsnap\n");

		if (strcmp(cmd.xrange,"autoscale"))
			gd.x_autoscale=FALSE;
		if (!(sscanf(cmd.xrange, "%lf:%lf", &gd.xmin, &gd.xmax) == 2))
			xrangeflag=0;
		if ( !gd.x_autoscale && xrangeflag==0) 
			error("\nstartrun_parameterfile: xrange must be autoscale or in the form xmin:xmax\n\n");
		if (xrangeflag != 0 && gd.xmin >= gd.xmax)
			error("\nstartrun_parameterfile: xmin must be < xmax : (xmin,xmax)=(%lf,%lf)\n\n",gd.xmin,gd.xmax);

//
	if (scanopt(cmd.data_analysis_type, "rhotheta-anim")) {
		if (strcmp(cmd.yrange,"autoscale")) {
			gd.y_autoscale=FALSE;
			yrangeflag=0;
			if (!strnull(cmd.yrange)) {
				strcpy(yrangeSetstmp,cmd.yrange);
				fprintf(stdout,"\nyrange : Splitting string \"%s\" in tokens:\n",yrangeSetstmp);
				gd.nyrangeSets=0;
				pch = strtok(yrangeSetstmp," ,");
				while (pch != NULL) {
					if (!(sscanf(pch, "%lf:%lf", 
						&gd.yrangeMin[gd.nyrangeSets], 
						&gd.yrangeMax[gd.nyrangeSets]) == 2))
						error("\nstartrun_parameterfile: yrange set must be in the form f1:f2\n\n");
					fprintf(stdout,"%g %g\n",
					gd.yrangeMin[gd.nyrangeSets],gd.yrangeMax[gd.nyrangeSets]);
					++gd.nyrangeSets;
					pch = strtok (NULL, " ,");
				}
				fprintf(stdout,"num. of yrange sets in yrange %s =%d\n",
					cmd.yrange,gd.nyrangeSets);
			}
			if (gd.nyrangeSets != MAXYRANGESETS)
				error("\n\nstartrun_parameterfile: yrangeSets must be = 6\n\n");
			for (i=0; i<MAXYRANGESETS; i++) {
				if (gd.yrangeMin[i]>gd.yrangeMax[i])
					error("\n\nstartrun_parameterfile: yrangeMin > yrangeMax \n\n");
			}
			yrangeflag=1;
		}
		if ( !gd.y_autoscale && yrangeflag==0) 
			error("\nstartrun_parameterfile: yrange must be autoscale or in the form of a set of ymin:ymax\n\n");

	} else {
		if (strcmp(cmd.yrange,"autoscale"))
			gd.y_autoscale=FALSE;
		if (!(sscanf(cmd.yrange, "%lf:%lf", &gd.ymin, &gd.ymax) == 2))
			yrangeflag=0;
		if ( !gd.y_autoscale && yrangeflag==0) 
			error("\nstartrun_parameterfile: yrange must be autoscale or in the form ymin:ymax\n\n");
		if (yrangeflag != 0 && gd.ymin >= gd.ymax)
			error("\nstartrun_parameterfile: ymin must be < ymax : (ymin,ymax)=(%lf,%lf)\n\n",gd.ymin,gd.ymax);
	}
//

/*
		if (strcmp(cmd.yrange,"autoscale"))
			gd.y_autoscale=FALSE;
		if (!(sscanf(cmd.yrange, "%lf:%lf", &gd.ymin, &gd.ymax) == 2))
			yrangeflag=0;
		if ( !gd.y_autoscale && yrangeflag==0) 
			error("\nGet_Parameters: yrange must be autoscale or in the form ymin:ymax\n\n");
		if (yrangeflag != 0 && gd.ymin >= gd.ymax)
			error("\nGet_Parameters: ymin must be < ymax : (ymin,ymax)=(%lf,%lf)\n\n",gd.ymin,gd.ymax);
*/
		if (strcmp(cmd.zrange,"autoscale"))
			gd.z_autoscale=FALSE;
		if (!(sscanf(cmd.zrange, "%lf:%lf", &gd.zmin, &gd.zmax) == 2))
			zrangeflag=0;
		if ( !gd.z_autoscale && zrangeflag==0) 
			error("\nstartrun_parameterfile: zrange must be autoscale or in the form zmin:zmax\n\n");
		if (zrangeflag != 0 && gd.zmin >= gd.zmax)
			error("\nstartrun_parameterfile: zmin must be < zmax : (zmin,zmax)=(%lf,%lf)\n\n",gd.zmin,gd.zmax);

		if (scanopt(cmd.data_analysis_type,"thermo-avg") && scanopt(cmd.options, "errorbar")) {
			if (!(sscanf(cmd.usingcolumns, "%d:%d:%d", 
				&gd.column1, &gd.column2, &gd.column3) == 3))
				error("\nstartrun_parameterfile: usingcolumns must be in the form c1:c2:c3\n\n");
		} else
			if (!(sscanf(cmd.usingcolumns, "%d:%d", &gd.column1, &gd.column2) == 2))
				error("\nstartrun_parameterfile: usingcolumns must be in the form c1:c2\n\n");

		if (!(sscanf(cmd.usingrows, "%d:%d", &gd.row1, &gd.row2) == 2))
			error("\nstartrun_parameterfile: usingrows must be in the form r1:r2\n\n");

        PrintParameterFile(cmd.paramfile);
}

local void startrun_cmdline(void)
{
	char *pch;
	char bodiesIDtmp[256];					// 8 x MAXBODIESID < 256
	char bodiesSetstmp[256];					// 8 x MAXBODIESSETS < 256
	char yrangeSetstmp[256];
	real dt1, dt2;
	int xrangeflag=1, yrangeflag=1, zrangeflag=1;
	int i;

    if(!(gd.outlog=fopen(logfile,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",logfile);
        exit(0);
    }

// Seran activados cuando se necesite calcular la fuerza entre particulas.
//	cmd.usequad = GetbParam("usequad");
//	cmd.theta = GetdParam("theta");
//	cmd.eps = GetdParam("eps");

    cmd.in = GetParam("in");
    cmd.infmt = GetParam("infmt");
    cmd.out = GetParam("out");
    cmd.outfmt = GetParam("outfmt");
	cmd.basedir = GetParam("basedir");
	cmd.isnap = GetiParam("isnap");
	cmd.fsnap = GetiParam("fsnap");
	if (cmd.isnap<0 || cmd.fsnap<0)
		error("\n\nstartrun: isnap and fsnap must be >= 0\n");
	if (cmd.isnap > cmd.fsnap)
		error("\n\nstartrun: isnap must be < fsnap\n");

	gd.headerfmt = "snap-tlj-ascii";

	cmd.stepAvgRhoTheta = GetiParam("stepAvgRhoTheta");
	cmd.sizeHistRhoTheta = GetiParam("sizeHistRhoTheta");
	cmd.RhoDeltaZ = GetdParam("RhoDeltaZ");
	cmd.RhoR = GetdParam("RhoR");
	cmd.RhoDeltaR = GetdParam("RhoDeltaR");
	cmd.RMax = GetdParam("RMax");

	cmd.xmin = GetdParam("xmin");
	cmd.xmax = GetdParam("xmax");
	cmd.ymin = GetdParam("ymin");
	cmd.ymax = GetdParam("ymax");
	cmd.zmin = GetdParam("zmin");
	cmd.zmax = GetdParam("zmax");

	cmd.ThetaMin = GetdParam("ThetaMin");
	cmd.ThetaMax = GetdParam("ThetaMax");

	cmd.stepAvgRho = GetiParam("stepAvgRho");
	cmd.sizeHistRho = GetiParam("sizeHistRho");

	cmd.stepAvgVcR = GetiParam("stepAvgVcR");
	cmd.sizeHistVcR = GetiParam("sizeHistVcR");
	cmd.rangeR = GetdParam("rangeR");

//	cmd.computeTransport = GetbParam("computeTransport");
//	cmd.stepAcf = GetiParam("stepAcf");
	cmd.limitAcfAv = GetiParam("limitAcfAv");
	cmd.nBuffAcf = GetiParam("nBuffAcf");
	cmd.nValAcf = GetiParam("nValAcf");

	cmd.options = GetParam("options");
//	nbody = GetiParam("nbody");         
//	if (nbody < 1)                      
//		error("startrun: absurd value for nbody\n");
	cmd.data_analysis_type = GetParam("analysis_type");
	if (strnull(cmd.data_analysis_type))
		error("\nrun from command line as:\n%s",
			"\t'analysis_md analysis_type=analysis_type_name'\n\n");
//
	cmd.bodiesID = GetParam("bodiesID");
	if (!strnull(cmd.bodiesID)) {
		strcpy(bodiesIDtmp,cmd.bodiesID);
		fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",bodiesIDtmp);
		gd.nbodiesID=0;
		pch = strtok(bodiesIDtmp," ,");
		while (pch != NULL) {
			gd.bodyID[gd.nbodiesID] = atoi(pch);
			fprintf(stdout,"%d\n",gd.bodyID[gd.nbodiesID]);
			++gd.nbodiesID;
			pch = strtok (NULL, " ,");
		}
		fprintf(stdout,"num. of bodies ID in bodiesID %s =%d\n",
			cmd.bodiesID,gd.nbodiesID);
	} //else
//		if (scanopt(cmd.data_analysis_type, "bodies-anim"))
//			error("\n\nstartrun_cmdline: analysis_type=bodies-anim must be run with nbodiesID not null\n\n");
	if (gd.nbodiesID > MAXBODIESID)
		error("\n\nstartrun_parameterfile: nbodiesID must be <= MAXBODIESID\n\n");
	for (i=0; i<gd.nbodiesID; i++)
		if (gd.bodyID[i]<0)
			error("\n\nstartrun_parameterfile: bodyID < zero\n\n");
//
	cmd.bodiesSets = GetParam("bodiesSets");
	if (!strnull(cmd.bodiesSets)) {
		strcpy(bodiesSetstmp,cmd.bodiesSets);
		fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",bodiesSetstmp);
		gd.nbodiesSets=0;
		pch = strtok(bodiesSetstmp," ,");
		while (pch != NULL) {
			if (!(sscanf(pch, "%d:%d", 
					&gd.bodyIDMin[gd.nbodiesSets], 
					&gd.bodyIDMax[gd.nbodiesSets]) == 2))
				error("\nstartrun_cmdline: bodies set must be in the form d1:d2\n\n");
			fprintf(stdout,"%d %d\n",
				gd.bodyIDMin[gd.nbodiesSets],gd.bodyIDMax[gd.nbodiesSets]);
			++gd.nbodiesSets;
			pch = strtok (NULL, " ,");
		}
		fprintf(stdout,"num. of bodies sets in bodiesSets %s =%d\n",
			cmd.bodiesSets,gd.nbodiesSets);
	}
	if (gd.nbodiesSets > MAXBODIESSETS)
		error("\n\nstartrun_cmdline: nbodiesSets must be <= MAXBODIESSETS\n\n");
	for (i=0; i<gd.nbodiesSets; i++) {
		if (gd.bodyIDMin[i]<0)
			error("\n\nstartrun_cmdline: bodyIDMin < zero\n\n");
		if (gd.bodyIDMax[i]<0)
			error("\n\nstartrun_cmdline: bodyIDMax < zero\n\n");
		if (gd.bodyIDMin[i]>gd.bodyIDMax[i])
			error("\n\nstartrun_cmdline: bodyIDMin > bodyIDMax \n\n");
	}
//
	if (scanopt(cmd.data_analysis_type, "bodies-anim")) {
		if (strnull(cmd.bodiesID)) {
			if (strnull(cmd.bodiesSets)) {
				error("\n\nstartrun_cmdline: analysis_type=bodies-anim %s\n\n", 
					"must be run with bodiesID or bodiesSets not null\n\n");
			}
		}
	}
//
	if (scanopt(cmd.data_analysis_type, "rhotheta-anim")) {
		if (strnull(cmd.bodiesSets)) {
			error("\n\nstartrun_cmdline: analysis_type=rhotheta-anim %s\n\n", 
				"must be run with bodiesSets not null\n\n");
		}
		if (gd.nbodiesSets!=3) {
			error("\n\nstartrun_cmdline: analysis_type=rhotheta-anim %s\n\n", 
				"must be run with bodiesSets equal to 3\n\n");
		}
	}
//
	if (scanopt(cmd.data_analysis_type, "locate-bodiesid")) {
		if (strnull(cmd.bodiesSets)) {
			error("\n\nstartrun_cmdline: analysis_type=locate-bodiesid %s\n\n", 
				"must be run with bodiesSets not null\n\n");
		}
		if (gd.nbodiesSets!=3) {
			error("\n\nstartrun_cmdline: analysis_type=locate-bodiesid %s\n\n", 
				"must be run with bodiesSets equal to 3\n\n");
		}
	}
//
	cmd.reductionFac = GetdParam("reductionFac");
	if (cmd.reductionFac <= 0.0)
		error("\nstartrun_parameterfile: reductionFac must be positive : \n\n",
			cmd.reductionFac);

	cmd.foffile = GetParam("foffile");
	if (scanopt(cmd.data_analysis_type,"groups_catalog") && strnull(cmd.foffile))
		error("\n\ngroups_catalog datanaly_type option must be accompained with a fof file\n");

	cmd.gravityConstant = GetdParam("gravityConstant");
//
	cmd.xrange = GetParam("xrange");
	if (strcmp(cmd.xrange,"autoscale"))
		gd.x_autoscale=FALSE;
	if (!(sscanf(cmd.xrange, "%lf:%lf", &gd.xmin, &gd.xmax) == 2))
		xrangeflag=0;
	if ( !gd.x_autoscale && xrangeflag==0) 
		error("\nGet_Parameters: xrange must be autoscale or in the form xmin:xmax\n\n");
	if (xrangeflag != 0 && gd.xmin >= gd.xmax)
		error("\nGet_Parameters: xmin must be < xmax : (xmin,xmax)=(%lf,%lf)\n\n",gd.xmin,gd.xmax);

//
	cmd.yrange = GetParam("yrange");
	if (scanopt(cmd.data_analysis_type, "rhotheta-anim") 
			|| scanopt(cmd.data_analysis_type, "rho-anim") ) {
		if (strcmp(cmd.yrange,"autoscale")) {
			gd.y_autoscale=FALSE;
			yrangeflag=0;
			if (!strnull(cmd.yrange)) {
				strcpy(yrangeSetstmp,cmd.yrange);
				fprintf(stdout,"\nyrange : Splitting string \"%s\" in tokens:\n",yrangeSetstmp);
				gd.nyrangeSets=0;
				pch = strtok(yrangeSetstmp," ,");
				while (pch != NULL) {
					if (!(sscanf(pch, "%lf:%lf", 
						&gd.yrangeMin[gd.nyrangeSets], 
						&gd.yrangeMax[gd.nyrangeSets]) == 2))
						error("\nstartrun_cmdline: yrange set must be in the form f1:f2\n\n");
					fprintf(stdout,"%g %g\n",
					gd.yrangeMin[gd.nyrangeSets],gd.yrangeMax[gd.nyrangeSets]);
					++gd.nyrangeSets;
					pch = strtok (NULL, " ,");
				}
				fprintf(stdout,"num. of yrange sets in yrange %s =%d\n",
					cmd.yrange,gd.nyrangeSets);
			}
			if (gd.nyrangeSets != MAXYRANGESETS 
				&& !scanopt(cmd.data_analysis_type, "rho-anim") )
				error("\n\nstartrun_cmdline: yrangeSets must be = 6\n\n");
			for (i=0; i<MAXYRANGESETS; i++) {
				if (gd.yrangeMin[i]>gd.yrangeMax[i])
					error("\n\nstartrun_cmdline: yrangeMin > yrangeMax \n\n");
			}
			yrangeflag=1;
		}
		if ( !gd.y_autoscale && yrangeflag==0) 
			error("\nstartrun_cmdline: yrange must be autoscale or in the form of a set of ymin:ymax\n\n");

//		if (!(sscanf(GetParam("yrange"), "%lf:%lf", &gd.ymin, &gd.ymax) == 2))
//			yrangeflag=0;
//		if ( !gd.y_autoscale && yrangeflag==0) 
//			error("\nGet_Parameters: yrange must be autoscale or in the form ymin:ymax\n\n");
//		if (yrangeflag != 0 && gd.ymin >= gd.ymax)
//			error("\nGet_Parameters: ymin must be < ymax : (ymin,ymax)=(%lf,%lf)\n\n",gd.ymin,gd.ymax);
	} else {
		if (strcmp(cmd.yrange,"autoscale"))
			gd.y_autoscale=FALSE;
		if (!(sscanf(GetParam("yrange"), "%lf:%lf", &gd.ymin, &gd.ymax) == 2))
			yrangeflag=0;
		if ( !gd.y_autoscale && yrangeflag==0) 
			error("\nGet_Parameters: yrange must be autoscale or in the form ymin:ymax\n\n");
		if (yrangeflag != 0 && gd.ymin >= gd.ymax)
			error("\nGet_Parameters: ymin must be < ymax : (ymin,ymax)=(%lf,%lf)\n\n",gd.ymin,gd.ymax);
	}
//

	cmd.zrange = GetParam("zrange");
	if (strcmp(cmd.zrange,"autoscale"))
		gd.z_autoscale=FALSE;
	if (!(sscanf(GetParam("zrange"), "%lf:%lf", &gd.zmin, &gd.zmax) == 2))
		zrangeflag=0;
	if ( !gd.z_autoscale && zrangeflag==0)
		error("\nGet_Parameters: zrange must be autoscale or in the form zmin:zmax\n\n");
	if (zrangeflag != 0 && gd.zmin >= gd.zmax)
		error("\nGet_Parameters: zmin must be < zmax : (zmin,zmax)=(%lf,%lf)\n\n",gd.zmin,gd.zmax);

	cmd.usingcolumns = GetParam("usingcolumns");
	if (scanopt(cmd.data_analysis_type,"thermo-avg") && scanopt(cmd.options, "errorbar")) {
		if (!(sscanf(cmd.usingcolumns, "%d:%d:%d", 
			&gd.column1, &gd.column2, &gd.column3) == 3))
			error("\nstartrun_cmdline: usingcolumns must be in the form c1:c2:c3\n\n");
	} else
		if (!(sscanf(cmd.usingcolumns, "%d:%d", &gd.column1, &gd.column2) == 2))
			error("\nstartrun_cmdline: usingcolumns must be in the form c1:c2\n\n");

	cmd.usingrows = GetParam("usingrows");
	if (!(sscanf(cmd.usingrows, "%d:%d", &gd.row1, &gd.row2) == 2))
		error("\nGet_Parameters: usingrows must be in the form r1:r2\n\n");

	cmd.xlabel = GetParam("xlabel");
	cmd.ylabel = GetParam("ylabel");
	cmd.plotlabel = GetParam("plotlabel");
	cmd.labelfontsize = GetdParam("labelfontsize");
	cmd.withdots = GetbParam("withdots");
	cmd.withsymbols = GetbParam("withsymbols");
	cmd.symboltype	= GetiParam("symboltype");
	cmd.symbolsize = GetdParam("symbolsize");

	cmd.labelfontweight = GetiParam("labelfontweight");
	cmd.symbolweight	= GetiParam("symbolweight");
	cmd.nlsize			= GetdParam("nlsize");
	cmd.linewidth		= GetiParam("linewidth");
	cmd.axeswidth		= GetiParam("axeswidth");
	cmd.symbolcolor		= GetiParam("symbolcolor");

	cmd.pl_a = GetParam("a");
	cmd.pl_dev = GetParam("dev");
	cmd.pl_geo = GetParam("geo");
	cmd.pl_ori = GetParam("ori");
	cmd.pl_bg = GetParam("bg");
	cmd.pl_ncol0 = GetParam("ncol0");
	cmd.pl_ncol1 = GetParam("ncol1");

	PrintParameterFile(parameter_null);
}

#undef parameter_null
#undef logfile

local void Cmd_to_Global(void)
{
// Seran activados cuando se necesite calcular la fuerza entre particulas.
//	gdtreegrav.usequad = cmd.usequad;
//	gdtreegrav.theta = cmd.theta;
//	gdtreegrav.eps = cmd.eps;
}

// Seran activados cuando se necesite calcular la fuerza entre particulas.
local void ComputeParameters(void)
{
//    gdtreegrav.eps2 = gdtreegrav.eps * gdtreegrav.eps;
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

// Seran activados cuando se necesite calcular la fuerza entre particulas.
//	BPName(cmd.usequad,"usequad");
//	RPName(cmd.theta,"theta");
//	RPName(cmd.eps,"eps");

	SPName(cmd.in,"in",100);
	SPName(cmd.infmt,"infmt",100);
	SPName(cmd.out,"out",100);
	SPName(cmd.outfmt,"outfmt",100);
	SPName(cmd.basedir,"basedir",100);

	IPName(cmd.isnap,"isnap");
	IPName(cmd.fsnap,"fsnap");
//	RPName(dtime,"dtime");
//	RPName(cmd.density,"density");
//	IPName(cmd.stepEquil,"stepEquil");

/*
	IPName(cmd.stepAvgRhoAxes,"stepAvgRhoAxes");
	IPName(cmd.sizeHistRhoAxes,"sizeHistRhoAxes");
	IPName(cmd.stepAvgRhoAxes,"stepAvgNFrecAxes");
	IPName(cmd.sizeHistRhoAxes,"sizeHistNFrecAxes");

	IPName(cmd.stepAvgVel,"stepAvgVel");
	IPName(cmd.sizeHistVel,"sizeHistVel");
	RPName(cmd.rangeVel,"rangeVel");

	IPName(cmd.stepAvgRdf,"stepAvgRdf");
	IPName(cmd.sizeHistRdf,"sizeHistRdf");
	RPName(cmd.rangeRdf,"rangeRdf");
*/

	IPName(cmd.stepAvgRhoTheta,"stepAvgRhoTheta");
	IPName(cmd.sizeHistRhoTheta,"sizeHistRhoTheta");
	RPName(cmd.RhoDeltaZ,"RhoDeltaZ");
	RPName(cmd.RhoR,"RhoR");
	RPName(cmd.RhoDeltaR,"RhoDeltaR");
	RPName(cmd.RMax,"RMax");

	RPName(cmd.xmin,"xmin");
	RPName(cmd.xmax,"xmax");
	RPName(cmd.ymin,"ymin");
	RPName(cmd.ymax,"ymax");
	RPName(cmd.zmin,"zmin");
	RPName(cmd.zmax,"zmax");

	RPName(cmd.ThetaMin,"ThetaMin");
	RPName(cmd.ThetaMax,"ThetaMax");

	IPName(cmd.stepAvgRho,"stepAvgRho");
	IPName(cmd.sizeHistRho,"sizeHistRho");

	IPName(cmd.stepAvgVcR,"stepAvgVcR");
	IPName(cmd.sizeHistVcR,"sizeHistVcR");
	RPName(cmd.rangeR,"rangeR");

//	BPName(cmd.computeTransport,"computeTransport");
//	IPName(cmd.stepAcf,"stepAcf");
	IPName(cmd.limitAcfAv,"limitAcfAv");
	IPName(cmd.nBuffAcf,"nBuffAcf");
	IPName(cmd.nValAcf,"nValAcf");

	SPName(cmd.options,"options",100);
//	RPName(dtout,"dtout");
//	IPName(nbody,"nbody");
	SPName(cmd.data_analysis_type,"analysis_type",100);

	SPName(cmd.bodiesID,"bodiesID",256);
	SPName(cmd.bodiesSets,"bodiesSets",256);

	RPName(cmd.reductionFac,"reductionFac");
	SPName(cmd.foffile,"foffile",100);
	RPName(cmd.gravityConstant,"gravityConstant");

	SPName(cmd.xrange,"xrange",100);
	SPName(cmd.yrange,"yrange",100);
	SPName(cmd.zrange,"zrange",100);

	SPName(cmd.usingcolumns,"usingcolumns",100);
	SPName(cmd.usingrows,"usingrows",100);
	SPName(cmd.xlabel,"xlabel",100);
	SPName(cmd.ylabel,"ylabel",100);
	RPName(cmd.labelfontsize,"labelfontsize");
	SPName(cmd.plotlabel,"plotlabel",100);
	BPName(cmd.withsymbols,"withsymbols");
	BPName(cmd.withdots,"withdots");
	IPName(cmd.symboltype,"symboltype");
	RPName(cmd.symbolsize,"symbolsize");

	IPName(cmd.labelfontweight,"labelfontweight");
	RPName(cmd.nlsize,"nlsize");
	IPName(cmd.linewidth,"linewidth");
	IPName(cmd.axeswidth,"axeswidth");
	IPName(cmd.symbolweight,"symbolweight");
	IPName(cmd.symbolcolor,"symbolcolor");

	SPName(cmd.pl_a,"a",1);
	SPName(cmd.pl_dev,"dev",1);
	SPName(cmd.pl_geo,"geo",1);
	SPName(cmd.pl_ori,"ori",1);
	SPName(cmd.pl_bg,"bg",1);
	SPName(cmd.pl_ncol0,"ncol0",1);
	SPName(cmd.pl_ncol1,"ncol1",1);

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
//                *buf2=NULL;							// Original statement
//                *buf2=(int) NULL;
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
        }					// end of while
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1); 
    }						// end of the if((fd=fopen...))
  
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

    sprintf(buf,"%s%s",fname,"-analysis_galaxy");
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

// Seran activados cuando se necesite calcular la fuerza entre particulas.
//        fprintf(fdout,FMTT,"usequad",cmd.usequad ? "true" : "false");
//        fprintf(fdout,FMTR,"theta",cmd.theta);
//        fprintf(fdout,FMTR,"eps",cmd.eps);

        fprintf(fdout,FMTT,"in",cmd.in);
        fprintf(fdout,FMTT,"infmt",cmd.infmt);
        fprintf(fdout,FMTT,"out",cmd.out);
        fprintf(fdout,FMTT,"outfmt",cmd.outfmt);
        fprintf(fdout,FMTT,"basedir",cmd.basedir);

        fprintf(fdout,FMTI,"isnap",cmd.isnap);
        fprintf(fdout,FMTI,"fsnap",cmd.fsnap);

//        fprintf(fdout,FMTI,"stepEquil",stepEquil);

//        fprintf(fdout,FMTR,"dtime",dtime);

//        fprintf(fdout,FMTR,"dtout",dtout);

        fprintf(fdout,FMTI,"stepAvgRhoTheta",cmd.stepAvgRhoTheta);
        fprintf(fdout,FMTI,"sizeHistRhoTheta",cmd.sizeHistRhoTheta);
        fprintf(fdout,FMTR,"RhoDeltaZ",cmd.RhoDeltaZ);
        fprintf(fdout,FMTR,"RhoR",cmd.RhoR);
        fprintf(fdout,FMTR,"RhoDeltaR",cmd.RhoDeltaR);
        fprintf(fdout,FMTR,"RMax",cmd.RMax);

        fprintf(fdout,FMTR,"xmin",cmd.xmin);
        fprintf(fdout,FMTR,"xmax",cmd.xmax);
        fprintf(fdout,FMTR,"ymin",cmd.ymin);
        fprintf(fdout,FMTR,"ymax",cmd.ymax);
        fprintf(fdout,FMTR,"zmin",cmd.zmin);
        fprintf(fdout,FMTR,"zmax",cmd.zmax);

        fprintf(fdout,FMTR,"ThetaMin",cmd.ThetaMin);
        fprintf(fdout,FMTR,"ThetaMax",cmd.ThetaMax);

        fprintf(fdout,FMTI,"stepAvgRho",cmd.stepAvgRho);
        fprintf(fdout,FMTI,"sizeHistRho",cmd.sizeHistRho);

        fprintf(fdout,FMTI,"stepAvgVcR",cmd.stepAvgVcR);
        fprintf(fdout,FMTI,"sizeHistVcR",cmd.sizeHistVcR);
        fprintf(fdout,FMTR,"rangeR",cmd.rangeR);

//        fprintf(fdout,FMTT,"computeTransport",
//				cmd.computeTransport ? "true" : "false");
//        fprintf(fdout,FMTI,"stepAcf",cmd.stepAcf);
        fprintf(fdout,FMTI,"limitAcfAv",cmd.limitAcfAv);
        fprintf(fdout,FMTI,"nBuffAcf",cmd.nBuffAcf);
        fprintf(fdout,FMTI,"nValAcf",cmd.nValAcf);

        fprintf(fdout,FMTT,"options",cmd.options);
//        fprintf(fdout,FMTI,"nbody",nbody);
        fprintf(fdout,FMTT,"analysis_type",cmd.data_analysis_type);

        fprintf(fdout,FMTT,"bodiesID",cmd.bodiesID);
        fprintf(fdout,FMTT,"bodiesSets",cmd.bodiesSets);

        fprintf(fdout,FMTR,"reductionFac",cmd.reductionFac);
        fprintf(fdout,FMTT,"foffile",cmd.foffile);
        fprintf(fdout,FMTR,"gravityConstant",cmd.gravityConstant);

        fprintf(fdout,FMTT,"xrange",cmd.xrange);
        fprintf(fdout,FMTT,"yrange",cmd.yrange);
        fprintf(fdout,FMTT,"zrange",cmd.zrange);

        fprintf(fdout,FMTT,"usingcolumns",cmd.usingcolumns);
        fprintf(fdout,FMTT,"usingrows",cmd.usingrows);

        fprintf(fdout,FMTT,"xlabel",cmd.xlabel);
        fprintf(fdout,FMTT,"ylabel",cmd.ylabel);
        fprintf(fdout,FMTT,"plotlabel",cmd.plotlabel);
        fprintf(fdout,FMTR,"labelfontsize",cmd.labelfontsize);
        fprintf(fdout,FMTT,"withdots",cmd.withdots ? "true" : "false");
        fprintf(fdout,FMTT,"withsymbols",cmd.withsymbols ? "true" : "false");
        fprintf(fdout,FMTI,"symboltype",cmd.symboltype);
        fprintf(fdout,FMTR,"symbolsize",cmd.symbolsize);

        fprintf(fdout,FMTI,"labelfontweight",cmd.labelfontweight);
        fprintf(fdout,FMTR,"nlsize",cmd.nlsize);
        fprintf(fdout,FMTI,"linewidth",cmd.linewidth);
        fprintf(fdout,FMTI,"axeswidth",cmd.axeswidth);
        fprintf(fdout,FMTI,"symbolweight",cmd.symbolweight);
        fprintf(fdout,FMTI,"symbolcolor",cmd.symbolcolor);

        fprintf(fdout,FMTT,"a",cmd.pl_a);
        fprintf(fdout,FMTT,"dev",cmd.pl_dev);
        fprintf(fdout,FMTT,"geo",cmd.pl_geo);
        fprintf(fdout,FMTT,"ori",cmd.pl_ori);
        fprintf(fdout,FMTT,"bg",cmd.pl_bg);
        fprintf(fdout,FMTT,"ncol0",cmd.pl_ncol0);
        fprintf(fdout,FMTT,"ncol1",cmd.pl_ncol1);

        fprintf(fdout,"\n\n");
							// is important "\n\n" so can be read as
    }						// input parameter file
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR


local void CheckInOutFmt(void)
{
	gd.in_long_fmt = 0;
//	gd.out_long_fmt = 0;
	
	if (strcmp(cmd.infmt,"snap-ascii-long") == 0)
		gd.in_long_fmt = 1;
//	if (strcmp(cmd.outfmt,"snap-ascii-long") == 0)
//		gd.out_long_fmt = 1;
	if (strcmp(cmd.infmt,"gadget11-ascii-long") == 0)
		gd.in_long_fmt = 1;
//	if (strcmp(cmd.outfmt,"gadget11-normal-body-ascii-long") == 0)
//		gd.out_long_fmt = 1;
	if (strcmp(cmd.infmt,"heitmann-ascii-long") == 0)
		gd.in_long_fmt = 1;
//	if (strcmp(cmd.outfmt,"powmes-ascii-long") == 0)
//		gd.out_long_fmt = 1;
	
	//	if (strcmp(cmd.outfmt,"powmes-ascii-long") == 0 && 
	//		strcmp(cmd.infmt,"gadget11-ascii-long") != 0)
	//		error("\nCheckParameters: you should suply an input long-file\n\n");
	
	//	if (strcmp(cmd.outfmt,"snap-ascii-long") == 0 && 
	//		strcmp(cmd.infmt,"gadget11-ascii-long") != 0)
	//		error("\nCheckParameters: you should suply an input long-file\n\n");
	
//	if (gd.in_long_fmt != gd.out_long_fmt)
//		error("\nCheckInOutFmt: you have to give infmt and ofmt of the same type\n\n");
}
