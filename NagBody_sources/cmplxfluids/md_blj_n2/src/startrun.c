/*==============================================================================
	MODULE: startrun.c				[md_blj_n2]
	Written by: M.A. Rodriguez-Meza
	Starting date: Enero 2005
	Purpose: Initialize md_blj_n2
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

	Major revisions: December 2006; November 2008;
	Copyright: (c) 2005-2012 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#include "globaldefs.h"

local void Header_to_Global(void);
local void Compute_nbody(void);
local void Compute_Box_Units(void);
local void Compute_Parameters(void);
local void Print_Parameters(void);
//local void testdata_sp_01(void);
//local void testdata_sp_02(void);
//local void testdata_sp_03(void);
//local void testdata_sp_04(void);
//local void testdata_sp_05(void);
//local void testdata_sp_06(void);
//local void testdata(void);
//local void MinimumDistanceBBodies(real);
local void SetRcut(void);
local void SetBodyChemPot(bodyptr *);
local void ReadParameterFile(char *);
local void PrintParameterFile(char *);
local void forcecalc_method_string_to_int(string,int *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_restorefile(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

local void InputForcePotTable(void);

														// CHECK 2D --- OK!!!
void StartRun(string head0, string head1, string head2, string head3)
{
    double cpustart;
	bodyptr p;

    cpustart = cputime();                       

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

//    LicDriver(gd.headline0);

	gd.stopflag = 0;
	gd.cputotforce = 0;
	gd.cputotout = 0.;
	gd.cputotal = 0.;

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

	StartOutput();
	fprintf(stdout,"\n\nStartRun CPU time: %g\n",cputime()-cpustart);
}

local void startrun_parameterfile(void)					// CHECK 2D --- OK!!!
{
	ReadParameterFile(cmd.paramfile);
	if (strnull(cmd.restorefile))
		startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-md_blj_n2"

local void startrun_cmdline(void)						// CHECK 2D --- OK!!!
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)					// CHECK 2D --- OK!!!
{
    cmd.icfile = GetParam("icfile");
    cmd.icfilefmt = GetParam("icfilefmt");
    cmd.snapoutfile = GetParam("snapout");
    cmd.snapoutfilefmt = GetParam("snapoutfmt");
    cmd.statefile = GetParam("statefile");
	cmd.stepState = GetiParam("stepState");
    cmd.restorefile = GetParam("restorefile");

	cmd.temperature = GetdParam("temperature");
	cmd.adjustTemperature = GetbParam("adjustTemperature");
	cmd.stepAdjustTemperature = GetiParam("stepAdjustTemperature");
	cmd.adjustCenterOfMass = GetbParam("adjustCenterOfMass");
	cmd.density = GetdParam("density");
	cmd.stepEquil = GetiParam("stepEquil");
	cmd.stepAvg = GetiParam("stepAvg");

// borrar cuando ya este seguro de que no sirven...
//	cmd.printSnap = GetbParam("printSnap");
//	cmd.stepSnap = GetiParam("stepSnap");

	cmd.computeRhoAxes = GetbParam("computeRhoAxes");
	cmd.stepAvgRhoAxes = GetiParam("stepAvgRhoAxes");
	cmd.stepRhoAxes = GetiParam("stepRhoAxes");
	cmd.sizeHistRhoAxes = GetiParam("sizeHistRhoAxes");

	cmd.computeNFrecAxes = GetbParam("computeNFrecAxes");
	cmd.stepAvgNFrecAxes = GetiParam("stepAvgNFrecAxes");
	cmd.stepNFrecAxes = GetiParam("stepNFrecAxes");
	cmd.sizeHistNFrecAxes = GetiParam("sizeHistNFrecAxes");

	cmd.computeVelDist = GetbParam("computeVelDist");
	cmd.stepAvgVel = GetiParam("stepAvgVel");
	cmd.stepVel = GetiParam("stepVel");
	cmd.sizeHistVel = GetiParam("sizeHistVel");
	cmd.rangeVel = GetdParam("rangeVel");

	cmd.computeRdf = GetbParam("computeRdf");
	cmd.stepAvgRdf = GetiParam("stepAvgRdf");
	cmd.stepRdf = GetiParam("stepRdf");
	cmd.sizeHistRdf = GetiParam("sizeHistRdf");
	cmd.rangeRdf = GetdParam("rangeRdf");

	cmd.computePressAxes = GetbParam("computePressAxes");
	cmd.stepAvgPressAxes = GetiParam("stepAvgPressAxes");
	cmd.stepPressAxes = GetiParam("stepPressAxes");
	cmd.sizeHistPressAxes = GetiParam("sizeHistPressAxes");

	cmd.computeChemPot = GetbParam("computeChemPot");
	cmd.stepAvgChemPot = GetiParam("stepAvgChemPot");
	cmd.stepChemPot = GetiParam("stepChemPot");
	cmd.sizeHistChemPot = GetiParam("sizeHistChemPot");
	cmd.numTestBodies = GetiParam("numTestBodies");

	cmd.computeDiffusion = GetbParam("computeDiffusion");
	cmd.stepDiffuse = GetiParam("stepDiffuse");
	cmd.stepAvgDiffuse = GetiParam("stepAvgDiffuse");
	cmd.nBuffDiffuse = GetiParam("nBuffDiffuse");
	cmd.nValDiffuse = GetiParam("nValDiffuse");

	cmd.computeVelAcf = GetbParam("computeVelAcf");
	cmd.computeBulkViscosity = GetbParam("computeBulkViscosity");

	cmd.computeTransport = GetbParam("computeTransport");
	cmd.stepAcf = GetiParam("stepAcf");
	cmd.stepAvgAcf = GetiParam("stepAvgAcf");
	cmd.nBuffAcf = GetiParam("nBuffAcf");
	cmd.nValAcf = GetiParam("nValAcf");

	cmd.computeSTCorr = GetbParam("computeSTCorr");
	cmd.stepCorr = GetiParam("stepCorr");
	cmd.stepAvgCorr = GetiParam("stepAvgCorr");
	cmd.nBuffCorr = GetiParam("nBuffCorr");
	cmd.nValCorr = GetiParam("nValCorr");
	cmd.nFunCorr = GetiParam("nFunCorr");

	cmd.lattCorr_kx = GetdParam("lattCorr_kx");
	cmd.lattCorr_ky = GetdParam("lattCorr_ky");
#ifdef THREEDIM
	cmd.lattCorr_kz = GetdParam("lattCorr_kz");
#endif

	cmd.dtimestr = GetParam("dtime");
	cmd.icModel = GetiParam("icModel");
	cmd.intMethod = GetiParam("intMethod");

	cmd.tstop = GetdParam("tstop");

	cmd.dtoutstr = GetParam("dtout");

	cmd.dtoutinfostr = GetParam("dtoutinfo");

	cmd.forcecalc_method = GetParam("forcecalc_method");
	cmd.potType = GetiParam("potentialType");
    cmd.fnamePot = GetParam("fnamePot");
	cmd.options = GetParam("options");

	cmd.eps11 = GetdParam("eps11");
	cmd.eps12 = GetdParam("eps12");
	cmd.eps22 = GetdParam("eps22");
	cmd.sigma11 = GetdParam("sigma11");
	cmd.sigma12 = GetdParam("sigma12");
	cmd.sigma22 = GetdParam("sigma22");
	cmd.Rcut11 = GetdParam("Rcut11");
	cmd.Rcut12 = GetdParam("Rcut12");
	cmd.Rcut22 = GetdParam("Rcut22");

	cmd.seed=GetiParam("seed");

	cmd.unitCells = GetParam("unitCells");

	cmd.nbodyprop = GetParam("nbodyprop");
	cmd.massprop = GetParam("massprop");

#ifdef THREEDIM
	cmd.LxLyprop = GetParam("LxLyprop");
#else
	cmd.Lx = GetParam("Lx");
#endif

}

#undef parameter_null

#define logfile			"md_blj.log"

local void startrun_Common(void)						// CHECK 2D --- OK!!!
{
	real dt1, dt2;
	bool exist_snap;
	int ndim, step=0;

	if (strnull(cmd.restorefile))
		strcpy(gd.mode,"w");
	else
		strcpy(gd.mode,"a");

	if(!(gd.outlog=fopen(logfile,gd.mode)))
		error("\nstart_Common: error opening file '%s' \n",logfile);

//	if (strnull(cmd.icfilefmt))
		gd.headerfmt = "snap-blj-ascii";
//	else
//		gd.infilefmt = cmd.icfilefmt;

	if (strnull(cmd.restorefile)) {
		gd.dtime = (sscanf(cmd.dtimestr, "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtimestr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtime : dt2 must be finite\n");

		gd.dtout = (sscanf(cmd.dtoutstr, "%lf/%lf", &dt1, &dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtoutstr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtout : dt2 must be finite\n");

		gd.dtoutinfo=(sscanf(cmd.dtoutinfostr,"%lf/%lf",&dt1,&dt2) == 2 ?
				dt1/dt2 : atof(cmd.dtoutinfostr));
		if ( dt2 == 0. )
			error("\n\nstartrun_Common: dtoutinfo :dt2 must be finite\n");

		if (! strnull(cmd.icfile)) {
			inputdata(cmd.icfile,cmd.icfilefmt,gd.headerfmt,step,&gd.nbody,
				&ndim,&gd.tnow,&exist_snap,&hdr,cmd.options,gd.model_comment);
			Header_to_Global();
			CheckParameters();
			if (cmd.potType==POTTYPE_FILE) InputForcePotTable();
			Compute_Box_Units();
			Compute_Parameters();
		} else {

#ifdef THREEDIM
			if (!(sscanf(cmd.unitCells,"%ld:%ld:%ld", 
				&gd.numUcell[0], &gd.numUcell[1], &gd.numUcell[2]) == 3))
				error("\nstartrun_cmdline: %s",
					"unitCells must be in the form of ncx:ncy:ncz\n\n");
#else
			if (!(sscanf(cmd.unitCells,"%ld:%ld", 
				&gd.numUcell[0], &gd.numUcell[1]) == 2))
				error("\nstartrun_cmdline: %s",
					"unitCells must be in the form of ncx:ncy\n\n");
#endif
			if (!(sscanf(cmd.nbodyprop,"%ld/%ld", &gd.nbody1, &gd.nbody2) == 2))
				error("\nstartrun_parameterfile: %s",
					"nbodyprop must be in the form of nbody1/nbody2\n\n");

			if (cmd.icModel==3 || cmd.icModel==6) {
//			if (cmd.icModel==3) {
				gd.nbc1 = gd.nbody1;
				gd.nbc2 = gd.nbody2;
			}

			if (!(sscanf(cmd.massprop, "%lf/%lf", &gd.mass1, &gd.mass2) == 2))
				error("\nstartrun_parameterfile: %s",
					"massprop must be in the form of mass1/mass2\n\n");

#ifdef THREEDIM
			if (!(sscanf(cmd.LxLyprop, "%lf/%lf", &gd.Lx, &gd.Ly) == 2))
				error("\nstartrun_parameterfile: %s",
						"LxLyprop must be in the form of Lx/Ly\n\n");
#else
			if (!(sscanf(cmd.Lx, "%lf", &gd.Lx) == 1))
				error("\nstartrun_parameterfile: %s",
						"Lx must be in the form of Lx\n\n");
#endif
			idum=cmd.seed;
			xsrandom(idum);
			CheckParameters();
			if (cmd.potType==POTTYPE_FILE) InputForcePotTable();
			Compute_nbody();
			Compute_Box_Units();
			testdata();
			Compute_Parameters();
			gd.tnow = 0.0;
		}
		gdtree.rsize = 1.0;                            
		gd.nstep = 0;

		gd.nstepOld = 0;					// In order to restore work ...
		gd.nstepNew = 0;

		gd.tout = gd.tnow;
		gd.toutinfo = gd.tnow;
	} else
		startrun_restorefile();

	SetRcut();
	if (cmd.computeChemPot) SetBodyChemPot(&bodytab);

	forcecalc_method_string_to_int(cmd.forcecalc_method, &gd.forcemethod_int);
	Print_Parameters();
}

#undef logfile

local void startrun_restorefile(void)					// CHECK 2D --- OK!!!
{
	fprintf(gd.outlog,"\n\nAdded after restart from restart file\n");
	restorestate(cmd.restorefile);
	fprintf(gd.outlog,"\nnbody=%d\n",gd.nbody);
	fprintf(gd.outlog,"\nnbody1=%d\n",gd.nbody1);
	fprintf(gd.outlog,"\nnbody2=%d\n",gd.nbody2);
	fprintf(gd.outlog,"\nnbc1, nbc2 = %d, %d\n",gd.nbc1, gd.nbc2);
	fprintf(gd.outlog,"\ntemperature & density : %g %g\n",
		cmd.temperature,cmd.density);

	startrun_ParamStat();
	CheckParameters();

	if (scanopt(cmd.options, "new-tout"))
		gd.tout = gd.tnow + gd.dtout;
	if (scanopt(cmd.options, "new-toutinfo"))
		gd.toutinfo = gd.tnow + gd.dtoutinfo;

	if (cmd.potType==POTTYPE_FILE) InputForcePotTable();
	Compute_Box_Units();
	Compute_Parameters();
}

local void startrun_ParamStat(void)						// CHECK 2D --- OK!!!
{
	real dt1, dt2;

	if (GetParamStat("forcecalc_method") & ARGPARAM) {
		cmd.forcecalc_method = GetParam("forcecalc_method");
		fprintf(gd.outlog,"\n\nrunning now %s force calc. method ...\n",
				cmd.forcecalc_method);
	}

	if (GetParamStat("potentialType") & ARGPARAM) { 
		cmd.potType = GetiParam("potentialType");
		fprintf(gd.outlog,"\n\nrunning now force model type: %d ...\n",
				cmd.potType);
		if (cmd.potType==POTTYPE_FILE)
			if (GetParamStat("fnamePot") & ARGPARAM)
				cmd.fnamePot = GetParam("fnamePot");
	}

	if (GetParamStat("snapout") & ARGPARAM) 
        cmd.snapoutfile = GetParam("snapout");
	if (GetParamStat("snapoutfmt") & ARGPARAM) 
        cmd.snapoutfilefmt = GetParam("snapoutfmt");

	if (GetParamStat("statefile") & ARGPARAM) 
        cmd.statefile = GetParam("statefile");
	if (GetParamStat("stepState") & ARGPARAM)
		cmd.stepState = GetiParam("stepState");    
    
    
	if (GetParamStat("temperature") & ARGPARAM) 
		cmd.temperature = GetdParam("temperature");
	if (GetParamStat("adjustTemperature") & ARGPARAM) 
		cmd.adjustTemperature = GetbParam("adjustTemperature");
	if (GetParamStat("stepAdjustTemperature") & ARGPARAM)
		cmd.stepAdjustTemperature = GetiParam("stepAdjustTemperature");

	if (GetParamStat("adjustCenterOfMass") & ARGPARAM) 
		cmd.adjustCenterOfMass = GetbParam("adjustCenterOfMass");
	if (GetParamStat("density") & ARGPARAM) 
		cmd.density = GetdParam("density");

	if (GetParamStat("Rcut11") & ARGPARAM) 
		cmd.Rcut11 = GetdParam("Rcut11");
	if (GetParamStat("Rcut12") & ARGPARAM) 
		cmd.Rcut12 = GetdParam("Rcut12");
	if (GetParamStat("Rcut22") & ARGPARAM) 
		cmd.Rcut22 = GetdParam("Rcut22");

	if (GetParamStat("stepEquil") & ARGPARAM)
		cmd.stepEquil = GetiParam("stepEquil");

	gd.stepAvgFlag = FALSE;						// In order to restore work ...
	gd.stepAvgStatus = FALSE;						// In order to restore work ...
	if (GetParamStat("stepAvg") & ARGPARAM) {
		gd.stepAvgOld = cmd.stepAvg;			// In order to restore work ...
		cmd.stepAvg = GetiParam("stepAvg");
		gd.stepAvgFlag = TRUE;					// In order to restore work ...
		gd.stepAvgStatus = TRUE;					// In order to restore work ...
	}

// borrar cuando ya este seguro de que no sirven...
//	if (GetParamStat("printSnap") & ARGPARAM)
//		cmd.printSnap = GetbParam("printSnap");
//	if (GetParamStat("stepSnap") & ARGPARAM)
//		cmd.stepSnap = GetiParam("stepSnap");

	if (GetParamStat("stepAvgRhoAxes") & ARGPARAM)
		cmd.stepAvgRhoAxes = GetiParam("stepAvgRhoAxes");         
	if (GetParamStat("stepRhoAxes") & ARGPARAM)
		cmd.stepRhoAxes = GetiParam("stepRhoAxes");         
	if (GetParamStat("sizeHistRhoAxes") & ARGPARAM)
		cmd.sizeHistRhoAxes = GetiParam("sizeHistRhoAxes");

	if (GetParamStat("stepAvgNFrecAxes") & ARGPARAM)
		cmd.stepAvgNFrecAxes = GetiParam("stepAvgNFrecAxes");         
	if (GetParamStat("stepNFrecAxes") & ARGPARAM)
		cmd.stepNFrecAxes = GetiParam("stepNFrecAxes");         
	if (GetParamStat("sizeHistNFrecAxes") & ARGPARAM)
		cmd.sizeHistNFrecAxes = GetiParam("sizeHistNFrecAxes");

	if (GetParamStat("computeVelDist") & ARGPARAM)
		cmd.computeVelDist = GetbParam("computeVelDist");
	if (GetParamStat("stepAvgVel") & ARGPARAM)
		cmd.stepAvgVel = GetiParam("stepAvgVel");         
	if (GetParamStat("stepVel") & ARGPARAM)
		cmd.stepVel = GetiParam("stepVel");         
	if (GetParamStat("sizeHistVel") & ARGPARAM)
		cmd.sizeHistVel = GetiParam("sizeHistVel");
	if (GetParamStat("rangeVel") & ARGPARAM)
		cmd.rangeVel = GetdParam("rangeVel");

	if (GetParamStat("computeRdf") & ARGPARAM)
		cmd.computeRdf = GetbParam("computeRdf");
	if (GetParamStat("stepAvgRdf") & ARGPARAM)
		cmd.stepAvgRdf = GetiParam("stepAvgRdf");         
	if (GetParamStat("stepRdf") & ARGPARAM)
		cmd.stepRdf = GetiParam("stepRdf");         
	if (GetParamStat("sizeHistRdf") & ARGPARAM)
		cmd.sizeHistRdf = GetiParam("sizeHistRdf");
	if (GetParamStat("rangeRdf") & ARGPARAM)
		cmd.rangeRdf = GetdParam("rangeRdf");

	if (GetParamStat("stepAvgChemPot") & ARGPARAM)
		cmd.stepAvgChemPot = GetiParam("stepAvgChemPot");
	if (GetParamStat("stepChemPot") & ARGPARAM)
		cmd.stepChemPot = GetiParam("stepChemPot");
	if (GetParamStat("sizeHistChemPot") & ARGPARAM)
		cmd.sizeHistChemPot = GetiParam("sizeHistChemPot");
	if (GetParamStat("computeChemPot") & ARGPARAM)
		cmd.computeChemPot = GetbParam("computeChemPot");
	if (GetParamStat("numTestBodies") & ARGPARAM)
		cmd.numTestBodies = GetiParam("numTestBodies");

	if (GetParamStat("computeDiffusion") & ARGPARAM)
		cmd.computeDiffusion = GetbParam("computeDiffusion");
	if (GetParamStat("stepDiffuse") & ARGPARAM)
		cmd.stepDiffuse = GetiParam("stepDiffuse");
	if (GetParamStat("stepAvgDiffuse") & ARGPARAM)
		cmd.stepAvgDiffuse = GetiParam("stepAvgDiffuse");
	if (GetParamStat("nBuffDiffuse") & ARGPARAM)
		cmd.nBuffDiffuse = GetiParam("nBuffDiffuse");
	if (GetParamStat("nValDiffuse") & ARGPARAM)
		cmd.nValDiffuse = GetiParam("nValDiffuse");

	if (GetParamStat("computeVelAcf") & ARGPARAM)
		cmd.computeVelAcf = GetbParam("computeVelAcf");

	if (GetParamStat("computeBulkViscosity") & ARGPARAM)
		cmd.computeBulkViscosity = GetbParam("computeBulkViscosity");

	if (GetParamStat("computeTransport") & ARGPARAM)
		cmd.computeTransport = GetbParam("computeTransport");
	if (GetParamStat("stepAcf") & ARGPARAM)
		cmd.stepAcf = GetiParam("stepAcf");
	if (GetParamStat("stepAvgAcf") & ARGPARAM)
		cmd.stepAvgAcf = GetiParam("stepAvgAcf");
	if (GetParamStat("nBuffAcf") & ARGPARAM)
		cmd.nBuffAcf = GetiParam("nBuffAcf");
	if (GetParamStat("nValAcf") & ARGPARAM)
		cmd.nValAcf = GetiParam("nValAcf");

	if (GetParamStat("computeSTCorr") & ARGPARAM)
		cmd.computeSTCorr = GetbParam("computeSTCorr");
	if (GetParamStat("stepCorr") & ARGPARAM)
		cmd.stepCorr = GetiParam("stepCorr");
	if (GetParamStat("stepAvgCorr") & ARGPARAM)
		cmd.stepAvgCorr = GetiParam("stepAvgCorr");
	if (GetParamStat("nBuffCorr") & ARGPARAM)
		cmd.nBuffCorr = GetiParam("nBuffCorr");
	if (GetParamStat("nValCorr") & ARGPARAM)
		cmd.nValCorr = GetiParam("nValCorr");
	if (GetParamStat("nFunCorr") & ARGPARAM)
		cmd.nFunCorr = GetiParam("nFunCorr");

	if (GetParamStat("lattCorr_kx") & ARGPARAM)
		cmd.lattCorr_kx = GetdParam("lattCorr_kx");
	if (GetParamStat("lattCorr_ky") & ARGPARAM)
		cmd.lattCorr_ky = GetdParam("lattCorr_ky");
#ifdef THREEDIM
	if (GetParamStat("lattCorr_kz") & ARGPARAM)
		cmd.lattCorr_kz = GetdParam("lattCorr_kz");
#endif

	if (GetParamStat("options") & ARGPARAM)
		cmd.options = GetParam("options");
	if (GetParamStat("tstop") & ARGPARAM)
		cmd.tstop = GetdParam("tstop");

	if (GetParamStat("dtout") & ARGPARAM) {
		cmd.dtoutstr = GetParam("dtout");
		gd.dtout = (sscanf(cmd.dtoutstr, "%lf/%lf", &dt1, &dt2) == 2 ?
			dt1/dt2 : atof(cmd.dtoutstr));
		if ( dt2 == 0. )
			error("\n\nstartrun_ParamStat: dtout : dt2 must be finite\n");
	}
	if (GetParamStat("dtoutinfo") & ARGPARAM) {
		cmd.dtoutinfostr = GetParam("dtoutinfo");
		gd.dtoutinfo = (sscanf(cmd.dtoutinfostr,"%lf/%lf",&dt1,&dt2) == 2 ?
			dt1/dt2 : atof(cmd.dtoutinfostr));
		if ( dt2 == 0. )
			error("\n\nstartrun_ParamStat: dtoutinfo : dt2 must be finite\n");
	}
}

local void CheckParameters(void)						// CHECK 2D --- OK!!!
{
	if (!strnull(cmd.restorefile) && !strnull(cmd.icfile))
//		error("\nstartrun_parameterfile: %s\n\n",
//			"Use only one of the options restorefile or icfile");
		fprintf(stdout,"\nCheckParameters: Warning! : %s\n\n",
			"You are using options restorefile and icfile at the same time");

	if (cmd.numTestBodies==0)
		error("\n\nstartrun_cmdline: numTestBodies must be != 0\n");

	if ( cmd.stepAdjustTemperature <= 0 )
		error("\n\nstartrun_parameterfile: stepAdjustTemperature must be > 0\n");
	if ( cmd.stepAvg <= 0 )
		error("\n\nstartrun_parameterfile: stepAvg must be > 0\n");

// borrar cuando ya este seguro de que no sirven...
//	if ( cmd.stepSnap <= 0 )
//		error("\n\nstartrun_parameterfile: stepSnap must be > 0\n");

	if ( cmd.stepState <= 0 )
		error("\n\nstartrun_parameterfile: stepState must be > 0\n");

	if ( cmd.stepAvgRhoAxes == 0 )
		error("\n\nstartrun_parameterfile: stepAvgRhoAxes must be != 0\n");
	if ( cmd.stepAvgNFrecAxes == 0 )
		error("\n\nstartrun_parameterfile: stepAvgNFrecAxes must be != 0\n");
	if ( cmd.stepAvgVel == 0 )
		error("\n\nstartrun_parameterfile: stepAvgVel must be != 0\n");
	if ( cmd.stepAvgRdf == 0 )
		error("\n\nstartrun_parameterfile: stepAvgRdf must be != 0\n");
	if ( cmd.stepAvgPressAxes == 0 )
		error("\n\nstartrun_parameterfile: stepAvgPressAxes must be != 0\n");
	if ( cmd.stepAvgChemPot == 0 )
		error("\n\nstartrun_parameterfile: stepAvgChemPot must be != 0\n");
	if ( cmd.stepAvgAcf == 0 )
		error("\n\nstartrun_parameterfile: stepAvgAcf must be != 0\n");
	if ( cmd.stepAvgDiffuse == 0 )
		error("\n\nstartrun_parameterfile: stepAvgDiffuse must be != 0\n");
	if ( cmd.stepAvgCorr == 0 )
		error("\n\nstartrun_parameterfile: stepAvgCorr must be != 0\n");

	if ( gd.dtime == 0. )
		error("\n\nstartrun_parameterfile: dtime must be != 0\n");

	if (gd.nbody1==0 && gd.nbody2==0)
		error("\nstartrun_parameterfile: %s",
			"nbody1 and nbody2 must not be zero at the same time\n\n");
	if (gd.nbody1 < 0)
		error("startrun_parameterfile: absurd value for nbody1\n");
	if (gd.nbody2 < 0)
		error("startrun_parameterfile: absurd value for nbody2\n");
	if (gd.nbc1 < 0)
		error("startrun_parameterfile: absurd value for nbc1\n");
	if (gd.nbc2 < 0)
		error("startrun_parameterfile: absurd value for nbc2\n");
#ifdef THREEDIM
	if (cmd.icModel==6)
		if (gd.nbc1+gd.nbc2!=4)
			error("\n\nCheckParameters: icModel=6 : nbodyprop=n1/n2 and n1+n2 must be 4!: %d %d\n",
				gd.nbc1, gd.nbc2);
#else
	if (cmd.icModel==6)
		if (gd.nbc1+gd.nbc2!=2)
			error("\n\nCheckParameters: icModel=6 : nbodyprop=n1/n2 and n1+n2 must be 2!: %d %d\n",
				gd.nbc1, gd.nbc2);
#endif

#ifdef THREEDIM
	if (gd.Lx <= 0.)
		error("startrun_parameterfile: absurd value for Lx\n");
	if (gd.Ly <= 0.)
		error("startrun_parameterfile: absurd value for Ly\n");
#else
	if (gd.Lx <= 0.)
		error("startrun_parameterfile: absurd value for Lx\n");
#endif

	if (gd.mass1==0. && gd.mass2==0.)
		error("\nstartrun_parameterfile: %s",
				"mass1 and mass2 must not be zero at the same time\n\n");
	if (gd.mass1 < 0.)
		error("startrun_parameterfile: absurd value for mass1\n");
	if (gd.mass2 < 0.)
		error("startrun_parameterfile: absurd value for mass2\n");

	if (cmd.potType == POTTYPE_FILE && strnull(cmd.fnamePot))
		error("CheckParameters: You should give potential filename\n");
}

local void Header_to_Global(void)						// CHECK 2D --- OK!!!
{
	gd.nbody		= hdr.nbody;
	gd.nbody1		= hdr.nbody1;
	gd.nbody2		= hdr.nbody2;
	gd.tnow			= hdr.tnow;
	cmd.temperature	= hdr.temperature;
	cmd.density		= hdr.density;
	gd.mass1		= hdr.mass1;
	gd.mass2		= hdr.mass2;
	gd.Lx			= hdr.Lx;
	gd.Ly			= hdr.Ly;
#ifdef THREEDIM
	gd.Lz			= hdr.Lz;
#endif
	cmd.eps11		= hdr.eps11;
	cmd.eps12		= hdr.eps12;
	cmd.eps22		= hdr.eps22;
	cmd.sigma11		= hdr.sigma11;
	cmd.sigma12		= hdr.sigma12;
	cmd.sigma22		= hdr.sigma22;
	cmd.Rcut11		= hdr.Rcut11;
	cmd.Rcut12		= hdr.Rcut12;
	cmd.Rcut22		= hdr.Rcut22;
}

local void Compute_Parameters(void)						// CHECK 2D --- OK!!!
{
    bodyptr p;
	real rri, rri3, rtmp;
	int k;
	real Rcut, fphi, ssq, fa;

	fprintf(gd.outlog,"\n\nCompute_Parameters ...");

	gdforce.Rcut11Max = MAX(cmd.Rcut11,cmd.Rcut12);
	gdforce.Rcut22Max = MAX(cmd.Rcut22,cmd.Rcut12);
	gd.Rcut11Min = MIN(cmd.Rcut11,cmd.Rcut12);
	gd.Rcut22Min = MIN(cmd.Rcut22,cmd.Rcut12);

	gd.RcutAllMax = MAX(gdforce.Rcut11Max,gdforce.Rcut22Max);
	gd.RcutAllMin = MIN(gd.Rcut11Min,gd.Rcut22Min);

	gdforce.RcutSq11 = cmd.Rcut11*cmd.Rcut11;
	gdforce.RcutSq12 = cmd.Rcut12*cmd.Rcut12;
	gdforce.RcutSq22 = cmd.Rcut22*cmd.Rcut22;
	
	gd.RcutSq11Max = gdforce.Rcut11Max*gdforce.Rcut11Max;
	gd.RcutSq11Min = gd.Rcut11Min*gd.Rcut11Min;
	gd.RcutSq22Max = gdforce.Rcut22Max*gdforce.Rcut22Max;
	gd.RcutSq22Min = gd.Rcut22Min*gd.Rcut22Min;
	gd.RcutSqAllMax = gd.RcutAllMax*gd.RcutAllMax;
	gd.RcutSqAllMin = gd.RcutAllMin*gd.RcutAllMin;

	Rcut = gd.RcutAllMax;

//	gd.stepAvg = (int) (gd.dtout/gd.dtime);
//	gd.stepSnap = gd.stepAvg;
//	gd.stepSnapInit = 0;	// borrar cuando ya este seguro de que no sirven...

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	gd.kinEnergySave = 0.0;

	fphi = 4.0*gd.eps;
	ssq = gd.sigma*gd.sigma;
	fa = 48.0*gd.eps/ssq;

	gdforce.potType = cmd.potType;
	PotentialParameters(gdforce.potType, 
		cmd.sigma11, cmd.sigma12, cmd.sigma22,
		cmd.eps11, cmd.eps12, cmd.eps22, 
		cmd.Rcut11, cmd.Rcut12, cmd.Rcut22, &gdforce);
	if (cmd.potType==POTTYPE_FILE && 
		gd.RcutAllMin>rPos(forcepottab+nforcepot-1))
		error("\n\nCompute_Parameters: RcutAllMin > rPotMax [%g, %g]\n\n",
			gd.RcutAllMin,rPos(forcepottab+nforcepot-1));

	MULVS(gdforce.cells, gdforce.Box, 1.0/Rcut);
	AllocMem(gdforce.cellList, VProd (gdforce.cells)
				+ gd.nbody, int);
#if (NDIM==3)
	if (scanopt(cmd.forcecalc_method,"cells") 
			&& (gdforce.cells[0]==0||gdforce.cells[1]==0||gdforce.cells[2]==0))
#else
	if (scanopt(cmd.forcecalc_method,"cells") 
			&& (gdforce.cells[0]==0||gdforce.cells[1]==0))
#endif
		error("\n\nCompute_Parameters: invalid number of cells\n\n");

#if (NDIM==2)
	gd.LBoxMin=MIN(gdforce.Box[0],gdforce.Box[1]);
#else
	rtmp=MIN(gdforce.Box[0],gdforce.Box[1]);
	gd.LBoxMin=MIN(gdforce.Box[2],rtmp);
#endif

	gd.VelMax=0.;
	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		rtmp = 0.0;
		DO_COORD(k)
			rtmp += rsqr(Vel(p)[k]);
		if (gd.VelMax<rtmp) gd.VelMax=rtmp;
	}
	gd.VelMax=rsqrt(gd.VelMax);
	gd.TimeMin = gd.LBoxMin/gd.VelMax;

	gd.lAvg = 0.; DO_COORD(k) gd.lAvg += gdforce.Box[k];
	gd.lAvg = (gd.lAvg/NDIM) / rpow(gd.nbody, 1./(real)NDIM);
	gd.dtcrit = gd.lAvg/gd.vMag;
	rri=ssq/rsqr(gd.lAvg); rri3=rri*rri*rri;

	gd.fAvg = fa*rri3*(rri3-0.5)*rri*gd.lAvg;			// USAR LA RUTINA GENERAL...
	gd.phiAvg = fphi*(rri3-1.0)*rri3;

	gdforce.computeTransport = cmd.computeTransport;

	fprintf(gd.outlog,"done.\n");
}

local void Compute_nbody(void)							// CHECK 2D --- OK!!!
{
	fprintf(gd.outlog,"\n\nCompute_nbody ...");
	if (cmd.icModel==3 || cmd.icModel==6) {
//	if (cmd.icModel==3) {
		gd.nbody1 = VProd(gd.numUcell);
		gd.nbody2 = gd.nbody1;
//		if (cmd.icModel==6) {
//			gd.nbc1 *= 4;
//			gd.nbc2 *= 4;
//		}
		gd.nbody1 = gd.nbody1*gd.nbc1;
		gd.nbody2 = gd.nbody2*gd.nbc2;
	}
	gd.nbody=gd.nbody1+gd.nbody2;

	fprintf(gd.outlog,"\nInitial nbody : %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);
	fprintf(gd.outlog,"\ndone.\n");
}

local void Compute_Box_Units(void)						// CHECK 2D --- OK!!!
{
	fprintf(gd.outlog,"\n\nCompute_Box_Units ...");
	gd.eps = 1.0; gd.sigma = 1.0;	gd.mass = 1.0;
	gd.kB = 1.0;

	if (cmd.icModel==3 || cmd.icModel==6) {
//	if (cmd.icModel==3) {
#if (NDIM==3)
		if (strnull(cmd.icfile) && strnull(cmd.restorefile)) {
			gd.Lx=gd.numUcell[0]*rpow((gd.nbc1*gd.mass1
						+gd.nbc2*gd.mass2)/cmd.density, 1./3.);
			gd.Ly=gd.numUcell[1]*rpow((gd.nbc1*gd.mass1
						+gd.nbc2*gd.mass2)/cmd.density, 1./3.);
			gd.Lz=gd.numUcell[2]*rpow((gd.nbc1*gd.mass1
						+gd.nbc2*gd.mass2)/cmd.density, 1./3.);
			gdforce.Box[0]=gd.Lx; gdforce.Box[1]=gd.Ly; gdforce.Box[2]=gd.Lz;
		} else
			gdforce.Box[0]=gd.Lx; gdforce.Box[1]=gd.Ly; gdforce.Box[2]=gd.Lz;
#else
		if (strnull(cmd.icfile) && strnull(cmd.restorefile)) {
			gd.Lx=gd.numUcell[0]*rsqrt((gd.nbc1*gd.mass1
						+gd.nbc2*gd.mass2)/cmd.density);
			gd.Ly=gd.numUcell[1]*rsqrt((gd.nbc1*gd.mass1
						+gd.nbc2*gd.mass2)/cmd.density);
			gdforce.Box[0]=gd.Lx; gdforce.Box[1]=gd.Ly;
		} else
			gdforce.Box[0]=gd.Lx; gdforce.Box[1]=gd.Ly;
#endif
	} else {
#if (NDIM==3)
		if (strnull(cmd.icfile) && strnull(cmd.restorefile)) {
			gd.Lz = ( (real) gd.nbody1 * gd.mass1 + 
					(real) gd.nbody2 * gd.mass2 )/(gd.Lx*gd.Ly*cmd.density);
			gdforce.Box[0]=gd.Lx; gdforce.Box[1]=gd.Ly; gdforce.Box[2]=gd.Lz;
		} else
			gdforce.Box[0]=gd.Lx; gdforce.Box[1]=gd.Ly; gdforce.Box[2]=gd.Lz;
#else
		if (strnull(cmd.icfile) && strnull(cmd.restorefile)) {
			gd.Ly = ( (real) gd.nbody1 * gd.mass1 + 
					(real) gd.nbody2 * gd.mass2 )/(gd.Lx*cmd.density);
			gdforce.Box[0] = gd.Lx; gdforce.Box[1] = gd.Ly;
		} else
			gdforce.Box[0] = gd.Lx; gdforce.Box[1] = gd.Ly;
#endif
	}

	fprintf(gd.outlog,"done.\n");
}

local void Print_Parameters(void)						// CHECK 2D --- OK!!!
{
	fprintf(gd.outlog,"\nParameters:\n");
#if (NDIM==3)
	fprintf(gd.outlog,"\nBox: %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\ncells: %d %d %d",
		gdforce.cells[0],gdforce.cells[1],gdforce.cells[2]);
#else
#if (NDIM==2)
	fprintf(gd.outlog,"\nBox: %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\ncells: %d %d",gdforce.cells[0],gdforce.cells[1]);
#endif
#endif
	fprintf(gd.outlog,"\neps= %g",gd.eps);
	fprintf(gd.outlog,"\nsigma= %g",gd.sigma);
	fprintf(gd.outlog,"\ntemperature= %g",cmd.temperature);
	fprintf(gd.outlog,"\nmass= %g",gd.mass);
	fprintf(gd.outlog,"\ndensity= %g",cmd.density);

	fprintf(gd.outlog,"\nnum. of cells allocated (cellList) = %d",
		VProd (gdforce.cells)+gd.nbody);

	fprintf(gd.outlog,"\nstepAvg= %d\n\n",cmd.stepAvg);
	fprintf(gd.outlog,"\nlAvg= %g",gd.lAvg);
	fprintf(gd.outlog,"\ndtcrit= %g",gd.dtcrit);
	fprintf(gd.outlog,"\nfAvg= %g",gd.fAvg);
	fprintf(gd.outlog,"\nphiAvg= %g",gd.phiAvg);

	fprintf(gd.outlog,"\nRcut= %g",cmd.Rcut11);
	fprintf(gd.outlog,"\nRcutSq= %g",gdforce.RcutSq11);
	fprintf(gd.outlog,"\nRcut= %g",cmd.Rcut12);
	fprintf(gd.outlog,"\nRcutSq= %g",gdforce.RcutSq12);
	fprintf(gd.outlog,"\nRcut= %g",cmd.Rcut22);
	fprintf(gd.outlog,"\nRcutSq= %g",gdforce.RcutSq22);
	fprintf(gd.outlog,"\neps11= %g",cmd.eps11);
	fprintf(gd.outlog,"\neps12= %g",cmd.eps12);
	fprintf(gd.outlog,"\neps22= %g",cmd.eps22);
	fprintf(gd.outlog,"\nsigma11= %g",cmd.sigma11);
	fprintf(gd.outlog,"\nsigma12= %g",cmd.sigma12);
	fprintf(gd.outlog,"\nsigma22= %g",cmd.sigma22);
	fprintf(gd.outlog,"\nmass1= %g",gd.mass1);
	fprintf(gd.outlog,"\nmass2= %g",gd.mass2);

	fprintf(gd.outlog,"\nfphi11= %g",gdforce.fphi11);
	fprintf(gd.outlog,"\nfphi12= %g",gdforce.fphi12);
	fprintf(gd.outlog,"\nfphi22= %g",gdforce.fphi22);
	fprintf(gd.outlog,"\nfa11= %g",gdforce.fa11);
	fprintf(gd.outlog,"\nfa12= %g",gdforce.fa12);
	fprintf(gd.outlog,"\nfa22= %g",gdforce.fa22);

	fprintf(gd.outlog,"\n\nLBoxMin, VelMax : %g %g",gd.LBoxMin,gd.VelMax);
	fprintf(gd.outlog,"\nTimeMin, dtime : %g %g",gd.TimeMin, gd.dtime);
	fprintf(gd.outlog,"\ndtout, dtoutinfo : %g %g",gd.dtout, gd.dtoutinfo);

	fflush(gd.outlog);
}

/*
#define NULLMETHOD	0
#define TESTDATA1	1
#define TESTDATA2	2
#define TESTDATA3	3
#define TESTDATA4	4
#define TESTDATA5	5
#define TESTDATA6	6

local void testdata(void)								// CHECK 2D --- OK!!!
{
    switch(cmd.icModel) {
        case TESTDATA1:
            testdata_sp_01(); break;
        case TESTDATA2:
            testdata_sp_02(); break;
        case TESTDATA3:
            testdata_sp_03(); break;
        case TESTDATA4:
            testdata_sp_04(); break;
        case TESTDATA5:
            testdata_sp_05(); break;
        case TESTDATA6:
            testdata_sp_06(); break;
        case NULLMETHOD:
            printf("\n\tdefault IC ...\n");
            testdata_sp_02(); break;
        default:
            printf("\n\tUnknown IC method...");
            printf("\n\tdefault IC method...\n"); 
            testdata_sp_02(); break;
    }
}

#undef NULLMETHOD
#undef TESTDATA1
#undef TESTDATA2
#undef TESTDATA3
#undef TESTDATA4
#undef TESTDATA5
#undef TESTDATA6

local void testdata_sp_01(void)
{
    vector vcm, b, l;
	bodyptr p;
	real ls, Ekin, tmass, velsq;
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real lmin, lmax;
	vector dr;
    double cpustart;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_01: Parallelepiped : smaller parallelepipeds\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [1] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2,	gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

#if (NDIM==3)
	ls=rpow(gd.nbody,1./3.);
#else
	ls=rsqrt(gd.nbody);
#endif
	DIVVS(l,gdforce.Box,ls);

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g %g",
		l[0],l[1],l[2],ls);
#else
	fprintf(stdout,"\nSimulation box : %g %g",	gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g",
		l[0],l[1],ls);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

#if (NDIM==3)
	b[2]=0.;
	while(b[2]<gdforce.Box[2]) {
#endif
		b[1]=0.;
		while(b[1]<gdforce.Box[1]) {
			b[0]=0.;
			while(b[0]<gdforce.Box[0]) {
				++ibox;
				++i; ++i1;
				if (i1<=gd.nbody1) {
					DO_COORD(k)
						Pos(p)[k] = b[k]+0.001*l[k];
					Id(p) = p-bodytab+1;
					Type(p) = BODY1;
					Mass(p) = gd.mass1;
					DO_COORD(k)
						Vel(p)[k] = grandom(0.0,1.0);
					tmass += Mass(p);
					ADDMULVS(vcm, Vel(p), Mass(p));
					++p;
				} else { 
					--i1; 
					--i;
					++i; 
					++i2;
					if (i2<=gd.nbody2) {
						DO_COORD(k)
							Pos(p)[k] = b[k]+0.001*l[k];
						Id(p) = p-bodytab+1;
						Type(p) = BODY2;
						Mass(p) = gd.mass2;
						DO_COORD(k)
							Vel(p)[k] = grandom(0.0,1.0);
						tmass += Mass(p);
						ADDMULVS(vcm, Vel(p), Mass(p));
						++p;
					} else { 
						--i2; 
						--i; 
					}
				}

				b[0] += l[0];
			}
			b[1] += l[1];
		}
#if (NDIM==3)
		b[2] += l[2];
	}
#endif

	fprintf(gd.outlog,
		"\nTotal number of particles created: %d %d in %d invoqued boxes\n",
		p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin
	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	DO_BODY(p, bodytab, bodytab+gd.nbody)
		VVSAdd(Pos(p), -0.5, gdforce.Box);		// Box center is at (0,0,0)

	Diagnose();

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 01
*/

/*
local void testdata_sp_02(void)
{
    vector vcm, b, l;
	bodyptr p;
	real ls, Ekin, tmass, velsq;
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real lmin, lmax;
	real vol;
    double cpustart;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_02: Parallelepiped : smaller cubic boxes\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [2] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

#if (NDIM==3)
	vol = gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]/gd.nbody;
	ls=rpow(vol,1./3.);
	l[0]=l[1]=l[2]=ls;
#else
	vol = gdforce.Box[0]*gdforce.Box[1]/gd.nbody
	ls=rsqrt(vol);
	l[0]=l[1]=ls;
#endif

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g %g",
		l[0],l[1],l[2],ls);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g",
		l[0],l[1],ls);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

#if (NDIM==3)
	b[2]=0.;
	while(b[2]<gdforce.Box[2]) {
#endif
		b[1]=0.;
		while(b[1]<gdforce.Box[1]) {
			b[0]=0.;
			while(b[0]<gdforce.Box[0]) {
				++ibox;
				++i; ++i1;
				if (i1<=gd.nbody1) {
					DO_COORD(k) {
						Pos(p)[k] = b[k]+0.5*l[k];
						if(Pos(p)[k]>gdforce.Box[k])
							Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
					}
					Id(p) = p-bodytab+1;
					Type(p) = BODY1;
					Mass(p) = gd.mass1;
					DO_COORD(k)
						Vel(p)[k] = grandom(0.0,1.0);
					tmass += Mass(p);
					ADDMULVS(vcm, Vel(p), Mass(p));
					++p;
				} else { 
					--i1; 
					--i;
					++i; 
					++i2;
					if (i2<=gd.nbody2) {
						DO_COORD(k) {
							Pos(p)[k] = b[k]+0.5*l[k];
							if(Pos(p)[k]>gdforce.Box[k])
								Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
						}
						Id(p) = p-bodytab+1;
						Type(p) = BODY2;
						Mass(p) = gd.mass2;
						DO_COORD(k)
							Vel(p)[k] = grandom(0.0,1.0);
						tmass += Mass(p);
						ADDMULVS(vcm, Vel(p), Mass(p));
						++p;
					} else { 
						--i2; 
						--i; 
					}
				}

				b[0] += l[0];
			}
			b[1] += l[1];
		}
#if (NDIM==3)
		b[2] += l[2];
	}
#endif

	fprintf(gd.outlog,
		"\nTotal number of particles created: %d %d in %d invoqued boxes\n",
		p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	DO_BODY(p, bodytab, bodytab+gd.nbody)
		VVSAdd(Pos(p), -0.5, gdforce.Box);		// Box center is at (0,0,0)

	Diagnose();

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 02
*/
/*
local void testdata_sp_05(void)
{
    vector vcm, b, l;
	bodyptr p, q;
	real ls, Ekin, tmass, velsq; 
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real dist, distMin, lmin, lmax, drpq, drpq2;
	vector dr;
	int ido;
    double cpustart;
	int iran, ranmax=200;

    cpustart = cputime();                       
	fprintf(gd.outlog,
		"\ntestdata_sp_05: Parallelepiped : Spatially random positions\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [3] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

#if (NDIM==3)
	ls=rpow(gd.nbody,1./3.);
#else
	ls=rsqrt(gd.nbody);
#endif
	DIVVS(l,gdforce.Box,ls);

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
#endif
	lmin *= 1.0;   // Set appropriately. Smaller the scale faster the routine...
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

	for (i=1; i<=gd.nbody; i++) {
		++i1;
		if (i1<=gd.nbody1) {
			iran=0;
			do {
				DO_COORD(k)
					Pos(p)[k] = xrandom(-0.5*gdforce.Box[k],0.5*gdforce.Box[k]);
				distMin = lmax;
				for (q=bodytab; q<p; q++) {
					DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
					DO_COORD(k)
						dr[k]=dr[k] - 
							((real)(nint(dr[k]/gdforce.Box[k])))*gdforce.Box[k];
					DOTVP(drpq2, dr, dr);
					drpq = rsqrt(drpq2);
					distMin = MIN(distMin, drpq);
				}
				++iran; 
				if (iran>ranmax) 
					error("\ntestdata_sp_05: [1] body %d : %s",
						p-bodytab,"Maximum number of iter reached\n");
			} while (distMin < lmin);
			Id(p) = p-bodytab+1;
			Type(p) = BODY1;
			Mass(p) = gd.mass1;
			DO_COORD(k)
				Vel(p)[k] = grandom(0.0,1.0);
			tmass += Mass(p);
			ADDMULVS(vcm, Vel(p), Mass(p));
			++p;
		} else {
			--i1; ++i2;
			if (i2<=gd.nbody2) {
				iran=0;
				do {
					DO_COORD(k)
						Pos(p)[k] = 
							xrandom(-0.5*gdforce.Box[k], 0.5*gdforce.Box[k]);
					distMin = lmax;
					for (q=bodytab; q<p; q++) {
						DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
						DO_COORD(k)
							dr[k]=dr[k] - 
							((real)(nint(dr[k]/gdforce.Box[k])))*gdforce.Box[k];
						DOTVP(drpq2, dr, dr);
						drpq = rsqrt(drpq2);
						distMin = MIN(distMin, drpq);
					}
					++iran;
					if (iran>ranmax)
						error("\ntestdata_sp_05: [2] body %d : %s",
							p-bodytab,"Maximum number of iter reached\n");
				} while (distMin < lmin);
				Id(p) = p-bodytab+1;
				Type(p) = BODY2;
				Mass(p) = gd.mass2;
				DO_COORD(k)
					Vel(p)[k] = grandom(0.0,1.0);
				tmass += Mass(p);
				ADDMULVS(vcm, Vel(p), Mass(p));
				++p;
			} else {
				--i2;
			}
		}
	}

	fprintf(gd.outlog,
		"\nTotal number of particles created: %d %d in %d invoqued boxes\n",
		p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	Diagnose();

#if (NDIM==3)
		fprintf(gd.outlog,"\ndensity = %g %g\n",
			tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
		fprintf(gd.outlog,"\ndensity = %g %g\n",
			tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 05
*/

/*
local void testdata_sp_04(void)
{
    vector vcm, b, l;
	bodyptr p;
	real ls, Ekin, tmass, velsq, L;
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real lmin, lmax;
	real vol;
    double cpustart;

    cpustart = cputime();                       

	fprintf(gd.outlog,"\ntestdata_sp_04: Cubic box : smaller cubic boxes\n");

    strcpy(gd.model_comment, "Simple Cubic Box Model [4] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

	vol = (((real)gd.nbody1)*gd.mass1+((real)gd.nbody2)*gd.mass2)/cmd.density;
#if (NDIM==3)
	L=rpow(vol,1./3.);
	gdforce.Box[0]=gdforce.Box[1]=gdforce.Box[2]=L;
	gd.Lx=gd.Ly=gd.Lz=L;
	ls=L/rpow(gd.nbody,1./3.);
	l[0]=l[1]=l[2]=ls;
#else
	L=rsqrt(vol);
	gdforce.Box[0]=gdforce.Box[1]=L;
	gd.Lx=gd.Ly=L;
	ls=L/rsqrt(gd.nbody);
	l[0]=l[1]=ls;
#endif

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g %g",
		l[0],l[1],l[2],ls);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g",
		l[0],l[1],ls);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

#if (NDIM==3)
	b[2]=0.;
	while(b[2]<gdforce.Box[2]) {
#endif
		b[1]=0.;
		while(b[1]<gdforce.Box[1]) {
			b[0]=0.;
			while(b[0]<gdforce.Box[0]) {
				++ibox;
				++i; ++i1;
				if (i1<=gd.nbody1) {
					DO_COORD(k) {
						Pos(p)[k] = b[k]+0.5*l[k];
						if(Pos(p)[k]>gdforce.Box[k])
							Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
					}
					Id(p) = p-bodytab+1;
					Type(p) = BODY1;
					Mass(p) = gd.mass1;
					DO_COORD(k)
						Vel(p)[k] = grandom(0.0,1.0);
					tmass += Mass(p);
					ADDMULVS(vcm, Vel(p), Mass(p));
					++p;
				} else { 
					--i1; 
					--i;
					++i; 
					++i2;
					if (i2<=gd.nbody2) {
						DO_COORD(k) {
							Pos(p)[k] = b[k]+0.5*l[k];
							if(Pos(p)[k]>gdforce.Box[k])
								Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
						}
						Id(p) = p-bodytab+1;
						Type(p) = BODY2;
						Mass(p) = gd.mass2;
						DO_COORD(k)
							Vel(p)[k] = grandom(0.0,1.0);
						tmass += Mass(p);
						ADDMULVS(vcm, Vel(p), Mass(p));
						++p;
					} else { 
						--i2; 
						--i; 
					}
				}

				b[0] += l[0];
			}
			b[1] += l[1];
		}
#if (NDIM==3)
		b[2] += l[2];
	}
#endif

	fprintf(gd.outlog,"\n%s: %d %d in %d invoqued boxes\n",
		"Total number of particles created",p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	DO_BODY(p, bodytab, bodytab+gd.nbody)
		VVSAdd(Pos(p), -0.5, gdforce.Box);		// Box center is at (0,0,0)

	Diagnose();

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 04
*/

/*
local void testdata_sp_03(void)
{
	vector c, gap, vSum, l, sc, sgap, ctmp;
	vectorI scell;
	int nx, ny, nz, i, i1, i2, si1, si2, j, nxs, nys, nzs;
    double cpustart;
	real tmass, lmin, lmax, ls;
	bodyptr p;
	real velMag;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_03: Simple parallelepiped box : smaller cubic boxes\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [3] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

	VDiv (gap, gdforce.Box, gd.numUcell);

#if (NDIM==3)
	ls=rpow((gd.nbc1*gd.mass1+gd.nbc2*gd.mass2)/cmd.density, 1./3.);
	l[0]=l[1]=l[2]=ls;
#else
	ls=gd.numUcell[0]*rsqrt((gd.nbc1*gd.mass1+gd.nbc2*gd.mass2)/cmd.density);
	l[0]=l[1]=ls;
#endif

	j=1;
	while (gd.nbc1+gd.nbc2 > Cube(j)) ++j;
#if (NDIM==3)
	scell[0]=scell[1]=scell[2]=j;
#else
	scell[0]=scell[1]=j;
#endif

	VDiv (sgap, l, scell);
#if (NDIM==3)
	fprintf(stdout,"\nsGap : %g %g %g", sgap[0], sgap[1], sgap[2]);
#else
	fprintf(stdout,"\nsGap : %g %g",sgap[0],sgap[1]);
#endif

	lmin = MIN(sgap[0],sgap[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(sgap[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nUnit cell : l : %g %g %g",l[0],l[1],l[2]);
	fprintf(stdout,"\nGap : %g %g %g", gap[0], gap[1], gap[2]);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nUnit cell : l : %g %g",l[0],l[1]);
	fprintf(stdout,"\nGap : %g %g",gap[0],gap[1]);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

	fprintf(stdout,"\nNumber of bodies in %d smaller cubes : %d %d",
			Cube(j), gd.nbc1, gd.nbc2);

	p = bodytab; i=i1=i2=0;
	tmass = 0.;

#if (NDIM==3)
	for (nz = 0; nz < gd.numUcell[2]; nz++) {
#endif
		for (ny = 0; ny < gd.numUcell[1]; ny++) {
			for (nx = 0; nx < gd.numUcell[0]; nx++) {
				VSet (c, nx, ny, nz);
				VMul (c, c, gap);
				si1=si2=0;
#if (NDIM==3)
				for (nzs = 0; nzs < scell[2]; nzs++) {
#endif
					for (nys = 0; nys < scell[1]; nys++) {
						for (nxs = 0; nxs < scell[0]; nxs++) {
							VSet (sc, nxs + 0.5, nys + 0.5, nzs + 0.5);
							VMul (sc, sc, sgap);
							VAdd(ctmp, c, sc);
							VVSAdd(ctmp, -0.5, gdforce.Box);
							++i; ++si1;
							if (si1<=gd.nbc1) {
								SETV(Pos(p), ctmp);
								Id(p) = p-bodytab+1;
								Type(p) = BODY1;
								Mass(p) = gd.mass1;
								tmass += Mass(p);
								++p; ++i1;
							} else { 
								--si1; 
								--i;
								++i; 
								++si2;
								if (si2<=gd.nbc2) {
									SETV(Pos(p), ctmp);
									Id(p) = p-bodytab+1;
									Type(p) = BODY2;
									Mass(p) = gd.mass2;
									tmass += Mass(p);
									++p; ++i2;
								} else { 
									--si2; 
									--i; 
								}
							}
						}
					}
#if (NDIM==3)
				}
#endif
			}
		}
#if (NDIM==3)
	}
#endif

	velMag = sqrt (NDIM * (1. - 1. / gd.nbody) * cmd.temperature);
	VZero(vSum);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		VRand(Vel(p));
		VScale(Vel(p), velMag);
		VVAdd(vSum, Vel(p));
	}
	DO_BODY(p, bodytab, bodytab+gd.nbody) 
		VVSAdd (Vel(p), - 1. / gd.nbody, vSum);

	Diagnose();

	fprintf(gd.outlog,"\n%s: %d in %d invoqued boxes\n",
		"Total number of particles created",p-bodytab, 
		VProd(gd.numUcell));
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
}  // End icModel 03
*/

/*
local void testdata_sp_06(void)
{
	vector c, gap, vSum, l, sc, sgap, ctmp;
	vectorI scell;
	int nx, ny, nz, i, i1, i2, si1, si2, j, nxs, nys, nzs;
    double cpustart;
	real tmass, lmin, lmax, ls;
	bodyptr p;
	real velMag;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_03: Simple parallelepiped box : smaller cubic boxes (FCC)\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [6] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

	VDiv (gap, gdforce.Box, gd.numUcell);

#if (NDIM==3)
	ls=rpow((gd.nbc1*gd.mass1+gd.nbc2*gd.mass2)/cmd.density, 1./3.);
	l[0]=l[1]=l[2]=ls;
#else
	ls=gd.numUcell[0]*rsqrt((gd.nbc1*gd.mass1+gd.nbc2*gd.mass2)/cmd.density);
	l[0]=l[1]=ls;
#endif

	j=1;
	while (gd.nbc1+gd.nbc2 > Cube(j)) ++j;
#if (NDIM==3)
	scell[0]=scell[1]=scell[2]=j;
#else
	scell[0]=scell[1]=j;
#endif

	VDiv (sgap, l, scell);
#if (NDIM==3)
	fprintf(stdout,"\nsGap : %g %g %g", sgap[0], sgap[1], sgap[2]);
#else
	fprintf(stdout,"\nsGap : %g %g",sgap[0],sgap[1]);
#endif

	lmin = MIN(sgap[0],sgap[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(sgap[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nUnit cell : l : %g %g %g",l[0],l[1],l[2]);
	fprintf(stdout,"\nGap : %g %g %g", gap[0], gap[1], gap[2]);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nUnit cell : l : %g %g",l[0],l[1]);
	fprintf(stdout,"\nGap : %g %g",gap[0],gap[1]);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

	fprintf(stdout,"\nNumber of bodies in %d smaller cubes : %d %d",
			Cube(j), gd.nbc1, gd.nbc2);

	p = bodytab; i=i1=i2=0;
	tmass = 0.;

#if (NDIM==3)
	for (nz = 0; nz < gd.numUcell[2]; nz++) {
#endif
		for (ny = 0; ny < gd.numUcell[1]; ny++) {
			for (nx = 0; nx < gd.numUcell[0]; nx++) {
				VSet (c, nx+0.25, ny+0.25, nz+0.25);
				VMul (c, c, gap);
				VVSAdd (c, -0.5, gdforce.Box);
				if (i1<=gd.nbody1) {
					for (j = 0; j < 4; j ++) {
						SETV(Pos(p), c);
						if (j != 3) {
							if (j != 0) Pos(p)[0] += 0.5 * gap[0];
							if (j != 1) Pos(p)[1] += 0.5 * gap[1];
							if (j != 2) Pos(p)[2] += 0.5 * gap[2];
						}
						Id(p) = p-bodytab+1;
						Type(p) = BODY1;
						Mass(p) = gd.mass1;
						tmass += Mass(p);
						++p; ++i1;
					}
				}
//				si1=si2=0;
*/
/*
#if (NDIM==3)
				for (nzs = 0; nzs < scell[2]; nzs++) {
#endif
					for (nys = 0; nys < scell[1]; nys++) {
						for (nxs = 0; nxs < scell[0]; nxs++) {
							VSet (sc, nxs + 0.5, nys + 0.5, nzs + 0.5); // OK
							VMul (sc, sc, sgap);  // OK
							VAdd(ctmp, c, sc); 
							VVSAdd(ctmp, -0.5, gdforce.Box); // OK
							++i; ++si1;
							if (si1<=gd.nbc1) {
								SETV(Pos(p), ctmp);
								Id(p) = p-bodytab+1;
								Type(p) = BODY1;
								Mass(p) = gd.mass1;
								tmass += Mass(p);
								++p; ++i1;
							} else { 
								--si1; 
								--i;
								++i;					// Solo una ESPECIE...
								++si2;
								if (si2<=gd.nbc2) {
									SETV(Pos(p), ctmp);
									Id(p) = p-bodytab+1;
									Type(p) = BODY2;
									Mass(p) = gd.mass2;
									tmass += Mass(p);
									++p; ++i2;
								} else { 
									--si2; 
									--i; 
								}
							}
						}
					}
#if (NDIM==3)
				}
#endif
*/
/*			}
		}
#if (NDIM==3)
	}
#endif

	velMag = sqrt (NDIM * (1. - 1. / gd.nbody) * cmd.temperature);
	VZero(vSum);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		VRand(Vel(p));
		VScale(Vel(p), velMag);
		VVAdd(vSum, Vel(p));
	}
	DO_BODY(p, bodytab, bodytab+gd.nbody) 
		VVSAdd (Vel(p), - 1. / gd.nbody, vSum);

	Diagnose();

	fprintf(gd.outlog,"\n%s: %d in %d invoqued boxes\n",
		"Total number of particles created",p-bodytab, 
		VProd(gd.numUcell));
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 6
*/


local void ReadParameterFile(char *fname)				// CHECK 2D --- OK!!!
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

  nt=0;

	SPName(cmd.icfile,"icfile",100);
	SPName(cmd.icfilefmt,"icfilefmt",100);
	SPName(cmd.snapoutfile,"snapout",100);
	SPName(cmd.snapoutfilefmt,"snapoutfmt",100);

	IPName(cmd.icModel,"icModel");
	IPName(cmd.intMethod,"intMethod");
	SPName(cmd.dtimestr,"dtime",64);

	RPName(cmd.temperature,"temperature");
	BPName(cmd.adjustTemperature,"adjustTemperature");
	IPName(cmd.stepAdjustTemperature,"stepAdjustTemperature");
	BPName(cmd.adjustCenterOfMass,"adjustCenterOfMass");
	RPName(cmd.density,"density");

	IPName(cmd.stepEquil,"stepEquil");
	IPName(cmd.stepAvg,"stepAvg");

// borrar cuando ya este seguro de que no sirven...
//	BPName(cmd.printSnap,"printSnap");
//	IPName(cmd.stepSnap,"stepSnap");

	BPName(cmd.computeRhoAxes,"computeRhoAxes");
	IPName(cmd.stepAvgRhoAxes,"stepAvgRhoAxes");
	IPName(cmd.stepRhoAxes,"stepRhoAxes");
	IPName(cmd.sizeHistRhoAxes,"sizeHistRhoAxes");

	BPName(cmd.computeNFrecAxes,"computeNFrecAxes");
	IPName(cmd.stepAvgNFrecAxes,"stepAvgNFrecAxes");
	IPName(cmd.stepNFrecAxes,"stepNFrecAxes");
	IPName(cmd.sizeHistNFrecAxes,"sizeHistNFrecAxes");

	BPName(cmd.computeVelDist,"computeVelDist");
	IPName(cmd.stepAvgVel,"stepAvgVel");
	IPName(cmd.stepVel,"stepVel");
	IPName(cmd.sizeHistVel,"sizeHistVel");
	RPName(cmd.rangeVel,"rangeVel");

	BPName(cmd.computeRdf,"computeRdf");
	IPName(cmd.stepAvgRdf,"stepAvgRdf");
	IPName(cmd.stepRdf,"stepRdf");
	IPName(cmd.sizeHistRdf,"sizeHistRdf");
	RPName(cmd.rangeRdf,"rangeRdf");

	BPName(cmd.computePressAxes,"computePressAxes");
	IPName(cmd.stepAvgPressAxes,"stepAvgPressAxes");
	IPName(cmd.stepPressAxes,"stepPressAxes");
	IPName(cmd.sizeHistPressAxes,"sizeHistPressAxes");

	BPName(cmd.computeChemPot,"computeChemPot");
	IPName(cmd.stepAvgChemPot,"stepAvgChemPot");
	IPName(cmd.stepChemPot,"stepChemPot");
	IPName(cmd.sizeHistChemPot,"sizeHistChemPot");
	IPName(cmd.numTestBodies,"numTestBodies");

//	BPName(cmd.computeTransport,"computeDiffusion");		// Linea original parece
	BPName(cmd.computeDiffusion,"computeDiffusion");	// que hubo error CHECAR!!!!!
	IPName(cmd.stepDiffuse,"stepDiffuse");
	IPName(cmd.stepAvgDiffuse,"stepAvgDiffuse");
	IPName(cmd.nBuffDiffuse,"nBuffDiffuse");
	IPName(cmd.nValDiffuse,"nValDiffuse");

	BPName(cmd.computeVelAcf,"computeVelAcf");

	BPName(cmd.computeBulkViscosity,"computeBulkViscosity");

	BPName(cmd.computeTransport,"computeTransport");
	IPName(cmd.stepAcf,"stepAcf");
	IPName(cmd.stepAvgAcf,"stepAvgAcf");
	IPName(cmd.nBuffAcf,"nBuffAcf");
	IPName(cmd.nValAcf,"nValAcf");

	BPName(cmd.computeSTCorr,"computeSTCorr");
	IPName(cmd.stepCorr,"stepCorr");
	IPName(cmd.stepAvgCorr,"stepAvgCorr");
	IPName(cmd.nBuffCorr,"nBuffCorr");
	IPName(cmd.nValCorr,"nValCorr");
	IPName(cmd.nFunCorr,"nFunCorr");

	RPName(cmd.lattCorr_kx,"lattCorr_kx");
	RPName(cmd.lattCorr_ky,"lattCorr_ky");
#ifdef THREEDIM
	RPName(cmd.lattCorr_kz,"lattCorr_kz");
#endif

	SPName(cmd.options,"options",100);

	SPName(cmd.forcecalc_method,"forcecalc_method",100);
	IPName(cmd.potType,"potentialType");
	SPName(cmd.fnamePot,"fnamePot",100);

	RPName(cmd.tstop,"tstop");
	SPName(cmd.dtoutstr,"dtout",64);
	SPName(cmd.dtoutinfostr,"dtoutinfo",64);

	SPName(cmd.unitCells,"unitCells",100);
	SPName(cmd.nbodyprop,"nbodyprop",100);
	SPName(cmd.massprop,"massprop",100);
#if (NDIM==3)
	SPName(cmd.LxLyprop,"LxLyprop",100);
#else
	SPName(cmd.Lx,"Lx",100);
#endif

	RPName(cmd.eps11,"eps11");
	RPName(cmd.eps12,"eps12");
	RPName(cmd.eps22,"eps22");
	RPName(cmd.sigma11,"sigma11");
	RPName(cmd.sigma12,"sigma12");
	RPName(cmd.sigma22,"sigma22");
	RPName(cmd.Rcut11,"Rcut11");
	RPName(cmd.Rcut12,"Rcut12");
	RPName(cmd.Rcut22,"Rcut22");

	IPName(cmd.seed,"seed");

	SPName(cmd.statefile,"statefile",100);
	IPName(cmd.stepState,"stepState");

	SPName(cmd.restorefile,"restorefile",100);

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
						if (strchr("tTyY1", *buf2) != NULL)
                            *((bool*)addr[j])=TRUE;
                        else 
                            if (strchr("fFnN0", *buf2) != NULL)
                                *((bool*)addr[j])=FALSE;
							else
                                error("ReadParameterFile: %s=%s not bool\n", 
									buf1, buf2);
						break;
                }
			} else {
                error("ReadParameterFile: Error in file %s: Tag '%s' %s",
					fname,buf1,"not allowed or multiple defined.\n");
            }
        }
        fclose(fd);
    } else {
        error("Parameter file %s not found.\n", fname);
    }

    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            error("ReadParameterFile: %s '%s'\n\t\t %s '%s'\n\n",
				"Error. I miss a value for tag",tag[i],
				"in parameter file",fname);
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

local void PrintParameterFile(char *fname)				// CHECK 2D --- OK!!!
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
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",
			gd.headline1,gd.headline2,"%", gd.headline3);
        fprintf(fdout,"%s\n%s\n",
		"%-------------------------------------------------------------------",
		"%");
        fprintf(fdout,FMTT,"forcecalc_method",cmd.forcecalc_method);
        fprintf(fdout,FMTI,"potentialType",cmd.potType);
        fprintf(fdout,FMTT,"fnamePot",cmd.fnamePot);
        fprintf(fdout,FMTT,"icfile",cmd.icfile);
        fprintf(fdout,FMTT,"icfilefmt",cmd.icfilefmt);
        fprintf(fdout,FMTT,"snapout",cmd.snapoutfile);
        fprintf(fdout,FMTT,"snapoutfmt",cmd.snapoutfilefmt);
        fprintf(fdout,FMTT,"statefile",cmd.statefile);
        fprintf(fdout,FMTI,"stepState",cmd.stepState);

        fprintf(fdout,FMTT,"restorefile",cmd.restorefile);
        fprintf(fdout,FMTR,"temperature",cmd.temperature);
        fprintf(fdout,FMTT,"adjustTemperature",
				cmd.adjustTemperature ? "true" : "false");
        fprintf(fdout,FMTI,"stepAdjustTemperature",cmd.stepAdjustTemperature);
        fprintf(fdout,FMTT,"adjustCenterOfMass",
				cmd.adjustCenterOfMass ? "true" : "false");
        fprintf(fdout,FMTR,"density",cmd.density);
        fprintf(fdout,FMTI,"stepEquil",cmd.stepEquil);
        fprintf(fdout,FMTI,"stepAvg",cmd.stepAvg);

// borrar cuando ya este seguro de que no sirven...
//        fprintf(fdout,FMTT,"printSnap",
//				cmd.printSnap ? "true" : "false");
//        fprintf(fdout,FMTI,"stepSnap",cmd.stepSnap);

        fprintf(fdout,FMTT,"computeRhoAxes",
				cmd.computeRhoAxes ? "true" : "false");
        fprintf(fdout,FMTI,"stepAvgRhoAxes",cmd.stepAvgRhoAxes);
        fprintf(fdout,FMTI,"stepRhoAxes",cmd.stepRhoAxes);
        fprintf(fdout,FMTI,"sizeHistRhoAxes",cmd.sizeHistRhoAxes);

        fprintf(fdout,FMTT,"computeNFrecAxes",
				cmd.computeNFrecAxes ? "true" : "false");
        fprintf(fdout,FMTI,"stepAvgNFrecAxes",cmd.stepAvgNFrecAxes);
        fprintf(fdout,FMTI,"stepNFrecAxes",cmd.stepNFrecAxes);
        fprintf(fdout,FMTI,"sizeHistNFrecAxes",cmd.sizeHistNFrecAxes);

        fprintf(fdout,FMTT,"computeVelDist",
				cmd.computeVelDist ? "true" : "false");
        fprintf(fdout,FMTI,"stepAvgVel",cmd.stepAvgVel);
        fprintf(fdout,FMTI,"stepVel",cmd.stepVel);
        fprintf(fdout,FMTI,"sizeHistVel",cmd.sizeHistVel);
        fprintf(fdout,FMTR,"rangeVel",cmd.rangeVel);

        fprintf(fdout,FMTT,"computeRdf",
				cmd.computeRdf ? "true" : "false");
        fprintf(fdout,FMTI,"stepAvgRdf",cmd.stepAvgRdf);
        fprintf(fdout,FMTI,"stepRdf",cmd.stepRdf);
        fprintf(fdout,FMTI,"sizeHistRdf",cmd.sizeHistRdf);
        fprintf(fdout,FMTR,"rangeRdf",cmd.rangeRdf);

        fprintf(fdout,FMTT,"computePressAxes",
				cmd.computePressAxes ? "true" : "false");
        fprintf(fdout,FMTI,"stepAvgPressAxes",cmd.stepAvgPressAxes);
        fprintf(fdout,FMTI,"stepPressAxes",cmd.stepPressAxes);
        fprintf(fdout,FMTI,"sizeHistPressAxes",cmd.sizeHistPressAxes);

        fprintf(fdout,FMTT,"computeChemPot",
				cmd.computeChemPot ? "true" : "false");
        fprintf(fdout,FMTI,"stepAvgChemPot",cmd.stepAvgChemPot);
        fprintf(fdout,FMTI,"stepChemPot",cmd.stepChemPot);
        fprintf(fdout,FMTI,"sizeHistChemPot",cmd.sizeHistChemPot);
        fprintf(fdout,FMTI,"numTestBodies",cmd.numTestBodies);

        fprintf(fdout,FMTT,"computeDiffusion",
				cmd.computeDiffusion ? "true" : "false");
        fprintf(fdout,FMTI,"stepDiffuse",cmd.stepDiffuse);
        fprintf(fdout,FMTI,"stepAvgDiffuse",cmd.stepAvgDiffuse);
        fprintf(fdout,FMTI,"nBuffDiffuse",cmd.nBuffDiffuse);
        fprintf(fdout,FMTI,"nValDiffuse",cmd.nValDiffuse);

        fprintf(fdout,FMTT,"computeVelAcf",
				cmd.computeVelAcf ? "true" : "false");

        fprintf(fdout,FMTT,"computeBulkViscosity",
				cmd.computeBulkViscosity ? "true" : "false");

        fprintf(fdout,FMTT,"computeTransport",
				cmd.computeTransport ? "true" : "false");
        fprintf(fdout,FMTI,"stepAcf",cmd.stepAcf);
        fprintf(fdout,FMTI,"stepAvgAcf",cmd.stepAvgAcf);
        fprintf(fdout,FMTI,"nBuffAcf",cmd.nBuffAcf);
        fprintf(fdout,FMTI,"nValAcf",cmd.nValAcf);

        fprintf(fdout,FMTT,"computeSTCorr",
				cmd.computeSTCorr ? "true" : "false");
        fprintf(fdout,FMTI,"stepCorr",cmd.stepCorr);
        fprintf(fdout,FMTI,"stepAvgCorr",cmd.stepAvgCorr);
        fprintf(fdout,FMTI,"nBuffCorr",cmd.nBuffCorr);
        fprintf(fdout,FMTI,"nValCorr",cmd.nValCorr);
        fprintf(fdout,FMTI,"nFunCorr",cmd.nFunCorr);

        fprintf(fdout,FMTR,"lattCorr_kx",cmd.lattCorr_kx);
        fprintf(fdout,FMTR,"lattCorr_ky",cmd.lattCorr_ky);
#ifdef THREEDIM
        fprintf(fdout,FMTR,"lattCorr_kz",cmd.lattCorr_kz);
#endif

        fprintf(fdout,FMTI,"icModel",cmd.icModel);
        fprintf(fdout,FMTI,"intMethod",cmd.intMethod);
        fprintf(fdout,FMTT,"dtime",cmd.dtimestr);
        fprintf(fdout,FMTR,"tstop",cmd.tstop);
        fprintf(fdout,FMTT,"dtout",cmd.dtoutstr);
        fprintf(fdout,FMTT,"dtoutinfo",cmd.dtoutinfostr);
        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,FMTT,"unitCells",cmd.unitCells);
        fprintf(fdout,FMTT,"nbodyprop",cmd.nbodyprop);
        fprintf(fdout,FMTT,"massprop",cmd.massprop);
#if (NDIM==3)
        fprintf(fdout,FMTT,"LxLyprop",cmd.LxLyprop);
#else
        fprintf(fdout,FMTT,"Lx",cmd.Lx);
#endif
        fprintf(fdout,FMTR,"eps11",cmd.eps11);
        fprintf(fdout,FMTR,"eps12",cmd.eps12);
        fprintf(fdout,FMTR,"eps22",cmd.eps22);
        fprintf(fdout,FMTR,"sigma11",cmd.sigma11);
        fprintf(fdout,FMTR,"sigma12",cmd.sigma12);
        fprintf(fdout,FMTR,"sigma22",cmd.sigma22);
        fprintf(fdout,FMTR,"Rcut11",cmd.Rcut11);
        fprintf(fdout,FMTR,"Rcut12",cmd.Rcut12);
        fprintf(fdout,FMTR,"Rcut22",cmd.Rcut22);
        fprintf(fdout,FMTI,"seed",cmd.seed);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR
														// CHECK 2D --- OK!!!
local void forcecalc_method_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"barnes") == 0)	*method_int = 0;
    if (strnull(method_str))				*method_int = 1;
    if (strcmp(method_str,"normal") == 0)	*method_int = 2;
    if (strcmp(method_str,"nblist") == 0)	*method_int = 3;
    if (strcmp(method_str,"direct") == 0)	*method_int = 4;
    if (strcmp(method_str,"cells") == 0)	*method_int = 5;
    if (strcmp(method_str,"direct2") == 0)	*method_int = 6;
    if (strcmp(method_str,"normal2") == 0)	*method_int = 7;
    if (strcmp(method_str,"barnes2") == 0)	*method_int = 8;
    if (strcmp(method_str,"normal3") == 0)	*method_int = 9;
    if (strcmp(method_str,"cells3") == 0)	*method_int = 10;
}

/*
local void MinimumDistanceBBodies(real lmax)
{
	bodyptr p, q;
	int k;
	real distMin, drpq, drpq2;
	vector dr;

	fprintf(gd.outlog,"\nFinding minimum distances between bodies ... ");
	distMin = lmax;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		DO_BODY(q, bodytab, p) {
			DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
			DO_COORD(k)
				dr[k]=dr[k]-((real)(nint(dr[k]/gdforce.Box[k])))*gdforce.Box[k];
			DOTVP(drpq2, dr, dr);
			drpq = rsqrt(drpq2);
			distMin = MIN(distMin, drpq);
		}
	}
	fprintf(gd.outlog,"distMin=%g\n",distMin);
	fflush(gd.outlog);
}
*/

local void SetRcut(void)								// CHECK 2D --- OK!!!
{
	bodyptr p;

	DO_BODY(p, bodytab, bodytab+gd.nbody)
		if (Type(p)==BODY1)
			Rcut(p)=gdforce.Rcut11Max;
		else
			Rcut(p)=gdforce.Rcut22Max;
		
}

local void SetBodyChemPot(bodyptr *btab)				// CHECK 2D --- OK!!!
{
	bodyptr btabtmp, p, q;

    btabtmp = (bodyptr) allocate((gd.nbody) * sizeof(body));

	p=btabtmp;
	DO_BODY(q, *btab, *btab+gd.nbody) {
		Mass(p) = Mass(q);
		Id(p) = Id(q);
		Type(p) = Type(q);
		SETV(Pos(p), Pos(q));
		SETV(Vel(p), Vel(q));
		Rcut(p) = Rcut(q);
		++p;
	}
	free(*btab);
	
    *btab = (bodyptr) allocate((gd.nbody+1) * sizeof(body));

	p=*btab;
	DO_BODY(q, btabtmp, btabtmp+gd.nbody) {
		Mass(p) = Mass(q);
		Id(p) = Id(q);
		Type(p) = Type(q);
		SETV(Pos(p), Pos(q));
		SETV(Vel(p), Vel(q));
		Rcut(p) = Rcut(q);
		++p;
	}
	free(btabtmp);

	Mass(p) = Mass(p-1);
	Id(p) = p-*btab+1;
	Type(p) = TESTBODYMU;
	CLRV(Pos(p));
	CLRV(Vel(p));
	Update(p) = FALSE;
	Rcut(p) = gd.RcutAllMax;
	fprintf(gd.outlog,"\n\nAdditional body for testing chem pot:\n");
	fprintf(gd.outlog,"\nId, Type, Mass: %d %d %g", Id(p), Type(p), Mass(p));
#if (NDIM==3)
	fprintf(gd.outlog,"\nPos: %g %g %g", Pos(p)[0], Pos(p)[1], Pos(p)[2]);
	fprintf(gd.outlog,"\nVel: %g %g %g", Vel(p)[0], Vel(p)[1], Vel(p)[2]);
#else
	fprintf(gd.outlog,"\nPos: %g %g", Pos(p)[0], Pos(p)[1]);
	fprintf(gd.outlog,"\nVel: %g %g", Vel(p)[0], Vel(p)[1]);
#endif
	fprintf(gd.outlog,"\nUpdate, Rcut: %d %g\n", Update(p), Rcut(p));
}

#define fpfnametmp		"fptabletmp.pot"

local void InputForcePotTable(void)						// CHECK 2D --- OK!!!
{
    char namebuf[256];
    struct stat buf;
    stream instr;
	stream outstr;
    pointForcePotptr p;
	char gato[1], firstline[20];
	real *row;
	int ncol;
	real ri, rf, dr;

    sprintf(namebuf, cmd.fnamePot);
    if (stat(namebuf, &buf) != 0) {
		error("\n\nInputForcePotTable : file %s does not exist!\n\n",cmd.fnamePot);
    } else {
        instr = stropen(namebuf, "r");
	}

	fprintf(stdout,"\n\nReading pot & force from file %s...\n",cmd.fnamePot);	
	fgets(firstline,200,instr);
	fgets(firstline,200,instr);
	fscanf(instr,"%1s",gato);
    in_int(instr, &nforcepot);
    if (nforcepot < 1)
        error("\n\nInputForcePotTable: nforcepot = %d is absurd\n\n", nforcepot);

	in_real(instr, &ri);
	in_real(instr, &rf);
	in_real(instr, &dr);

	fscanf(instr,"%1s",gato);

	in_real(instr, &cmd.sigma11);
	in_real(instr, &cmd.sigma12);
	in_real(instr, &cmd.sigma22);
	in_real(instr, &cmd.eps11);
	in_real(instr, &cmd.eps12);
	in_real(instr, &cmd.eps22);

	in_real(instr, &gdforce.fphi11);
	in_real(instr, &gdforce.fphi12);
	in_real(instr, &gdforce.fphi22);
	in_real(instr, &gdforce.fa11);
	in_real(instr, &gdforce.fa12);
	in_real(instr, &gdforce.fa22);

	gdforce.ssq11 = cmd.sigma11*cmd.sigma11;
	gdforce.ssq12 = cmd.sigma12*cmd.sigma12;
	gdforce.ssq22 = cmd.sigma22*cmd.sigma22;

    forcepottab = (pointForcePotptr) allocate(nforcepot * sizeof(pointForcePot));

	fprintf(stdout,"nforcepot : %d\n", nforcepot);

	ncol=8;
	row = (realptr) allocate(ncol*sizeof(real));
	for (p=forcepottab; p<forcepottab+nforcepot; p++) {
		in_vector_ndim(instr, row, ncol);
		idPos(p) = row[0];
		rPos(p) = row[1];
		Pot11(p) = row[2];
		Force11(p) = row[3];
		Pot12(p) = row[4];
		Force12(p) = row[5];
		Pot22(p) = row[6];
		Force22(p) = row[7];
	}

	outstr = stropen(fpfnametmp,"w!");
	fprintf(outstr,"# nforcepot ri rf dr\n");
	fprintf(outstr,"# sigma11 sigma12 sigma22 eps11 eps12 eps22 %s\n",
					"fphi11 fphi12 fphi22 fa11 fa12 fa22");
	fprintf(outstr,"# %d %g %g %g\n", nforcepot, ri, rf, dr);
	fprintf(outstr,"# %g %g %g %g %g %g %g %g %g %g %g %g\n",
		cmd.sigma11,cmd.sigma12,cmd.sigma22,
		cmd.eps11,cmd.eps12,cmd.eps22,
		gdforce.fphi11,gdforce.fphi12,gdforce.fphi22,
		gdforce.fa11,gdforce.fa12,gdforce.fa22);
	for (p=forcepottab; p<forcepottab+nforcepot; p++) {
		fprintf(outstr,"%d %g %g %g %g %g %g %g\n",
		idPos(p),rPos(p),Pot11(p),Force11(p),Pot12(p),Force12(p),Pot22(p),Force22(p));
	}
	fclose(outstr);
}

#undef fpfnametmp

