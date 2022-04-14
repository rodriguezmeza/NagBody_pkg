/*==============================================================================
	MODULE: startrun.c			[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Initialize datanaly_md
	Language: C
	Use: 'startrun();'
	Routines and functions: testdata, ReadParameterFile, PrintParameterFile
	Modules, routines and external headers:
	Coments and notes:
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: January 2007; November 2008;
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
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

#define logfile			"analysis_md.log"


void startrun(string head0, string head1, string head2, string head3)
{
//	RestartFlag=0;					// De gadget ... se debe inicializar en 0...

	gd.x_autoscale = TRUE;
	gd.y_autoscale = TRUE;
#ifdef THREEDIM
	gd.z_autoscale = TRUE;
#endif

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

//    LicDriver("analysis_md");

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

#define parameter_null	"parameters_null"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
    cmd.in = GetParam("in");
    cmd.infmt = GetParam("infmt");
    cmd.out = GetParam("out");
    cmd.outfmt = GetParam("outfmt");
	cmd.basedir = GetParam("basedir");
	cmd.isnap = GetiParam("isnap");
	cmd.fsnap = GetiParam("fsnap");

//	cmd.temperature = GetdParam("temperature");
//	cmd.density = GetdParam("density");

//	cmd.eps11 = GetdParam("eps11");
//	cmd.eps12 = GetdParam("eps12");
//	cmd.eps22 = GetdParam("eps22");
//	cmd.sigma11 = GetdParam("sigma11");
//	cmd.sigma12 = GetdParam("sigma12");
//	cmd.sigma22 = GetdParam("sigma22");
//	cmd.Rcut11 = GetdParam("Rcut11");
//	cmd.Rcut12 = GetdParam("Rcut12");
//	cmd.Rcut22 = GetdParam("Rcut22");

//	cmd.nbodyprop = GetParam("nbodyprop");
//	cmd.massprop = GetParam("massprop");

//#ifdef THREEDIM
//		cmd.LxLyprop = GetParam("LxLyprop");
//#else
//		cmd.Lx = GetParam("Lx");
//#endif

	cmd.stepAvg = GetiParam("stepAvg");
	cmd.sizeHist = GetiParam("sizeHist");
	cmd.rangeVal = GetdParam("rangeVal");

//	cmd.stepAvgRhoAxes = GetiParam("stepAvgRhoAxes");
//	cmd.sizeHistRhoAxes = GetiParam("sizeHistRhoAxes");

//	cmd.stepAvgNFrecAxes = GetiParam("stepAvgNFrecAxes");
//	cmd.sizeHistNFrecAxes = GetiParam("sizeHistNFrecAxes");

//	cmd.stepAvgVel = GetiParam("stepAvgVel");
//	cmd.sizeHistVel = GetiParam("sizeHistVel");
//	cmd.rangeVel = GetdParam("rangeVel");

//	cmd.stepAvgRdf = GetiParam("stepAvgRdf");
//	cmd.sizeHistRdf = GetiParam("sizeHistRdf");
//	cmd.rangeRdf = GetdParam("rangeRdf");

	cmd.unitsset = GetParam("unitsset");
	cmd.options = GetParam("options");
	cmd.data_analysis_type = GetParam("analysis_type");
	cmd.reductionFac = GetdParam("reductionFac");
	cmd.bodiesID = GetParam("bodiesID");

	cmd.xrange = GetParam("xrange");
	cmd.yrange = GetParam("yrange");
#ifdef THREEDIM
	cmd.zrange = GetParam("zrange");
#endif

#ifdef THREEDIM
	cmd.xmin = GetdParam("xmin");
	cmd.xmax = GetdParam("xmax");
	cmd.ymin = GetdParam("ymin");
	cmd.ymax = GetdParam("ymax");
	cmd.zmin = GetdParam("zmin");
	cmd.zmax = GetdParam("zmax");
#endif

	cmd.usingcolumns = GetParam("usingcolumns");
	cmd.usingrows = GetParam("usingrows");

	cmd.xlabel = GetParam("xlabel");
	cmd.ylabel = GetParam("ylabel");
	cmd.plotlabel = GetParam("plotlabel");
	cmd.labelfontsize = GetdParam("labelfontsize");
	cmd.plotjoined = GetbParam("plotjoined");
	cmd.withdots = GetbParam("withdots");
	cmd.withsymbols = GetbParam("withsymbols");
	cmd.symbolsize = GetdParam("symbolsize");

	cmd.labelfontweight = GetiParam("labelfontweight");
	cmd.symbolweight	= GetiParam("symbolweight");
	cmd.nlsize			= GetdParam("nlsize");
	cmd.linewidth		= GetiParam("linewidth");
	cmd.axeswidth		= GetiParam("axeswidth");
	cmd.symbolcolor		= GetiParam("symbolcolor");

	cmd.pl_a = GetParam("aspectratio");
	cmd.pl_dev = GetParam("dev");
	cmd.pl_geo = GetParam("geometry");
	cmd.pl_ori = GetParam("ori");
	cmd.pl_bg = GetParam("background");
	cmd.pl_ncol0 = GetParam("ncol0");
	cmd.pl_ncol1 = GetParam("ncol1");
}

local void startrun_Common(void)
{
	char *pch;
	char bodiesIDtmp[256];					// 8 x MAXBODIESID < 256
	real dt1, dt2;
	int xrangeflag=1, yrangeflag=1;
#ifdef THREEDIM
	int zrangeflag=1;
#endif

    if(!(gd.outlog=fopen(logfile,"w")))
        error("error opening file '%s' \n",logfile);

	if (cmd.isnap<0 || cmd.fsnap<0)
		error("\n\nstartrun: isnap and fsnap must be >= 0\n");
	if (cmd.isnap > cmd.fsnap)
		error("\n\nstartrun: isnap must be < fsnap\n");

	gd.headerfmt = "snap-tlj-ascii";

//	if (!(sscanf(cmd.nbodyprop, "%ld/%ld", &gd.nbody1, &gd.nbody2) == 2))
//		error("\nstartrun : nbodyprop must be in the form of nbody1/nbody2\n\n");
//	if (gd.nbody1==0 && gd.nbody2==0)
//		error("\nstartrun : nbody1 and nbody2 must not be zero at the same time\n\n");
//	if (gd.nbody1 < 0)
//		error("startrun: absurd value for nbody1\n");
//	if (gd.nbody2 < 0)
//		error("startrun: absurd value for nbody2\n");

//	if (!(sscanf(cmd.massprop, "%lf/%lf", &gd.mass1, &gd.mass2) == 2))
//		error("\nstartrun : massprop must be in the form of mass1/mass2\n\n");
//	if (gd.mass1==0 && gd.mass2==0)
//		error("\nstartrun : mass1 and mass2 must not be zero at the same time\n\n");
//	if (gd.mass1 < 0)
//		error("startrun: absurd value for mass1\n");
//	if (gd.mass2 < 0)
//		error("startrun: absurd value for mass2\n");

//#ifdef THREEDIM
//		if (!(sscanf(cmd.LxLyprop, "%lf/%lf", &gd.Lx, &gd.Ly) == 2))
//			error("\nstartrun : LxLyprop must be in the form of Lx/Ly\n\n");
//		if (gd.Lx <= 0)
//			error("startrun: absurd value for Lx\n");
//		if (gd.Ly <= 0)
//			error("startrun: absurd value for Ly\n");
//#else
//		if (!(sscanf(cmd.Lx, "%lf", &gd.Lx) == 1))
//			error("\nstartrun : Lx must be in the form a real number\n\n");
//		if (gd.Lx <= 0)
//			error("startrun: absurd value for Lx\n");
//#endif

	if (!(sscanf(cmd.unitsset, "%lf:%lf:%lf", 
		&gd.unitLength, &gd.unitMass, &gd.unitEnergy) == 3))
		error("\nstartrun_Common : unitsset must be in the form of Length:Mass:Energy\n\n");

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
	} else
		if (scanopt(cmd.data_analysis_type, "bodies-anim") 
			|| scanopt(cmd.data_analysis_type, "bodies-animtrajectory"))
			error("\n\nstartrun_Common: analysis_type=bodies-anim or bodies-animtrajectory %s\n",
				"must be run with bodiesID not null\n");

	if (strcmp(cmd.xrange,"autoscale"))
		gd.x_autoscale=FALSE;
	if (!(sscanf(cmd.xrange, "%lf:%lf", &gd.xmin, &gd.xmax) == 2))
		xrangeflag=0;
	if ( !gd.x_autoscale && xrangeflag==0) 
		error("\nGet_Parameters: xrange must be autoscale or in the form xmin:xmax\n\n");
	if (xrangeflag != 0 && gd.xmin >= gd.xmax)
		error("\nGet_Parameters: xmin must be < xmax : (xmin,xmax)=(%lf,%lf)\n\n",gd.xmin,gd.xmax);

	if (strcmp(cmd.yrange,"autoscale"))
		gd.y_autoscale=FALSE;
	if (!(sscanf(GetParam("yrange"), "%lf:%lf", &gd.ymin, &gd.ymax) == 2))
		yrangeflag=0;
	if ( !gd.y_autoscale && yrangeflag==0) 
		error("\nGet_Parameters: yrange must be autoscale or in the form ymin:ymax\n\n");
	if (yrangeflag != 0 && gd.ymin >= gd.ymax)
		error("\nGet_Parameters: ymin must be < ymax : (ymin,ymax)=(%lf,%lf)\n\n",gd.ymin,gd.ymax);

#ifdef THREEDIM
	if (strcmp(cmd.zrange,"autoscale"))
		gd.z_autoscale=FALSE;
	if (!(sscanf(GetParam("zrange"), "%lf:%lf", &gd.zmin, &gd.zmax) == 2))
		zrangeflag=0;
	if ( !gd.z_autoscale && zrangeflag==0)
		error("\nGet_Parameters: zrange must be autoscale or in the form zmin:zmax\n\n");
	if (zrangeflag != 0 && gd.zmin >= gd.zmax)
		error("\nGet_Parameters: zmin must be < zmax : (zmin,zmax)=(%lf,%lf)\n\n",gd.zmin,gd.zmax);
#endif

	if (scanopt(cmd.data_analysis_type,"thermo-avg") && scanopt(cmd.options, "errorbar")) {
		if (!(sscanf(cmd.usingcolumns, "%d:%d:%d", 
			&gd.column1, &gd.column2, &gd.column3) == 3))
			error("\nstartrun_cmdline: usingcolumns must be in the form c1:c2:c3\n\n");
	} else
		if (!(sscanf(cmd.usingcolumns, "%d:%d", &gd.column1, &gd.column2) == 2))
			error("\nstartrun_cmdline: usingcolumns must be in the form c1:c2\n\n");

	if (!(sscanf(cmd.usingrows, "%d:%d", &gd.row1, &gd.row2) == 2))
		error("\nGet_Parameters: usingrows must be in the form r1:r2\n\n");

	CheckParameters();
}

local void startrun_ParamStat(void)						// CHECK 2D --- OK!!!
{
	if (GetParamStat("stepAvg") & ARGPARAM)
		cmd.stepAvg = GetiParam("stepAvg");         
	if (GetParamStat("sizeHist") & ARGPARAM)
		cmd.sizeHist = GetiParam("sizeHist");
	if (GetParamStat("rangeVal") & ARGPARAM)
		cmd.rangeVal = GetdParam("rangeVal");

//	if (GetParamStat("stepAvgVel") & ARGPARAM)
//		cmd.stepAvgVel = GetiParam("stepAvgVel");         
//	if (GetParamStat("sizeHistVel") & ARGPARAM)
//		cmd.sizeHistVel = GetiParam("sizeHistVel");
//	if (GetParamStat("rangeVel") & ARGPARAM)
//		cmd.rangeVel = GetdParam("rangeVel");
//	if (GetParamStat("stepAvgRdf") & ARGPARAM)
//		cmd.stepAvgRdf = GetiParam("stepAvgRdf");         
//	if (GetParamStat("sizeHistRdf") & ARGPARAM)
//		cmd.sizeHistRdf = GetiParam("sizeHistRdf");
//	if (GetParamStat("rangeRdf") & ARGPARAM)
//		cmd.rangeRdf = GetdParam("rangeRdf");

	if (GetParamStat("reductionFac") & ARGPARAM)
		cmd.reductionFac = GetdParam("reductionFac");

	if (GetParamStat("options") & ARGPARAM)
		cmd.options = GetParam("options");
}

local void CheckParameters(void)
{
	if ( cmd.stepAvg == 0 )
		error("\n\nCheckParameters: stepAvg must be != 0\n");
//	if ( cmd.stepAvgRhoAxes == 0 )
//		error("\n\nCheckParameters: stepAvgRhoAxes must be != 0\n");
//	if ( cmd.stepAvgNFrecAxes == 0 )
//		error("\n\nCheckParameters: stepAvgNFrecAxes must be != 0\n");
//	if ( cmd.stepAvgVel == 0 )
//		error("\n\nCheckParameters: stepAvgVel must be != 0\n");
//	if ( cmd.stepAvgRdf == 0 )
//		error("\n\nCheckParameters: stepAvgRdf must be != 0\n");

	if (gd.unitLength<=0 || gd.unitMass<=0 || gd.unitEnergy<=0)
		error("\nCheckParameters : Length, Mass or Energy scales must be >0\n\n");

	if (strnull(cmd.data_analysis_type))
		error("\nrun from command line as:\n%s",
			"\t'analysis_md analysis_type=analysis_type_name'\n\n");

	if ( cmd.reductionFac <= 0.0 || cmd.reductionFac > 1.0 )
		error("\n\nCheckParameters: reductionFac must be in the inverval (0,1]\n");
}


#undef parameter_null
#undef logfile

void Header_to_Global(void)
{
//    if ( (strcmp(cmd.infmt,"snap-blj-ascii") == 0)	||
//		(strcmp(cmd.infmt,"snap-blj-pv") == 0)		||
//		(strcmp(cmd.infmt,"snap-blj-bin") == 0) ) {

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
//	} else {
//#if (NDIM==3)
//		if (cmd.density!=0.0)
//			gd.Lz = ( (real) gd.nbody1 * gd.mass1 + 
//				(real) gd.nbody2 * gd.mass2 )/(gd.Lx*gd.Ly*cmd.density);
//		else
//			error("\nyou should give the density value\n");
//#else
//		if (cmd.density!=0.0)
//			gd.Ly = ( (real) gd.nbody1 * gd.mass1 + 
//				(real) gd.nbody2 * gd.mass2 )/(gd.Lx*cmd.density);
//		else
//			error("\nyou should give the density value\n");
//#endif
//	}

//	gd.nbody=gd.nbody1+gd.nbody2;
	if (gd.nbody != gd.nbody1+gd.nbody2)
		error("\nnbody not equal to nbody1+nbody2\n");

	gd.Box[0] = gd.Lx; gd.Box[1] = gd.Ly;
#ifdef THREEDIM
	gd.Box[2] = gd.Lz;
#endif
}

void Global_to_Header(void)
{
/*
#if (NDIM==3)
		if (cmd.density!=0.0)
			gd.Lz = ( (real) gd.nbody1 * gd.mass1 + 
				(real) gd.nbody2 * gd.mass2 )/(gd.Lx*gd.Ly*cmd.density);
		else
			error("\nyou should give the density value\n");
#else
		if (cmd.density!=0.0)
			gd.Ly = ( (real) gd.nbody1 * gd.mass1 + 
				(real) gd.nbody2 * gd.mass2 )/(gd.Lx*cmd.density);
		else
			error("\nyou should give the density value\n");
#endif

	gd.Box[0] = gd.Lx; gd.Box[1] = gd.Ly;
#ifdef THREEDIM
	gd.Box[2] = gd.Lz;
#endif
*/
//    if ( (strcmp(cmd.outfmt,"snap-blj-ascii") == 0)	||
//		(strcmp(cmd.outfmt,"snap-blj-pv") == 0)		||
//		(strcmp(cmd.outfmt,"snap-blj-bin") == 0) ) {

		hdr.nbody		= gd.nbody;
		hdr.nbody1		= gd.nbody1;
		hdr.nbody2		= gd.nbody2;
		hdr.tnow		= gd.tnow;
		hdr.temperature	= cmd.temperature;
		hdr.density		= cmd.density;
		hdr.mass1		= gd.mass1;
		hdr.mass2		= gd.mass2;
		hdr.Lx			= gd.Lx;
		hdr.Ly			= gd.Ly;
#ifdef THREEDIM
		hdr.Lz			= gd.Lz;
#endif
		hdr.eps11		= cmd.eps11;
		hdr.eps12		= cmd.eps12;
		hdr.eps22		= cmd.eps22;
		hdr.sigma11		= cmd.sigma11;
		hdr.sigma12		= cmd.sigma12;
		hdr.sigma22		= cmd.sigma22;
		hdr.Rcut11		= cmd.Rcut11;
		hdr.Rcut12		= cmd.Rcut12;
		hdr.Rcut22		= cmd.Rcut22;
//	}
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

	IPName(cmd.stepAvg,"stepAvg");
	IPName(cmd.sizeHist,"sizeHist");
	RPName(cmd.rangeVal,"rangeVal");

//	IPName(cmd.stepAvgRhoAxes,"stepAvgRhoAxes");
//	IPName(cmd.sizeHistRhoAxes,"sizeHistRhoAxes");
//	IPName(cmd.stepAvgNFrecAxes,"stepAvgNFrecAxes");
//	IPName(cmd.sizeHistNFrecAxes,"sizeHistNFrecAxes");

//	IPName(cmd.stepAvgVel,"stepAvgVel");
//	IPName(cmd.sizeHistVel,"sizeHistVel");
//	RPName(cmd.rangeVel,"rangeVel");

//	IPName(cmd.stepAvgRdf,"stepAvgRdf");
//	IPName(cmd.sizeHistRdf,"sizeHistRdf");
//	RPName(cmd.rangeRdf,"rangeRdf");

	SPName(cmd.unitsset,"unitsset",100);

	SPName(cmd.options,"options",100);
//	RPName(dtout,"dtout");
//	IPName(nbody,"nbody");
	SPName(cmd.data_analysis_type,"analysis_type",100);

	RPName(cmd.reductionFac,"reductionFac");

	SPName(cmd.bodiesID,"bodiesID",256);
	SPName(cmd.xrange,"xrange",100);
	SPName(cmd.yrange,"yrange",100);
#ifdef THREEDIM
	SPName(cmd.zrange,"zrange",100);
#endif

#ifdef THREEDIM
	RPName(cmd.xmin,"xmin");
	RPName(cmd.xmax,"xmax");
	RPName(cmd.ymin,"ymin");
	RPName(cmd.ymax,"ymax");
	RPName(cmd.zmin,"zmin");
	RPName(cmd.zmax,"zmax");
#endif

	SPName(cmd.usingcolumns,"usingcolumns",100);
	SPName(cmd.usingrows,"usingrows",100);
	SPName(cmd.xlabel,"xlabel",100);
	SPName(cmd.ylabel,"ylabel",100);
	RPName(cmd.labelfontsize,"labelfontsize");
	SPName(cmd.plotlabel,"plotlabel",100);
	BPName(cmd.plotjoined,"plotjoined");
	BPName(cmd.withsymbols,"withsymbols");
	BPName(cmd.withdots,"withdots");
	RPName(cmd.symbolsize,"symbolsize");

	IPName(cmd.labelfontweight,"labelfontweight");
	RPName(cmd.nlsize,"nlsize");
	IPName(cmd.linewidth,"linewidth");
	IPName(cmd.axeswidth,"axeswidth");
	IPName(cmd.symbolweight,"symbolweight");
	IPName(cmd.symbolcolor,"symbolcolor");

	SPName(cmd.pl_a,"aspectratio",1);
	SPName(cmd.pl_dev,"dev",1);
	SPName(cmd.pl_geo,"geometry",50);
	SPName(cmd.pl_ori,"ori",1);
	SPName(cmd.pl_bg,"background",50);
	SPName(cmd.pl_ncol0,"ncol0",1);
	SPName(cmd.pl_ncol1,"ncol1",1);

//	RPName(cmd.temperature,"temperature");
//	SPName(cmd.nbodyprop,"nbodyprop",100);
//	SPName(cmd.massprop,"massprop",100);
//#ifdef THREEDIM
//	SPName(cmd.LxLyprop,"LxLyprop",100);
//#else
//	SPName(cmd.Lx,"Lx",100);
//#endif

//	RPName(cmd.eps11,"eps11");
//	RPName(cmd.eps12,"eps12");
//	RPName(cmd.eps22,"eps22");
//	RPName(cmd.sigma11,"sigma11");
//	RPName(cmd.sigma12,"sigma12");
//	RPName(cmd.sigma22,"sigma22");
//	RPName(cmd.Rcut11,"Rcut11");
//	RPName(cmd.Rcut12,"Rcut12");
//	RPName(cmd.Rcut22,"Rcut22");

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

    sprintf(buf,"%s%s",fname,"-analysis_md");
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
        fprintf(fdout,FMTT,"basedir",cmd.basedir);

        fprintf(fdout,FMTI,"isnap",cmd.isnap);
        fprintf(fdout,FMTI,"fsnap",cmd.fsnap);

//        fprintf(fdout,FMTR,"temperature",cmd.temperature);
//        fprintf(fdout,FMTR,"density",cmd.density);

//        fprintf(fdout,FMTI,"stepEquil",stepEquil);

        fprintf(fdout,FMTI,"stepAvg",cmd.stepAvg);
        fprintf(fdout,FMTI,"sizeHist",cmd.sizeHist);
        fprintf(fdout,FMTR,"rangeVal",cmd.rangeVal);

//        fprintf(fdout,FMTI,"stepAvgRhoAxes",cmd.stepAvgRhoAxes);
//        fprintf(fdout,FMTI,"sizeHistRhoAxes",cmd.sizeHistRhoAxes);

//        fprintf(fdout,FMTI,"stepAvgNFrecAxes",cmd.stepAvgNFrecAxes);
//        fprintf(fdout,FMTI,"sizeHistNFrecAxes",cmd.sizeHistNFrecAxes);

//        fprintf(fdout,FMTI,"stepAvgVel",cmd.stepAvgVel);
//        fprintf(fdout,FMTI,"sizeHistVel",cmd.sizeHistVel);
//        fprintf(fdout,FMTR,"rangeVel",cmd.rangeVel);

//        fprintf(fdout,FMTI,"stepAvgRdf",cmd.stepAvgRdf);
//        fprintf(fdout,FMTI,"sizeHistRdf",cmd.sizeHistRdf);
//        fprintf(fdout,FMTR,"rangeRdf",cmd.rangeRdf);

//        fprintf(fdout,FMTR,"dtime",dtime);

//        fprintf(fdout,FMTR,"dtout",dtout);

        fprintf(fdout,FMTT,"unitsset",cmd.unitsset);

        fprintf(fdout,FMTT,"options",cmd.options);
//        fprintf(fdout,FMTI,"nbody",nbody);
        fprintf(fdout,FMTT,"analysis_type",cmd.data_analysis_type);

//        fprintf(fdout,FMTT,"nbodyprop",cmd.nbodyprop);
//        fprintf(fdout,FMTT,"massprop",cmd.massprop);
//#ifdef THREEDIM
//        fprintf(fdout,FMTT,"LxLyprop",cmd.LxLyprop);
//#else
//        fprintf(fdout,FMTT,"Lx",cmd.Lx);
//#endif
//        fprintf(fdout,FMTR,"eps11",cmd.eps11);
//        fprintf(fdout,FMTR,"eps12",cmd.eps12);
//        fprintf(fdout,FMTR,"eps22",cmd.eps22);
//        fprintf(fdout,FMTR,"sigma11",cmd.sigma11);
//        fprintf(fdout,FMTR,"sigma12",cmd.sigma12);
//        fprintf(fdout,FMTR,"sigma22",cmd.sigma22);
//        fprintf(fdout,FMTR,"Rcut11",cmd.Rcut11);
//        fprintf(fdout,FMTR,"Rcut12",cmd.Rcut12);
//        fprintf(fdout,FMTR,"Rcut22",cmd.Rcut22);

        fprintf(fdout,FMTR,"reductionFac",cmd.reductionFac);

        fprintf(fdout,FMTT,"bodiesID",cmd.bodiesID);

        fprintf(fdout,FMTT,"xrange",cmd.xrange);
        fprintf(fdout,FMTT,"yrange",cmd.yrange);
#ifdef THREEDIM
        fprintf(fdout,FMTT,"zrange",cmd.zrange);
#endif

#ifdef THREEDIM
        fprintf(fdout,FMTR,"xmin",cmd.xmin);
        fprintf(fdout,FMTR,"xmax",cmd.xmax);
        fprintf(fdout,FMTR,"ymin",cmd.ymin);
        fprintf(fdout,FMTR,"ymax",cmd.ymax);
        fprintf(fdout,FMTR,"zmin",cmd.zmin);
        fprintf(fdout,FMTR,"zmax",cmd.zmax);
#endif

        fprintf(fdout,FMTT,"usingcolumns",cmd.usingcolumns);
        fprintf(fdout,FMTT,"usingrows",cmd.usingrows);

        fprintf(fdout,FMTT,"xlabel",cmd.xlabel);
        fprintf(fdout,FMTT,"ylabel",cmd.ylabel);
        fprintf(fdout,FMTT,"plotlabel",cmd.plotlabel);
        fprintf(fdout,FMTR,"labelfontsize",cmd.labelfontsize);
        fprintf(fdout,FMTT,"plotjoined",cmd.plotjoined ? "true" : "false");
        fprintf(fdout,FMTT,"withdots",cmd.withdots ? "true" : "false");
        fprintf(fdout,FMTT,"withsymbols",cmd.withsymbols ? "true" : "false");
        fprintf(fdout,FMTR,"symbolsize",cmd.symbolsize);

        fprintf(fdout,FMTI,"labelfontweight",cmd.labelfontweight);
        fprintf(fdout,FMTR,"nlsize",cmd.nlsize);
        fprintf(fdout,FMTI,"linewidth",cmd.linewidth);
        fprintf(fdout,FMTI,"axeswidth",cmd.axeswidth);
        fprintf(fdout,FMTI,"symbolweight",cmd.symbolweight);
        fprintf(fdout,FMTI,"symbolcolor",cmd.symbolcolor);

        fprintf(fdout,FMTT,"aspectratio",cmd.pl_a);
        fprintf(fdout,FMTT,"dev",cmd.pl_dev);
        fprintf(fdout,FMTT,"geometry",cmd.pl_geo);
        fprintf(fdout,FMTT,"ori",cmd.pl_ori);
        fprintf(fdout,FMTT,"background",cmd.pl_bg);
        fprintf(fdout,FMTT,"ncol0",cmd.pl_ncol0);
        fprintf(fdout,FMTT,"ncol1",cmd.pl_ncol1);

        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR
