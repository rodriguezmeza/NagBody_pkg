/*==============================================================================
	MODULE: start_run.c				[nplot2d]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: routines to initialize the main code
	Language: C
	Use: 'StartRun();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#include "globaldefs.h"

local void StartRunParameterfile(void);
local void StartRunCmdline(void);
local void ReadParameterFile(char *);
local void PrintParameterFile(char *);
local void PrintParametersLog(void);

local void ReadParametersCmdline(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

local void scanbOption(string, bool *, int *, int, int, string);
local void scaniOption(string, int *, int *, int, int, string);

#define logfile			"nplot2d.log"

void StartRun(string head0, string head1, string head2, string head3)
{
    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

//    LicDriver("nplot2d");

	gd.x_autoscale = TRUE;
	gd.y_autoscale = TRUE;

    if(!(gd.outlog=fopen(logfile,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",logfile);
        exit(0);
    }

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		StartRunParameterfile();
	else
		StartRunCmdline();
}

local void StartRunParameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
	startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-nplot2d"

local void StartRunCmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
    cmd.inputfile			= GetParam("inputfile");
    cmd.legends				= GetParam("legends");
	cmd.legendspos			= GetiParam("legendspos");
    cmd.witherrorbars		= GetParam("witherrorbars");
	cmd.errorbarstype		= GetiParam("errorbarstype");
//    cmd.outputformat		= GetParam("outputformat");

	cmd.xrange				= GetParam("xrange");
	cmd.yrange				= GetParam("yrange");
	cmd.graphicsarray		= GetParam("graphicsarray");
	cmd.axesorigin			= GetParam("axesorigin");
	cmd.usingcolumns		= GetParam("usingcolumns");
	cmd.plottype			= GetiParam("plottype");

	cmd.axestype			= GetiParam("axestype");
	cmd.axeswidth			= GetiParam("axeswidth");
	cmd.axescolor			= GetiParam("axescolor");
	cmd.xlabel				= GetParam("xlabel");
	cmd.ylabel				= GetParam("ylabel");
	cmd.labelcolor			= GetiParam("labelcolor");

	cmd.nlxdigmax			= GetiParam("nlxdigmax");
	cmd.nlydigmax			= GetiParam("nlydigmax");
	cmd.nlsize				= GetdParam("nlsize");
	cmd.nlcolor				= GetiParam("nlcolor");

	cmd.legendsfontsize		= GetdParam("legendsfontsize");
	cmd.legendsfontweight	= GetiParam("legendsfontweight");

	cmd.fontset				= GetiParam("fontset");
	cmd.fontchr				= GetiParam("fontchr");
	cmd.labelfontsize		= GetdParam("labelfontsize");
	cmd.labelfontweight		= GetiParam("labelfontweight");
	cmd.plotlabel			= GetParam("plotlabel");
	cmd.plotjoined			= GetParam("plotjoined");
	cmd.linetype			= GetParam("linetype");
	cmd.linewidth			= GetiParam("linewidth");
	cmd.linecolor			= GetiParam("linecolor");
	cmd.withdots			= GetParam("withdots");
	cmd.withsymbols			= GetParam("withsymbols");
	cmd.symboltype			= GetParam("symboltype");
	cmd.symbolcolor			= GetiParam("symbolcolor");
	cmd.symbolweight		= GetiParam("symbolweight");
	cmd.symbolsize			= GetdParam("symbolsize");

	cmd.text1				= GetParam("text1");
	cmd.text1side			= GetParam("text1side");
	cmd.text1disp			= GetdParam("text1disp");
	cmd.text1pos			= GetdParam("text1pos");
	cmd.text1just			= GetdParam("text1just");
	cmd.text1size			= GetdParam("text1size");
	cmd.text1weight			= GetiParam("text1weight");
	cmd.text1color			= GetiParam("text1color");

	cmd.text2				= GetParam("text2");
	cmd.text2side			= GetParam("text2side");
	cmd.text2disp			= GetdParam("text2disp");
	cmd.text2pos			= GetdParam("text2pos");
	cmd.text2just			= GetdParam("text2just");
	cmd.text2size			= GetdParam("text2size");
	cmd.text2weight			= GetiParam("text2weight");
	cmd.text2color			= GetiParam("text2color");

	cmd.text3				= GetParam("text3");
	cmd.text3side			= GetParam("text3side");
	cmd.text3disp			= GetdParam("text3disp");
	cmd.text3pos			= GetdParam("text3pos");
	cmd.text3just			= GetdParam("text3just");
	cmd.text3size			= GetdParam("text3size");
	cmd.text3weight			= GetiParam("text3weight");
	cmd.text3color			= GetiParam("text3color");

	cmd.text4				= GetParam("text4");
	cmd.text4x				= GetdParam("text4x");
	cmd.text4y				= GetdParam("text4y");
	cmd.text4dx				= GetdParam("text4dx");
	cmd.text4dy				= GetdParam("text4dy");
	cmd.text4just			= GetdParam("text4just");
	cmd.text4size			= GetdParam("text4size");
	cmd.text4weight			= GetiParam("text4weight");
	cmd.text4color			= GetiParam("text4color");

	cmd.locatemode			= GetbParam("locatemode");

	cmd.xaxis				= GetbParam("xaxis");
	cmd.yaxis				= GetbParam("yaxis");
	cmd.frame				= GetbParam("frame");

	cmd.framestyle			= GetParam("framestyle");

	cmd.gridlines			= GetbParam("gridlines");
	cmd.epilog				= GetParam("epilog");

//---------MENU COMMAND LINE PLPLOT OPTIONS------------------------------------
// Option set to false == use default, they are bool

	cmd.pl_showall			= GetbParam("showall");
	cmd.pl_h				= GetbParam("h");
	cmd.pl_v				= GetbParam("v");
	cmd.pl_verbose			= GetbParam("verbose");
	cmd.pl_debug			= GetbParam("debug");
	cmd.pl_hack				= GetbParam("hack");
	cmd.pl_dev				= GetParam("dev");
	cmd.pl_o				= GetParam("outputfile");
	cmd.pl_display			= GetParam("display");

	cmd.pl_geometry			= GetParam("geometry");
	cmd.pl_wplt				= GetParam("wplt");

	cmd.pl_mar				= GetParam("mar");

	cmd.pl_a				= GetParam("aspectratio");
	cmd.pl_jx				= GetParam("jx");
	cmd.pl_jy				= GetParam("jy");
	cmd.pl_ori				= GetParam("ori");
	cmd.pl_freeaspect		= GetbParam("freeaspect");
	cmd.pl_portrait			= GetbParam("portrait");
	cmd.pl_width			= GetParam("width");
	cmd.pl_bg				= GetParam("background");
	cmd.pl_ncol0			= GetParam("ncol0");
	cmd.pl_ncol1			= GetParam("ncol1");
	cmd.pl_fam				= GetbParam("fam");
	cmd.pl_fsiz				= GetParam("fsiz");
	cmd.pl_fbeg				= GetParam("fbeg");
	cmd.pl_finc				= GetParam("finc");
	cmd.pl_fflen			= GetParam("fflen");
	cmd.pl_nopixmap			= GetbParam("nopixmap");
	cmd.pl_db				= GetbParam("db");
	cmd.pl_np				= GetbParam("np");
	cmd.pl_bufmax			= GetParam("bufmax");
//	cmd.pl_server_name		= GetParam("server_name");
//	cmd.pl_plserver			= GetParam("plserver");
//	cmd.pl_plwindow			= GetParam("plwindow");
//	cmd.pl_tcl_cmd			= GetParam("tcl_cmd");
//	cmd.pl_auto_path		= GetParam("auto_path");
//	cmd.pl_tk_file			= GetParam("tk_file");
	cmd.pl_dpi				= GetParam("dpi");
	cmd.pl_compression		= GetParam("compression");
	cmd.pl_drvopt			= GetParam("drvopt");
//------------------------------------------------------------------------------
}

local void startrun_Common(void)
{
	char *pch;
	int i, npcols, nitems;
	short flag;
	char *pusingcolumns[30], *framestyleitems[30];
	char inputfiletmp[200], legendstmp[100],
		usingcolumnstmp[100], framestyletmp[50];
	int xrangeflag=1, yrangeflag=1;

	if (strnull(cmd.inputfile)) {
		fprintf(stdout,"\nstartrun_Common: no inputfile was given making data ...\n");
		gd.nfiles=1;							// To test data...
	} else {
		strcpy(inputfiletmp,cmd.inputfile);
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
		fprintf(stdout,"num. of files in inputfile %s =%d\n",cmd.inputfile,gd.nfiles);
	}

	if (!strnull(cmd.legends)) {
		strcpy(legendstmp,cmd.legends);
		fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",legendstmp);
		gd.nlegends=0;
		pch = strtok(legendstmp," ,");
		while (pch != NULL) {
			gd.legendnames[gd.nlegends] = (string) malloc(30);
			strcpy(gd.legendnames[gd.nlegends],pch);
			++gd.nlegends;
			fprintf(stdout,"%s\n",gd.legendnames[gd.nlegends-1]);
			pch = strtok (NULL, " ,");
		}
		fprintf(stdout,"num. of legends in legends %s =%d\n",cmd.legends,gd.nlegends);
		if (gd.nlegends != gd.nfiles)
			error("\nStartRunParameterfile: nlegends must be equal to number of files\n\n");
	}

	scanbOption(cmd.witherrorbars, gd.errorbars, &gd.nwitherrorbars, gd.nfiles, 1,
		"witherrorbars");

	scanbOption(cmd.plotjoined, gd.plotjoined, &gd.nplotjoined, gd.nfiles, 1,
		"plotjoined");

	scanbOption(cmd.withsymbols, gd.withsymbols, &gd.nwithsymbols, gd.nfiles, 1,
		"withsymbols");
	scanbOption(cmd.withdots, gd.withdots, &gd.nwithdots, gd.nfiles, 1,
		"withdots");

	scaniOption(cmd.linetype, gd.linetype, &gd.nwithsymbols, gd.nfiles, 1,
		"linetype");
	scaniOption(cmd.symboltype, gd.symboltype, &gd.nwithsymbols, gd.nfiles, 1,
		"symboltype");

	if (strcmp(cmd.xrange,"autoscale"))
		gd.x_autoscale=FALSE;

	if (!(sscanf(cmd.xrange, "%lf:%lf", &gd.xmin, &gd.xmax) == 2))
		xrangeflag=0;
	if ( !gd.x_autoscale && xrangeflag==0) 
		error("\nStartRunParameterfile: xrange must be autoscale or in the form xmin:xmax\n\n");
	if (xrangeflag != 0 && gd.xmin >= gd.xmax)
		error("\nStartRunParameterfile: xmin must be < xmax : (xmin,xmax)=(%lf,%lf)\n\n",gd.xmin,gd.xmax);

	if (strcmp(cmd.yrange,"autoscale"))
		gd.y_autoscale=FALSE;
	if (!(sscanf(cmd.yrange, "%lf:%lf", &gd.ymin, &gd.ymax) == 2))
		yrangeflag=0;
	if ( !gd.y_autoscale && yrangeflag==0) 
		error("\nStartRunParameterfile: yrange must be autoscale or in the form ymin:ymax\n\n");
	if (yrangeflag != 0 && gd.ymin >= gd.ymax)
		error("\nStartRunParameterfile: ymin must be < ymax : (ymin,ymax)=(%lf,%lf)\n\n",gd.ymin,gd.ymax);

// ... se usaba antes %ld pero modificaba el valor de gd.x_autoscale, 
//	lo que hacia que no funcionara autoscale ...	
	if (!(sscanf(cmd.graphicsarray, "%dx%d", &gd.graphicsarray_nrows, &gd.graphicsarray_ncols) == 2))
		error("\nStartRunCmdline: graphicsarray must be in the form: nxm\n\n");

	printf("\n\ngraphics array %dx%d\n",gd.graphicsarray_nrows,gd.graphicsarray_ncols);

	if (gd.graphicsarray_nrows<=0 || gd.graphicsarray_ncols<=0)
		error("\nStartRunCmdline: graphicsarray nrows or ncols must be >= 0\n\n");
	if (gd.graphicsarray_nrows!=1 && gd.graphicsarray_ncols!=1)
		if (gd.graphicsarray_nrows*gd.graphicsarray_ncols != gd.nfiles)
			error("\nStartRunCmdline: graphicsarray nrows x ncols must be = nfiles\n\n");
	
	if (!(sscanf(cmd.axesorigin, "%lf,%lf", &gd.xorigin, &gd.yorigin) == 2))
		error("\nStartRunParameterfile: axesorigin must be in the form: xorigin,yorigin\n\n");
	fprintf(stdout,"\n\nOrigin's position is (%g,%g)\n",gd.xorigin,gd.yorigin);

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
		"num. of pairs of colmuns in usingcolumns %s =%d\n",cmd.usingcolumns,npcols);

	if (npcols != gd.nfiles)
		error("\nStartRunParameterfile: nusingcolumns must be equal to number of files\n\n");

	gd.vcol1 = (int *) allocate(npcols*sizeof(int));
	gd.vcol2 = (int *) allocate(npcols*sizeof(int));
	flag = 1;
	for (i=0; i<npcols; i++) {
		if (gd.errorbars[i] && flag==1) {
			flag=0;
			if (cmd.errorbarstype==1) {
				gd.vcol3 = (int *) allocate(npcols*sizeof(int));
				gd.vcol4 = (int *) allocate(npcols*sizeof(int));
			} else {
				gd.vcol3 = (int *) allocate(npcols*sizeof(int));
			}
		}
	}

	for (i=0; i<npcols; i++) {
		if (!strnull(cmd.witherrorbars) && gd.errorbars[i]==1) {
			if (cmd.errorbarstype==1) {
				if (!(sscanf(pusingcolumns[i], "%d:%d:%d:%d", 
					&gd.vcol1[i], &gd.vcol2[i], &gd.vcol3[i], &gd.vcol4[i]) == 4))
					error("\nStartRunCmdline: usingcolumns must be in the form c1:c2:c3:c4\n\n");
				fprintf(stdout,"pairs: %d %d %d %d\n",gd.vcol1[i],gd.vcol2[i],gd.vcol3[i],gd.vcol4[i]);
			} else {
				if (!(sscanf(pusingcolumns[i], "%d:%d:%d", 
					&gd.vcol1[i], &gd.vcol2[i], &gd.vcol3[i]) == 3))
					error("\nStartRunCmdline: usingcolumns must be in the form c1:c2:c3\n\n");
				fprintf(stdout,"pairs: %d %d %d\n",gd.vcol1[i],gd.vcol2[i],gd.vcol3[i]);
			}
		} else {
			if (!(sscanf(pusingcolumns[i], "%d:%d", &gd.vcol1[i], &gd.vcol2[i]) == 2))
				error("\nStartRunCmdline: usingcolumns must be in the form c1:c2\n\n");
			fprintf(stdout,"pairs: %d %d\n",gd.vcol1[i],gd.vcol2[i]);
		}
	}

	if (!strnull(cmd.framestyle)) {
		strcpy(framestyletmp,cmd.framestyle);
		fprintf(stdout,"\nSplitting string \"%s\" in tokens:\n",framestyletmp);
		nitems=0;
		pch = strtok(framestyletmp," ,");
		while (pch != NULL) {
			framestyleitems[nitems] = (string) malloc(10);
			strcpy(framestyleitems[nitems],pch);
			++nitems;
			fprintf(stdout,"%s\n",framestyleitems[nitems-1]);
			pch = strtok (NULL, " ,");
		}
		fprintf(stdout,
			"num. of frame items in framestyle %s = %d\n",cmd.framestyle,nitems);
		if (nitems != 6)
			error("\nStartRunCmdline: nframeitems must be equal to 6\n\n");
		gd.frame_xopt = (string) malloc(15);
		gd.frame_yopt = (string) malloc(15);
		strcpy(gd.frame_xopt,framestyleitems[0]);
		strcpy(gd.frame_yopt,framestyleitems[3]);
		if (!(sscanf(framestyleitems[1], "%lf", &gd.frame_xtick) == 1))
			error("\nStartRunCmdline: error getting frame_xtick\n\n");
		if (!(sscanf(framestyleitems[4], "%lf", &gd.frame_ytick) == 1))
			error("\nStartRunCmdline: error getting frame_ytick\n\n");
		if (!(sscanf(framestyleitems[2], "%d", &gd.frame_nxsub) == 1))
			error("\nStartRunCmdline: error getting frame_nxsub\n\n");
		if (!(sscanf(framestyleitems[5], "%d", &gd.frame_nysub) == 1))
			error("\nStartRunCmdline: error getting frame_nysub\n\n");
		fprintf(stdout,"\n\nframestyleitems parameters: %s %g %d %s %g %d\n\n",
			gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,
			gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
	}

	if (cmd.plottype<0 || cmd.plottype >3)
		error("\nStartRunParameterfile: plottype must be 0, 1, 2 or 3\n\n");

	PrintParametersLog();
}

local void startrun_ParamStat(void)
{
	if (GetParamStat("xlabel") & ARGPARAM)
		cmd.xlabel = GetParam("xlabel");
	if (GetParamStat("ylabel") & ARGPARAM)
		cmd.ylabel = GetParam("ylabel");
	if (GetParamStat("plotlabel") & ARGPARAM)
		cmd.plotlabel = GetParam("plotlabel");

	if (GetParamStat("gridlines") & ARGPARAM)
		cmd.gridlines = GetbParam("gridlines");


	if (GetParamStat("withsymbols") & ARGPARAM)
		cmd.withsymbols = GetParam("withsymbols");
	if (GetParamStat("symbolsize") & ARGPARAM)
		cmd.symbolsize = GetdParam("symbolsize");
	if (GetParamStat("symbolcolor") & ARGPARAM)
		cmd.symbolcolor = GetiParam("symbolcolor");
}

local void CheckParameters(void)
{
}


local void PrintParametersLog(void)
{
	fprintf(gd.outlog,"\nParameters:\n");
	fprintf(gd.outlog,"x_autoscale= %s",gd.x_autoscale ? "true" : "false");
	fprintf(gd.outlog,"y_autoscale= %s",gd.y_autoscale ? "true" : "false");
	fprintf(gd.outlog,"\nxmin= %g",gd.xmin);
	fprintf(gd.outlog,"\nxmax= %g",gd.xmax);
	fprintf(gd.outlog,"\nymin= %g",gd.ymin);
	fprintf(gd.outlog,"\nymax= %g",gd.ymax);
	fflush(gd.outlog);
}

#undef logfile
#undef parameter_null

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

	SPName(cmd.inputfile,"inputfile",200);
	SPName(cmd.legends,"legends",100);
	IPName(cmd.legendspos,"legendspos");
	SPName(cmd.witherrorbars,"witherrorbars",100);
	IPName(cmd.errorbarstype,"errorbarstype");
//	SPName(cmd.outputformat,"outputformat",100);
	SPName(cmd.axesorigin,"axesorigin",100);
	SPName(cmd.xrange,"xrange",1);
	SPName(cmd.yrange,"yrange",1);
	SPName(cmd.usingcolumns,"usingcolumns",100);
	SPName(cmd.xlabel,"xlabel",1);
	SPName(cmd.ylabel,"ylabel",1);
	SPName(cmd.graphicsarray,"graphicsarray",1);

	IPName(cmd.plottype,"plottype");
	IPName(cmd.labelcolor,"labelcolor");
	IPName(cmd.fontset,"fontset");
	IPName(cmd.fontchr,"fontchr");
	RPName(cmd.legendsfontsize,"legendsfontsize");
	IPName(cmd.legendsfontweight,"legendsfontweight");
	RPName(cmd.labelfontsize,"labelfontsize");
	IPName(cmd.labelfontweight,"labelfontweight");
	BPName(cmd.xaxis,"xaxis");
	BPName(cmd.yaxis,"yaxis");
	IPName(cmd.axestype,"axestype");
	IPName(cmd.axeswidth,"axeswidth");
	IPName(cmd.axescolor,"axescolor");
	SPName(cmd.plotlabel,"plotlabel",100);
	IPName(cmd.nlxdigmax,"nlxdigmax");
	IPName(cmd.nlydigmax,"nlydigmax");
	RPName(cmd.nlsize,"nlsize");
	IPName(cmd.nlcolor,"nlcolor");
	SPName(cmd.plotjoined,"plotjoined",100);
	SPName(cmd.linetype,"linetype",100);
	IPName(cmd.linewidth,"linewidth");
	IPName(cmd.linecolor,"linecolor");
	SPName(cmd.withsymbols,"withsymbols",100);
	SPName(cmd.withdots,"withdots",100);
	SPName(cmd.symboltype,"symboltype",100);
	IPName(cmd.symbolcolor,"symbolcolor");
	IPName(cmd.symbolweight,"symbolweight");
	RPName(cmd.symbolsize,"symbolsize");

	SPName(cmd.text1,"text1",100);
	SPName(cmd.text1side,"text1side",1);
	RPName(cmd.text1disp,"text1disp");
	RPName(cmd.text1pos,"text1pos");
	RPName(cmd.text1just,"text1just");
	RPName(cmd.text1size,"text1size");
	IPName(cmd.text1weight,"text1weight");
	IPName(cmd.text1color,"text1color");

	SPName(cmd.text2,"text2",100);
	SPName(cmd.text2side,"text2side",1);
	RPName(cmd.text2disp,"text2disp");
	RPName(cmd.text2pos,"text2pos");
	RPName(cmd.text2just,"text2just");
	RPName(cmd.text2size,"text2size");
	IPName(cmd.text2weight,"text2weight");
	IPName(cmd.text2color,"text2color");

	SPName(cmd.text3,"text3",100);
	SPName(cmd.text3side,"text3side",1);
	RPName(cmd.text3disp,"text3disp");
	RPName(cmd.text3pos,"text3pos");
	RPName(cmd.text3just,"text3just");
	RPName(cmd.text3size,"text3size");
	IPName(cmd.text3weight,"text3weight");
	IPName(cmd.text3color,"text3color");

	SPName(cmd.text4,"text4",100);
	RPName(cmd.text4x,"text4x");
	RPName(cmd.text4y,"text4y");
	RPName(cmd.text4dx,"text4dx");
	RPName(cmd.text4dy,"text4dy");
	RPName(cmd.text4just,"text4just");
	RPName(cmd.text4size,"text4size");
	IPName(cmd.text4weight,"text4weight");
	IPName(cmd.text4color,"text4color");

	BPName(cmd.locatemode,"locatemode");

	BPName(cmd.frame,"frame");
	SPName(cmd.framestyle,"framestyle",50);
	BPName(cmd.gridlines,"gridlines");

	SPName(cmd.epilog,"epilog",100);

//---------MENU COMMAND LINE PLPLOT OPTIONS------------------------------------
// Option set to false == use default, they are bool

	BPName(cmd.pl_showall,"showall");
	BPName(cmd.pl_h,"h");
	BPName(cmd.pl_v,"v");
	BPName(cmd.pl_verbose,"verbose");
	BPName(cmd.pl_debug,"debug");
	BPName(cmd.pl_hack,"hack");
	SPName(cmd.pl_dev,"dev",1);
	SPName(cmd.pl_o,"outputfile",100);
	SPName(cmd.pl_display,"display",1);
	SPName(cmd.pl_geometry,"geometry",50);
	SPName(cmd.pl_wplt,"wplt",1);
	SPName(cmd.pl_mar,"mar",1);
	SPName(cmd.pl_a,"aspectratio",50);
	SPName(cmd.pl_jx,"jx",1);
	SPName(cmd.pl_jy,"jy",1);
	SPName(cmd.pl_ori,"ori",1);
	BPName(cmd.pl_freeaspect,"freeaspect");
	BPName(cmd.pl_portrait,"portrait");
	SPName(cmd.pl_width,"width",1);
	SPName(cmd.pl_bg,"background",50);
	SPName(cmd.pl_ncol0,"ncol0",1);
	SPName(cmd.pl_ncol1,"ncol1",1);
	BPName(cmd.pl_fam,"fam");
	SPName(cmd.pl_fsiz,"fsiz",1);
	SPName(cmd.pl_fbeg,"fbeg",1);
	SPName(cmd.pl_finc,"finc",1);
	SPName(cmd.pl_fflen,"fflen",1);
	BPName(cmd.pl_nopixmap,"nopixmap");
	BPName(cmd.pl_db,"db");
	BPName(cmd.pl_np,"np");
	SPName(cmd.pl_bufmax,"bufmax",1);
//	SPName(cmd.pl_server_name,"server_name",1);
//	SPName(cmd.pl_plserver,"plserver",1);
//	SPName(cmd.pl_plwindow,"plwindow",1);
//	SPName(cmd.pl_tcl_cmd,"tcl_cmd",1);
//	SPName(cmd.pl_auto_path,"auto_path",1);
//	SPName(cmd.pl_tk_file,"tk_file",1);
	SPName(cmd.pl_dpi,"dpi",1);
	SPName(cmd.pl_compression,"compression",1);
	SPName(cmd.pl_drvopt,"drvopt",1);

//------------------------------------------------------------------------------

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
//                *buf2=(char)NULL;
// Removing the warning:
// warning: cast from pointer to integer of different size [-Wpointer-to-int-cast]
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
                                error("GetbParam: %s=%s not bool\n", buf1, buf2);
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
    
    sprintf(buf,"%s%s",fname,"-usedvalues");
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

        fprintf(fdout,FMTT,"inputfile",cmd.inputfile);
//        fprintf(fdout,FMTT,"outputformat",cmd.outputformat);
        fprintf(fdout,FMTT,"xrange",cmd.xrange);
        fprintf(fdout,FMTT,"yrange",cmd.yrange);
        fprintf(fdout,FMTT,"usingcolumns",cmd.usingcolumns);
        fprintf(fdout,FMTT,"xlabel",cmd.xlabel);
        fprintf(fdout,FMTT,"ylabel",cmd.ylabel);
        fprintf(fdout,FMTI,"labelcolor",cmd.labelcolor);

        fprintf(fdout,FMTT,"legends",cmd.legends);
        fprintf(fdout,FMTI,"legendspos",cmd.legendspos);
        fprintf(fdout,FMTR,"legendsfontsize",cmd.legendsfontsize);
        fprintf(fdout,FMTI,"legendsfontweight",cmd.legendsfontweight);

        fprintf(fdout,FMTT,"witherrorbars",cmd.witherrorbars);
        fprintf(fdout,FMTI,"errorbarstype",cmd.errorbarstype);
        fprintf(fdout,FMTI,"fontset",cmd.fontset);
        fprintf(fdout,FMTI,"fontchr",cmd.fontchr);
        fprintf(fdout,FMTR,"labelfontsize",cmd.labelfontsize);
        fprintf(fdout,FMTI,"labelfontweight",cmd.labelfontweight);
        fprintf(fdout,FMTT,"xaxis",cmd.xaxis ? "true" : "false");
        fprintf(fdout,FMTT,"yaxis",cmd.yaxis ? "true" : "false");

        fprintf(fdout,FMTT,"graphicsarray",cmd.graphicsarray);
        fprintf(fdout,FMTI,"plottype",cmd.plottype);

        fprintf(fdout,FMTT,"axesorigin",cmd.axesorigin);

        fprintf(fdout,FMTI,"axestype",cmd.axestype);
        fprintf(fdout,FMTI,"axeswidth",cmd.axeswidth);
        fprintf(fdout,FMTI,"axescolor",cmd.axescolor);
        fprintf(fdout,FMTT,"plotlabel",cmd.plotlabel);

        fprintf(fdout,FMTI,"nlxdigmax",cmd.nlxdigmax);
        fprintf(fdout,FMTI,"nlydigmax",cmd.nlydigmax);
        fprintf(fdout,FMTR,"nlsize",cmd.nlsize);
        fprintf(fdout,FMTI,"nlcolor",cmd.nlcolor);

        fprintf(fdout,FMTT,"plotjoined",cmd.plotjoined);
        fprintf(fdout,FMTT,"linetype",cmd.linetype);
        fprintf(fdout,FMTI,"linewidth",cmd.linewidth);
        fprintf(fdout,FMTI,"linecolor",cmd.linecolor);
        fprintf(fdout,FMTT,"withdots",cmd.withdots);
        fprintf(fdout,FMTT,"withsymbols",cmd.withsymbols);
        fprintf(fdout,FMTT,"symboltype",cmd.symboltype);
        fprintf(fdout,FMTI,"symbolcolor",cmd.symbolcolor);
        fprintf(fdout,FMTI,"symbolweight",cmd.symbolweight);
        fprintf(fdout,FMTR,"symbolsize",cmd.symbolsize);

        fprintf(fdout,FMTT,"text1",cmd.text1);
        fprintf(fdout,FMTT,"text1side",cmd.text1side);
        fprintf(fdout,FMTR,"text1disp",cmd.text1disp);
        fprintf(fdout,FMTR,"text1pos",cmd.text1pos);
        fprintf(fdout,FMTR,"text1just",cmd.text1just);
        fprintf(fdout,FMTR,"text1size",cmd.text1size);
        fprintf(fdout,FMTI,"text1weight",cmd.text1weight);
        fprintf(fdout,FMTI,"text1color",cmd.text1color);

        fprintf(fdout,FMTT,"text2",cmd.text2);
        fprintf(fdout,FMTT,"text2side",cmd.text2side);
        fprintf(fdout,FMTR,"text2disp",cmd.text2disp);
        fprintf(fdout,FMTR,"text2pos",cmd.text2pos);
        fprintf(fdout,FMTR,"text2just",cmd.text2just);
        fprintf(fdout,FMTR,"text2size",cmd.text2size);
        fprintf(fdout,FMTI,"text2weight",cmd.text2weight);
        fprintf(fdout,FMTI,"text2color",cmd.text2color);

        fprintf(fdout,FMTT,"text3",cmd.text3);
        fprintf(fdout,FMTT,"text3side",cmd.text3side);
        fprintf(fdout,FMTR,"text3disp",cmd.text3disp);
        fprintf(fdout,FMTR,"text3pos",cmd.text3pos);
        fprintf(fdout,FMTR,"text3just",cmd.text3just);
        fprintf(fdout,FMTR,"text3size",cmd.text3size);
        fprintf(fdout,FMTI,"text3weight",cmd.text3weight);
        fprintf(fdout,FMTI,"text3color",cmd.text3color);

        fprintf(fdout,FMTT,"text4",cmd.text4);
        fprintf(fdout,FMTR,"text4x",cmd.text4x);
        fprintf(fdout,FMTR,"text4y",cmd.text4y);
        fprintf(fdout,FMTR,"text4dx",cmd.text4dx);
        fprintf(fdout,FMTR,"text4dy",cmd.text4dy);
        fprintf(fdout,FMTR,"text4just",cmd.text4just);
        fprintf(fdout,FMTR,"text4size",cmd.text4size);
        fprintf(fdout,FMTI,"text4weight",cmd.text4weight);
        fprintf(fdout,FMTI,"text4color",cmd.text4color);

        fprintf(fdout,FMTT,"locatemode",cmd.locatemode ? "true" : "false");

        fprintf(fdout,FMTT,"frame",cmd.frame ? "true" : "false");
        fprintf(fdout,FMTT,"framestyle",cmd.framestyle);
        fprintf(fdout,FMTT,"gridlines",cmd.gridlines ? "true" : "false");
        fprintf(fdout,FMTT,"epilog",cmd.epilog);

//---------MENU COMMAND LINE PLPLOT OPTIONS------------------------------------
// Option set to false == use default, they are bool

        fprintf(fdout,FMTT,"showall",cmd.pl_showall ? "true" : "false");
        fprintf(fdout,FMTT,"h",cmd.pl_h ? "true" : "false");
        fprintf(fdout,FMTT,"v",cmd.pl_v ? "true" : "false");
        fprintf(fdout,FMTT,"verbose",cmd.pl_verbose ? "true" : "false");
        fprintf(fdout,FMTT,"debug",cmd.pl_debug ? "true" : "false");
        fprintf(fdout,FMTT,"hack",cmd.pl_hack ? "true" : "false");
        fprintf(fdout,FMTT,"dev",cmd.pl_dev);
        fprintf(fdout,FMTT,"outputfile",cmd.pl_o);
        fprintf(fdout,FMTT,"display",cmd.pl_display);
        fprintf(fdout,FMTT,"geometry",cmd.pl_geometry);
        fprintf(fdout,FMTT,"wplt",cmd.pl_wplt);
        fprintf(fdout,FMTT,"mar",cmd.pl_mar);
        fprintf(fdout,FMTT,"aspectratio",cmd.pl_a);
        fprintf(fdout,FMTT,"jx",cmd.pl_jx);
        fprintf(fdout,FMTT,"jy",cmd.pl_jy);
        fprintf(fdout,FMTT,"ori",cmd.pl_ori);
        fprintf(fdout,FMTT,"freeaspect",cmd.pl_freeaspect ? "true" : "false");
        fprintf(fdout,FMTT,"portrait",cmd.pl_portrait ? "true" : "false");
        fprintf(fdout,FMTT,"width",cmd.pl_width);
        fprintf(fdout,FMTT,"background",cmd.pl_bg);
        fprintf(fdout,FMTT,"ncol0",cmd.pl_ncol0);
        fprintf(fdout,FMTT,"ncol1",cmd.pl_ncol1);
        fprintf(fdout,FMTT,"fam",cmd.pl_fam ? "true" : "false");
        fprintf(fdout,FMTT,"fsiz",cmd.pl_fsiz);
        fprintf(fdout,FMTT,"fbeg",cmd.pl_fbeg);
        fprintf(fdout,FMTT,"finc",cmd.pl_finc);
        fprintf(fdout,FMTT,"fflen",cmd.pl_fflen);
        fprintf(fdout,FMTT,"nopixmap",cmd.pl_nopixmap ? "true" : "false");
        fprintf(fdout,FMTT,"db",cmd.pl_db ? "true" : "false");
        fprintf(fdout,FMTT,"np",cmd.pl_np ? "true" : "false");
        fprintf(fdout,FMTT,"bufmax",cmd.pl_bufmax);
//        fprintf(fdout,FMTT,"server_name",cmd.pl_server_name);
//        fprintf(fdout,FMTT,"plserver",cmd.pl_plserver);
//        fprintf(fdout,FMTT,"plwindow",cmd.pl_plwindow);
//        fprintf(fdout,FMTT,"tcl_cmd",cmd.pl_tcl_cmd);
//        fprintf(fdout,FMTT,"auto_path",cmd.pl_auto_path);
//        fprintf(fdout,FMTT,"tk_file",cmd.pl_tk_file);
        fprintf(fdout,FMTT,"dpi",cmd.pl_dpi);
        fprintf(fdout,FMTT,"compression",cmd.pl_compression);
        fprintf(fdout,FMTT,"drvopt",cmd.pl_drvopt);

//------------------------------------------------------------------------------

        fprintf(fdout,"\n\n");					// to read as input paramfile

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


local void scaniOption(string optionstr, int *option, int *noption, 
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
			option[i]=atoi(poptionstr[i]);
			fprintf(stdout,"option: %d\n",option[i]);
		}
		fprintf(stdout,"\nnoptions, nfiles: %d %d\n",*noption,nfiles);
		if (flag == 1) {
			if (*noption > nfiles)
				error("\nscanOption: noption = %d must be less or equal to number of files\n\n",*noption);
			else {
				for (i=*noption; i<nfiles; i++) {
					option[i]=option[i-1]+1;
					fprintf(stdout,"option: %d\n",option[i]);
				}
			}
		}
	} else {
		for (i=0; i<nfiles; i++) {
			option[i]=1;
			fprintf(stdout,"option: %d\n",option[i]);
		}
	}
}

