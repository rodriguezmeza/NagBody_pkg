/* =============================================================================
	MODULE: startrun.c			[galaxy_models]
	Written by: M.A. Rodriguez-Meza
	Fecha de creaci'on: May 2006
	Purpose: routines to initialize the main code
	Language: C
	Use: 'Start_Run();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: June 6, 2008;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved.
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrunCommon(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

local void SetUnits(void);
local void SetArrays(void);

void StartRun(string head0, string head1, string head2, string head3)
{
    double cpustart;

    cpustart = cputime();                       

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n", 
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);
    printf("\t-- Version %s --\n",getversion());
	gd.comment = "BDH Galaxy IC";
    printf("  \t -- %s --\n", gd.comment);
    printf("\n%s\n", copyright);

//    LicDriver(gd.headline0);

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

	StartOutput();

	fprintf(stdout,"\n\nStartRun CPU time: %g\n",cputime()-cpustart);
}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
	startrun_ParamStat();
	startrunCommon();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-galaxy_model"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrunCommon();
	PrintParameterFile(parameter_null);
}

#undef parameter_null

local void ReadParametersCmdline(void)
{
	cmd.outpteps = GetbParam("outpteps");
	cmd.filename = GetParam("filename");
	cmd.filenamefmt = GetParam("filenamefmt");
	cmd.statfile = GetParam("statfile");
	cmd.seed = GetiParam("seed");         

// ----Input bulge's parameters--------------
	cmd.usebulge = GetbParam("usebulge");
	cmd.nbulge = GetiParam("nbulge");
	cmd.bulgmass = GetdParam("bulgmass");
	cmd.abulge = GetdParam("abulge");
	cmd.selfgbul = GetbParam("selfgbul");
	cmd.rmaxbulg = GetdParam("rmaxbulg");
	cmd.epsbulge = GetdParam("epsbulge");
	cmd.axibulge = GetbParam("axibulge");
	cmd.cbulge = GetdParam("cbulge");
	cmd.zmaxbulg = GetdParam("zmaxbulg");
	cmd.nsimpson = GetiParam("nsimpson");
	cmd.bulgerot = GetbParam("bulgerot");
	cmd.brotfrac = GetdParam("brotfrac");

// ----Input disk's parameters--------------
	cmd.usedisk = GetbParam("usedisk");
	cmd.ndstars = GetiParam("ndstars");
	cmd.z0 = GetdParam("z0");
	cmd.rsolarstr = GetParam("rsolar");
	cmd.qsolar = GetdParam("qsolar");
	cmd.epsdisk = GetdParam("epsdisk");
	cmd.zmax = GetdParam("zmax");
	cmd.rmax = GetdParam("rmax");
	cmd.usegas = GetbParam("usegas");
	cmd.ndgas = GetiParam("ndgas");
	cmd.gasmass = GetbParam("gasmass");
	cmd.gastemp = GetdParam("gastemp");
	cmd.z0gas = GetdParam("z0gas");
	cmd.zmaxgas = GetdParam("zmaxgas");
	cmd.rmaxgas = GetdParam("rmaxgas");
	cmd.rmingas = GetdParam("rmingas");
	cmd.selfggas = GetbParam("selfggas");

// ----Input halo's parameters--------------
	cmd.usehalo = GetbParam("usehalo");
	cmd.nhalo = GetiParam("nhalo");
	cmd.halomass = GetdParam("halomass");
	cmd.gamhalo = GetdParam("gamhalo");
	cmd.selfghal = GetbParam("selfghal");
	cmd.rmaxhalo = GetdParam("rmaxhalo");
	cmd.epshalo = GetdParam("epshalo");
	cmd.halotype = GetParam("halotype");
	cmd.rthalo = GetdParam("rthalo");
	cmd.ahalo = GetdParam("ahalo");

// ----Input satellite's parameters----------
	cmd.usesat = GetbParam("usesat");
	cmd.nsat = GetiParam("nsat");
	cmd.satmass = GetdParam("satmass");
	cmd.selfgsat = GetbParam("selfgsat");
	cmd.rmaxsat = GetdParam("rmaxsat");
	cmd.epssat = GetdParam("epssat");
	cmd.asat = GetdParam("asat");
	cmd.xsat = GetdParam("xsat");
	cmd.ysat = GetdParam("ysat");
	cmd.zsat = GetdParam("zsat");
	cmd.vxsat = GetdParam("vxsat");
	cmd.vysat = GetdParam("vysat");
	cmd.vzsat = GetdParam("vzsat");

// ----Input addtwomod's parameters----------
	cmd.addmods = GetbParam("addmods");
	cmd.rp = GetdParam("rp");
	cmd.rsep = GetdParam("rsep");
	cmd.thetmod1 = GetdParam("thetmod1");
	cmd.phimod1 = GetdParam("phimod1");
	cmd.thetmod2 = GetdParam("thetmod2");
	cmd.phimod2 = GetdParam("phimod2");
}

#define logfile			"galaxy_models.log"

local void startrunCommon(void)
{
	real r1, r2;

    if(!(gd.outlog=fopen(logfile,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",logfile);
        exit(0);
    }

	idum=cmd.seed;							// Save the initial random seed
	xsrandom(idum);

	SetUnits();
	gd.nbodies=0;

// -----------------------Init bulge's parameters----------------------
	if (cmd.usebulge) {
		if (cmd.selfgbul) {
			if (cmd.axibulge) {
				if (cmd.cbulge > cmd.abulge)
					error("\n\ncbulge must be <= abulge\n");
				if (cmd.bulgerot) {
					if (cmd.brotfrac < 0.5 || cmd.brotfrac > 1.0)
						error("\n\nFrac. must be between 0.5 and 1\n");
				} else {
					cmd.bulgerot = FALSE;
					cmd.brotfrac = 0.0;
				}
			} else {
				cmd.bulgerot = FALSE;
				cmd.brotfrac = 0.0;
			}
		} else {
			cmd.nbulge = 0;
			cmd.axibulge = FALSE;
			cmd.bulgerot = FALSE;
			cmd.brotfrac = 0.0;
		}
	} else {
		cmd.nbulge = 0;
		cmd.bulgmass = 0.0;
	}

	gd.nbodies += cmd.nbulge;
	if (gd.nbodies<0 || gd.nbodies> nbodsmax)
		error("\n\nnbulge: Invalid number of bodies\n");

// -----------------------Init disk's parameters----------------------

	gd.rsolar = (sscanf(cmd.rsolarstr, "%lf/%lf", &r1, &r2) == 2 ?
			r1/r2 : atof(cmd.rsolarstr));
	if ( r2 == 0. )
		error("\n\nstartrunCommon: rsolar : r2 must be finite\n");

	if (!cmd.usedisk) {						// Checar que los valores en estos parametros
		cmd.ndstars=0;						// estan bien...
		cmd.zmax=0.0;
		cmd.rmax=0.0;
	}else{
		if (cmd.ndstars<=0) error("ndstars: Invalid number of bodies\n");
	}

	gd.nbodies += cmd.ndstars;
	if (gd.nbodies<0 || gd.nbodies> nbodsmax)
		error("\n\nndstars: Invalid number of bodies\n");

	gd.epsdisk2=cmd.epsdisk*cmd.epsdisk;

	if (!cmd.usegas) {
		cmd.ndgas=0;
		cmd.gasmass=0.0;
		cmd.zmaxgas=0.0;
		cmd.rmaxgas=cmd.rmax;
		cmd.rmingas=0.0;
	} else {
		if (cmd.ndgas<=0) error("ndgas: Invalid number of bodies\n");
	}

	gd.nbodies += cmd.ndgas;
	if (gd.nbodies<0 || gd.nbodies> nbodsmax)
		error("\n\nndgas: Invalid number of bodies\n");

	gd.ndisk=cmd.ndgas+cmd.ndstars;

// -----------------------Init halo's parameters----------------------

	if (!cmd.usehalo) {						// Checar que los valores en estos parametros
		cmd.nhalo=0;						// estan bien...
		cmd.halomass=0.0;
	}else{
		if (!cmd.selfghal)
			cmd.nhalo=0;

		if (cmd.nhalo<=0) error("nhalo: Invalid number of bodies\n");

		if (scanopt(cmd.halotype, "IS")) {
printf("\nHalo type is IS\n");
		} else {
printf("\nHalo type must be LH\n");
			if (!scanopt(cmd.halotype, "LH")) {
				error("halo: Halo type error\n");
			} else { 
//				ahalo=0.0;				// Checar que no se usa este parametro en este caso...
			}
		}
	}

	gd.nbodies += cmd.nhalo;
	if (gd.nbodies<0 || gd.nbodies> nbodsmax)
		error("\n\nnhalo: Invalid number of bodies\n");

// -----------------------Init satellite's parameters----------------------

	if (!cmd.usesat) {					// Checar que los valores en estos parametros
		cmd.nsat=0;						// estan bien...
		cmd.satmass=0.0;
	}else{
		if (!cmd.selfgsat) {
              cmd.nsat = 0; 
              error("\nRigid satellite option not implemented"); 
		}

		if (cmd.nsat<=0) error("nsat: Invalid number of bodies\n");
	}

	gd.nbodies += cmd.nsat;
	if (gd.nbodies<0 || gd.nbodies> nbodsmax)
		error("\n\nnsat: Invalid number of bodies\n");

	printf("\n\nnbodies, ndgas, ndstars, ndisk, nbulge, nhalo, nsat : %d %d %d %d %d %d %d\n\n",
	gd.nbodies, cmd.ndgas, cmd.ndstars, gd.ndisk, cmd.nbulge, cmd.nhalo, cmd.nsat);

	printf("\ngasmass, bulgmass, halomass, satmass : %g %g %g %g\n\n",
	cmd.gasmass, cmd.bulgmass, cmd.halomass, cmd.satmass);

	printf("\n\n zmaxgas, rmaxgas, rmingas, rmax: %g %g %g %g\n\n",
	cmd.zmaxgas, cmd.rmaxgas, cmd.rmingas, cmd.rmax);

	SetArrays();
}

#undef logfile

local void startrun_ParamStat(void)
{
// INCLUIR LOS PARAMETROS QUE PUEDAN SER UTILES DE AGREGAR 
// DESPUES DE UN ARCHIVO DE PARAMETROS ...
}

local void CheckParameters(void)
{
	if (strnull(cmd.filename)) 
		error("startrun: give a file name for the output of the galaxy model\n");
	if (strnull(cmd.statfile)) 
		error("startrun: give a file name for the statfile of the galaxy model\n");

	if (cmd.nbulge < 1)
		error("\nCheckParameters: absurd value for nbulge\n\n");
	if (cmd.ndstars < 1)
		error("\nCheckParameters: absurd value for ndstars\n\n");
	if (cmd.ndgas < 1)
		error("\nCheckParameters: absurd value for ndgas\n\n");
	if (cmd.nhalo < 1)
		error("\nCheckParameters: absurd value for nhalo\n\n");
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

	SPName(cmd.filename,"filename",100);
//	strcpy(tag[nt],"filename"); 
//  cmd.filename=(string) malloc(100);
//  addr[nt]=cmd.filename;
//  id[nt++]=STRING;

	SPName(cmd.filenamefmt,"filenamefmt",100);
//  strcpy(tag[nt],"filenamefmt"); 
//  cmd.filenamefmt=(string) malloc(100);
//  addr[nt]=cmd.filenamefmt;
//  id[nt++]=STRING;

	SPName(cmd.statfile,"statfile",100);
//  strcpy(tag[nt],"statfile"); 
//  cmd.statfile=(string) malloc(100);
//  addr[nt]=cmd.statfile;
//  id[nt++]=STRING;

 	BPName(cmd.outpteps,"outpteps");
//  strcpy(tag[nt],"outpteps"); 
//  addr[nt]=&cmd.outpteps;
//  id[nt++]=BOOLEAN;

// ----Input bulge's parameters--------------
 	BPName(cmd.usebulge,"usebulge");
//  strcpy(tag[nt],"usebulge");
//  addr[nt]=&cmd.usebulge;
//  id[nt++]=BOOLEAN;

 	IPName(cmd.nbulge,"nbulge");
//  strcpy(tag[nt],"nbulge");
//  addr[nt]=&cmd.nbulge;
//  id[nt++]=INT;

	RPName(cmd.bulgmass,"bulgmass");
//  strcpy(tag[nt],"bulgmass");
//  addr[nt]=&cmd.bulgmass;
//  id[nt++]=DOUBLE;

	RPName(cmd.abulge,"abulge");
//  strcpy(tag[nt],"abulge");
//  addr[nt]=&cmd.abulge;
//  id[nt++]=DOUBLE;

 	BPName(cmd.selfgbul,"selfgbul");
//  strcpy(tag[nt],"selfgbul");
//  addr[nt]=&cmd.selfgbul;
//  id[nt++]=BOOLEAN;

	RPName(cmd.rmaxbulg,"rmaxbulg");
//  strcpy(tag[nt],"rmaxbulg");
//  addr[nt]=&cmd.rmaxbulg;
//  id[nt++]=DOUBLE;

	RPName(cmd.epsbulge,"epsbulge");
//  strcpy(tag[nt],"epsbulge");
//  addr[nt]=&cmd.epsbulge;
//  id[nt++]=DOUBLE;

 	BPName(cmd.axibulge,"axibulge");
//  strcpy(tag[nt],"axibulge");
//  addr[nt]=&cmd.axibulge;
//  id[nt++]=BOOLEAN;

	RPName(cmd.cbulge,"cbulge");
//  strcpy(tag[nt],"cbulge");
//  addr[nt]=&cmd.cbulge;
//  id[nt++]=DOUBLE;

	RPName(cmd.zmaxbulg,"zmaxbulg");
//  strcpy(tag[nt],"zmaxbulg");
//  addr[nt]=&cmd.zmaxbulg;
//  id[nt++]=DOUBLE;

 	IPName(cmd.nsimpson,"nsimpson");
//  strcpy(tag[nt],"nsimpson");
//  addr[nt]=&cmd.nsimpson;
//  id[nt++]=INT;

 	BPName(cmd.bulgerot,"bulgerot");
//  strcpy(tag[nt],"bulgerot");
//  addr[nt]=&cmd.bulgerot;
//  id[nt++]=BOOLEAN;

	RPName(cmd.brotfrac,"brotfrac");
//  strcpy(tag[nt],"brotfrac");
//  addr[nt]=&cmd.brotfrac;
//  id[nt++]=DOUBLE;

// ----Input disk's parameters--------------
 	BPName(cmd.usedisk,"usedisk");
//  strcpy(tag[nt],"usedisk");
//  addr[nt]=&cmd.usedisk;
//  id[nt++]=BOOLEAN;

 	IPName(cmd.ndstars,"ndstars");
//  strcpy(tag[nt],"ndstars");
//  addr[nt]=&cmd.ndstars;
//  id[nt++]=INT;

	RPName(cmd.z0,"z0");
//  strcpy(tag[nt],"z0");
//  addr[nt]=&cmd.z0;
//  id[nt++]=DOUBLE;

	SPName(cmd.rsolarstr,"rsolar",100);
//  strcpy(tag[nt],"rsolar");
//  addr[nt]=&cmd.rsolar;
//  id[nt++]=DOUBLE;

	RPName(cmd.qsolar,"qsolar");
//  strcpy(tag[nt],"qsolar");
//  addr[nt]=&cmd.qsolar;
//  id[nt++]=DOUBLE;

	RPName(cmd.zmax,"zmax");
//  strcpy(tag[nt],"zmax");
//  addr[nt]=&cmd.zmax;
//  id[nt++]=DOUBLE;

	RPName(cmd.rmax,"rmax");
//  strcpy(tag[nt],"rmax");
//  addr[nt]=&cmd.rmax;
//  id[nt++]=DOUBLE;

	RPName(cmd.epsdisk,"epsdisk");
//  strcpy(tag[nt],"epsdisk");
//  addr[nt]=&cmd.epsdisk;
//  id[nt++]=DOUBLE;

 	BPName(cmd.usegas,"usegas");
//  strcpy(tag[nt],"usegas");
//  addr[nt]=&cmd.usegas;
//  id[nt++]=BOOLEAN;

 	IPName(cmd.ndgas,"ndgas");
//  strcpy(tag[nt],"ndgas");
//  addr[nt]=&cmd.ndgas;
//  id[nt++]=INT;

 	BPName(cmd.selfggas,"selfggas");
//  strcpy(tag[nt],"selfggas");
//  addr[nt]=&cmd.selfggas;
//  id[nt++]=BOOLEAN;

	RPName(cmd.gasmass,"gasmass");
//  strcpy(tag[nt],"gasmass");
//  addr[nt]=&cmd.gasmass;
//  id[nt++]=DOUBLE;

	RPName(cmd.gastemp,"gastemp");
//  strcpy(tag[nt],"gastemp");
//  addr[nt]=&cmd.gastemp;
//  id[nt++]=DOUBLE;

	RPName(cmd.z0gas,"z0gas");
//  strcpy(tag[nt],"z0gas");
//  addr[nt]=&cmd.z0gas;
//  id[nt++]=DOUBLE;

	RPName(cmd.zmaxgas,"zmaxgas");
//  strcpy(tag[nt],"zmaxgas");
//  addr[nt]=&cmd.zmaxgas;
//  id[nt++]=DOUBLE;

	RPName(cmd.rmaxgas,"rmaxgas");
//  strcpy(tag[nt],"rmaxgas");
//  addr[nt]=&cmd.rmaxgas;
//  id[nt++]=DOUBLE;

	RPName(cmd.rmingas,"rmingas");
//  strcpy(tag[nt],"rmingas");
//  addr[nt]=&cmd.rmingas;
//  id[nt++]=DOUBLE;

// ----Input halo's parameters--------------
 	BPName(cmd.usehalo,"usehalo");
//  strcpy(tag[nt],"usehalo");
//  addr[nt]=&cmd.usehalo;
//  id[nt++]=BOOLEAN;

 	IPName(cmd.nhalo,"nhalo");
//  strcpy(tag[nt],"nhalo");
//  addr[nt]=&cmd.nhalo;
//  id[nt++]=INT;

	RPName(cmd.halomass,"halomass");
//  strcpy(tag[nt],"halomass");
//  addr[nt]=&cmd.halomass;
//  id[nt++]=DOUBLE;

	RPName(cmd.gamhalo,"gamhalo");
//  strcpy(tag[nt],"gamhalo");
//  addr[nt]=&cmd.gamhalo;
//  id[nt++]=DOUBLE;

 	BPName(cmd.selfghal,"selfghal");
//  strcpy(tag[nt],"selfghal");
//  addr[nt]=&cmd.selfghal;
//  id[nt++]=BOOLEAN;

	RPName(cmd.rthalo,"rthalo");
//  strcpy(tag[nt],"rthalo");
//  addr[nt]=&cmd.rthalo;
//  id[nt++]=DOUBLE;

	RPName(cmd.epshalo,"epshalo");
//  strcpy(tag[nt],"epshalo");
//  addr[nt]=&cmd.epshalo;
//  id[nt++]=DOUBLE;

	RPName(cmd.ahalo,"ahalo");
//  strcpy(tag[nt],"ahalo");
//  addr[nt]=&cmd.ahalo;
//  id[nt++]=DOUBLE;

	RPName(cmd.rmaxhalo,"rmaxhalo");
//  strcpy(tag[nt],"rmaxhalo");
//  addr[nt]=&cmd.rmaxhalo;
//  id[nt++]=DOUBLE;

	SPName(cmd.halotype,"halotype",100);
//  strcpy(tag[nt],"halotype"); 
//  cmd.halotype=(string) malloc(100);
//  addr[nt]=cmd.halotype;
//  id[nt++]=STRING;

// ----Input satellite's parameters--------------
 	BPName(cmd.usesat,"usesat");
//  strcpy(tag[nt],"usesat");
//  addr[nt]=&cmd.usesat;
//  id[nt++]=BOOLEAN;

 	IPName(cmd.nsat,"nsat");
//  strcpy(tag[nt],"nsat");
//  addr[nt]=&cmd.nsat;
//  id[nt++]=INT;

	RPName(cmd.satmass,"satmass");
//  strcpy(tag[nt],"satmass");
//  addr[nt]=&cmd.satmass;
//  id[nt++]=DOUBLE;

 	BPName(cmd.selfgsat,"selfgsat");
//  strcpy(tag[nt],"selfgsat");
//  addr[nt]=&cmd.selfgsat;
//  id[nt++]=BOOLEAN;

	RPName(cmd.epssat,"epssat");
//  strcpy(tag[nt],"epssat");
//  addr[nt]=&cmd.epssat;
//  id[nt++]=DOUBLE;

	RPName(cmd.asat,"asat");
//  strcpy(tag[nt],"asat");
//  addr[nt]=&cmd.asat;
//  id[nt++]=DOUBLE;

	RPName(cmd.rmaxsat,"rmaxsat");
//  strcpy(tag[nt],"rmaxsat");
//  addr[nt]=&cmd.rmaxsat;
//  id[nt++]=DOUBLE;

	RPName(cmd.xsat,"xsat");
//  strcpy(tag[nt],"xsat");
//  addr[nt]=&cmd.xsat;
//  id[nt++]=DOUBLE;

	RPName(cmd.ysat,"ysat");
//  strcpy(tag[nt],"ysat");
//  addr[nt]=&cmd.ysat;
//  id[nt++]=DOUBLE;

	RPName(cmd.zsat,"zsat");
//  strcpy(tag[nt],"zsat");
//  addr[nt]=&cmd.zsat;
//  id[nt++]=DOUBLE;

	RPName(cmd.vxsat,"vxsat");
//  strcpy(tag[nt],"vxsat");
//  addr[nt]=&cmd.vxsat;
//  id[nt++]=DOUBLE;

	RPName(cmd.vysat,"vysat");
//  strcpy(tag[nt],"vysat");
//  addr[nt]=&cmd.vysat;
//  id[nt++]=DOUBLE;

	RPName(cmd.vzsat,"vzsat");
//  strcpy(tag[nt],"vzsat");
//  addr[nt]=&cmd.vzsat;
//  id[nt++]=DOUBLE;

// ----Input addtwomods's parameters--------------
 	BPName(cmd.addmods,"addmods");
//  strcpy(tag[nt],"addmods");
//  addr[nt]=&cmd.addmods;
//  id[nt++]=BOOLEAN;

	RPName(cmd.rp,"rp");
//  strcpy(tag[nt],"rp");
//  addr[nt]=&cmd.rp;
//  id[nt++]=DOUBLE;

	RPName(cmd.rsep,"rsep");
//  strcpy(tag[nt],"rsep");
//  addr[nt]=&cmd.rsep;
//  id[nt++]=DOUBLE;

	RPName(cmd.thetmod1,"thetmod1");
//  strcpy(tag[nt],"thetmod1");
//  addr[nt]=&cmd.thetmod1;
//  id[nt++]=DOUBLE;

	RPName(cmd.phimod1,"phimod1");
//  strcpy(tag[nt],"phimod1");
//  addr[nt]=&cmd.phimod1;
//  id[nt++]=DOUBLE;

	RPName(cmd.thetmod2,"thetmod2");
//  strcpy(tag[nt],"thetmod2");
//  addr[nt]=&cmd.thetmod2;
//  id[nt++]=DOUBLE;

	RPName(cmd.phimod2,"phimod2");
//  strcpy(tag[nt],"phimod2");
//  addr[nt]=&cmd.phimod2;
//  id[nt++]=DOUBLE;

// ------------------------------------------

 	IPName(cmd.seed,"seed");
//  strcpy(tag[nt],"seed"); 
//  addr[nt]=&cmd.seed;
//  id[nt++]=INT;

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2=(char)NULL;
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
                                error("getbparam: %s=%s not bool\n",
									buf1, buf2);
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
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n",
			"%",gd.headline1,gd.headline2,"%", gd.headline3);
        fprintf(fdout,"%s\n%s\n",
		"%-------------------------------------------------------------------",
		"%");

        fprintf(fdout,FMTT,"filename",cmd.filename);
        fprintf(fdout,FMTT,"filenamefmt",cmd.filenamefmt);
        fprintf(fdout,FMTT,"statfile",cmd.statfile);
        fprintf(fdout,FMTT,"outpteps",cmd.outpteps ? "true" : "false");

        fprintf(fdout,"\n%s\n",
"%-----------------------Print bulge's parameters----------------------");
        fprintf(fdout,FMTT,"usebulge",cmd.usebulge ? "true" : "false");
        fprintf(fdout,FMTI,"nbulge",cmd.nbulge);
        fprintf(fdout,FMTR,"bulgmass",cmd.bulgmass);
        fprintf(fdout,FMTR,"abulge",cmd.abulge);
        fprintf(fdout,FMTT,"selfgbul",cmd.selfgbul ? "true" : "false");
        fprintf(fdout,FMTR,"rmaxbulg",cmd.rmaxbulg);
        fprintf(fdout,FMTR,"epsbulge",cmd.epsbulge);
        fprintf(fdout,FMTT,"axibulge",cmd.axibulge ? "true" : "false");
        fprintf(fdout,FMTR,"cbulge",cmd.cbulge);
        fprintf(fdout,FMTR,"zmaxbulg",cmd.zmaxbulg);
        fprintf(fdout,FMTI,"nsimpson",cmd.nsimpson);
        fprintf(fdout,FMTT,"bulgerot",cmd.bulgerot ? "true" : "false");
        fprintf(fdout,FMTR,"brotfrac",cmd.brotfrac);

        fprintf(fdout,"\n%s\n",
"%-----------------------Print disk's parameters----------------------");
        fprintf(fdout,FMTT,"usedisk",cmd.usedisk ? "true" : "false");
        fprintf(fdout,FMTI,"ndstars",cmd.ndstars);
        fprintf(fdout,FMTR,"z0",cmd.z0);
        fprintf(fdout,FMTT,"rsolar",cmd.rsolarstr);
        fprintf(fdout,FMTR,"qsolar",cmd.qsolar);
        fprintf(fdout,FMTR,"epsdisk",cmd.epsdisk);
        fprintf(fdout,FMTR,"zmax",cmd.zmax);
        fprintf(fdout,FMTR,"rmax",cmd.rmax);
        fprintf(fdout,FMTT,"usegas",cmd.usegas ? "true" : "false");
        fprintf(fdout,FMTI,"ndgas",cmd.ndgas);
        fprintf(fdout,FMTR,"gasmass",cmd.gasmass);
        fprintf(fdout,FMTR,"gastemp",cmd.gastemp);
        fprintf(fdout,FMTR,"z0gas",cmd.z0gas);
        fprintf(fdout,FMTR,"zmaxgas",cmd.zmaxgas);
        fprintf(fdout,FMTR,"rmaxgas",cmd.rmaxgas);
        fprintf(fdout,FMTR,"rmingas",cmd.rmingas);
        fprintf(fdout,FMTT,"selfggas",cmd.selfggas ? "true" : "false");

        fprintf(fdout,"\n%s\n",
"%-----------------------Print halo's parameters----------------------");
        fprintf(fdout,FMTT,"usehalo",cmd.usehalo ? "true" : "false");
        fprintf(fdout,FMTT,"selfghal",cmd.selfghal ? "true" : "false");
        fprintf(fdout,FMTR,"rmaxhalo",cmd.rmaxhalo);
        fprintf(fdout,FMTI,"nhalo",cmd.nhalo);
        fprintf(fdout,FMTR,"epshalo",cmd.epshalo);
        fprintf(fdout,FMTT,"halotype",cmd.halotype);
        fprintf(fdout,FMTR,"halomass",cmd.halomass);
        fprintf(fdout,FMTR,"gamhalo",cmd.gamhalo);
        fprintf(fdout,FMTR,"rthalo",cmd.rthalo);
        fprintf(fdout,FMTR,"ahalo",cmd.ahalo);

        fprintf(fdout,"\n%s\n",
"%-----------------------Print satellite's parameters----------------------");
        fprintf(fdout,FMTT,"usesat",cmd.usesat ? "true" : "false");
        fprintf(fdout,FMTR,"satmass",cmd.satmass);
        fprintf(fdout,FMTR,"asat",cmd.asat);
        fprintf(fdout,FMTR,"xsat",cmd.xsat);
        fprintf(fdout,FMTR,"ysat",cmd.ysat);
        fprintf(fdout,FMTR,"zsat",cmd.zsat);
        fprintf(fdout,FMTR,"vxsat",cmd.vxsat);
        fprintf(fdout,FMTR,"vysat",cmd.vysat);
        fprintf(fdout,FMTR,"vzsat",cmd.vzsat);
        fprintf(fdout,FMTR,"rmaxsat",cmd.rmaxsat);
        fprintf(fdout,FMTT,"selfgsat",cmd.selfgsat ? "true" : "false");
        fprintf(fdout,FMTI,"nsat",cmd.nsat);
        fprintf(fdout,FMTR,"epssat",cmd.epssat);

        fprintf(fdout,"\n%s\n",
"%-----------------------Print add two models's parameters----------------------");
        fprintf(fdout,FMTT,"addmods",cmd.addmods ? "true" : "false");
        fprintf(fdout,FMTR,"rp",cmd.rp);
        fprintf(fdout,FMTR,"rsep",cmd.rsep);
        fprintf(fdout,FMTR,"thetmod1",cmd.thetmod1);
        fprintf(fdout,FMTR,"phimod1",cmd.phimod1);
        fprintf(fdout,FMTR,"thetmod2",cmd.thetmod2);
        fprintf(fdout,FMTR,"phimod2",cmd.phimod2);

// --------------------------------------------------------------------
        fprintf(fdout,FMTI,"seed",cmd.seed); 
		fprintf(fdout,"\n\n");			// Use "\n\n" to read as input paramfile
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

local void SetUnits(void)
{
	printf("\nSet_Units...\n");

	gd.h=1.0;
	gd.G=1.0;
	gd.diskmass=1.0;
}

local void SetArrays(void)
{
	printf("\nSet arrays...\n");

	surfd=nr_dvector(1,gd.nbodies);

	x=nr_dvector(1,gd.nbodies);
	y=nr_dvector(1,gd.nbodies);
	z=nr_dvector(1,gd.nbodies);
	vx=nr_dvector(1,gd.nbodies);
	vy=nr_dvector(1,gd.nbodies);
	vz=nr_dvector(1,gd.nbodies);
	pmass=nr_dvector(1,gd.nbodies);
	radcyl=nr_dvector(1,gd.nbodies);
	radsph=nr_dvector(1,gd.nbodies);
	ax=nr_dvector(1,gd.nbodies);
	ay=nr_dvector(1,gd.nbodies);
	az=nr_dvector(1,gd.nbodies);
	aradcyl=nr_dvector(1,gd.nbodies);
	kappa=nr_dvector(1,gd.nbodies);
	sigr=nr_dvector(1,gd.nbodies);
	sigphi=nr_dvector(1,gd.nbodies);
	sigz=nr_dvector(1,gd.nbodies);
	sigt=nr_dvector(1,gd.nbodies);
	pot=nr_dvector(1,gd.nbodies);
	rotcirc=nr_dvector(1,gd.nbodies);
	rotmean=nr_dvector(1,gd.nbodies);

	rhalo=nr_dvector(1,maxtabh);
	xmhalo=nr_dvector(1,maxtabh);
	uhalo=nr_dvector(1,maxtabh);

	dadrtab=nr_dvector(1,maxtab);
}


void FreeArrays(void)
{
	printf("\nFree arrays...\n");

	free_dvector(surfd,1,gd.nbodies);

	free_dvector(x,1,gd.nbodies);
	free_dvector(y,1,gd.nbodies);
	free_dvector(z,1,gd.nbodies);
	free_dvector(vx,1,gd.nbodies);
	free_dvector(vy,1,gd.nbodies);
	free_dvector(vz,1,gd.nbodies);
	free_dvector(pmass,1,gd.nbodies);
	free_dvector(radcyl,1,gd.nbodies);
	free_dvector(radsph,1,gd.nbodies);
	free_dvector(ax,1,gd.nbodies);
	free_dvector(ay,1,gd.nbodies);
	free_dvector(az,1,gd.nbodies);
	free_dvector(aradcyl,1,gd.nbodies);
	free_dvector(kappa,1,gd.nbodies);
	free_dvector(sigr,1,gd.nbodies);
	free_dvector(sigphi,1,gd.nbodies);
	free_dvector(sigz,1,gd.nbodies);
	free_dvector(sigt,1,gd.nbodies);
	free_dvector(pot,1,gd.nbodies);
	free_dvector(rotcirc,1,gd.nbodies);
	free_dvector(rotmean,1,gd.nbodies);

	free_dvector(rhalo,1,maxtabh);
	free_dvector(xmhalo,1,maxtabh);
	free_dvector(uhalo,1,maxtabh);

	free_dvector(dadrtab,1,maxtab);
}

