/*==============================================================================
   NAME: gbsphio.c					[gbsph]
	Written by: Mario Alberto Rodriguez-Meza
	Starting date: October 11, 1999
	Purpose: Drivers for I/O
	Language: C
	Info:	M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revision: July 24, 2007; October 04, 2007;
	Copyright: (c) 2005-2008 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "switches.h"

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/mathfns.h"
//#include "../../../General_libs/math/vectmath.h"
//#include "../../../General_libs/general/getparam.h"
//#include "../../../General_libs/general/strings.h"

#include "globaldefs.h"
#include "protodefs.h"

//#include <sys/types.h>
//#include <sys/stat.h>

//#include "../../../General_libs/io/inout.h"

local void PrintSummary(FILE *);

local void outputdata(void);                    
local void outputbindata(void);                   
local void outputpvdata(void);

local void outputdata_body(string, int, int);

local void savestate(string);                 

local void diagnostics(void);
local void output_diagnostics(void);

local void output_energy_statistics(void);

local void treereport(void);
local void reporttree(string);
local void nodereport(nodeptr);


local void inputdata_ascii(void);                    
local void inputdata_bin(void);                    

local real mtot;                                
local real etot[3];                             
local matrix keten;                             
local matrix peten;                             
local vector cmpos;                             
local vector cmvel;                             

#if defined(THREEDIM)
local vector amvec;                             
#else
local real amvec;
#endif

local void outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;
#define SNAP_ASCII_FMT	0
#define NULL_FMT		1
#define SNAP_PV_FMT		2
#define SNAP_BIN_FMT	3


void inputdata(void)
{
    if (strcmp(cmd.icfilefmt,"snap-bin")==0)
        inputdata_bin();
    else
        inputdata_ascii();
}

void inputdata_ascii(void)
{
    stream instr;
    int ndim;
    bodyptr p;

    gd.model_comment = "Input data file";
    instr = stropen(cmd.icfile, "r");               
	printf("Reading body data from file %s ...\n",cmd.icfile);
    in_int(instr, &cmd.nbody);                      
    if (cmd.nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd.nbody);
    in_int(instr, &ndim);                       
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real(instr, &gd.tnow);                      
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
                                                
    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        in_real(instr, &Mass(p));               
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        in_vector(instr, Pos(p));               
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        in_vector(instr, Vel(p));               
    fclose(instr);                              
    if (scanopt(cmd.options, "reset-time"))         
        gd.tnow = 0.0;                             
    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        Type(p) = BODY;                         
}

local void inputdata_bin(void)
{
    stream instr;
    int ndim;
    bodyptr p;

    instr = stropen(cmd.icfile, "r");               
    in_int_bin(instr, &cmd.nbody);                      
    if (cmd.nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd.nbody);
    in_int_bin(instr, &ndim);                       
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real_bin(instr, &gd.tnow);                      
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
                                               
    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        in_real_bin(instr, &Mass(p));               
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        in_vector_bin(instr, Pos(p));               
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        in_vector_bin(instr, Vel(p));               
    fclose(instr);                              
    if (scanopt(cmd.options, "reset-time"))         
        gd.tnow = 0.0;                             
    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        Type(p) = BODY;                         
}


#define outfile_energy		"energy.dat"
#define savestatetmp	"savestate-tmp"


void StartOutput(void)
{
    string force_models_input;

    stream outstr_energy;
	char buf[200];

	if (strnull(cmd.restorefile)) {
		sprintf(buf,"%s",outfile_energy);
		if(!(outstr_energy=fopen(buf,gd.mode)))
		{
			error("StartOutput: can't open file `%s`\n",buf);
		}
		fprintf(outstr_energy,"%1s%4s%8s%8s%8s%8s%10s%8s%6s%10s%8s\n",		// rg=GM^2/|Pot|
          "#","nstep","time","TotE","KinE","PotE","-K/PotE","2K+Pot","rg","|Vcom|","|Jtot|");
		fprintf(outstr_energy,"%1s%4s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n",
          "#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>","<10>");
		fclose(outstr_energy);
	}

    outfilefmt_string_to_int(cmd.snapoutfilefmt, &outfilefmt_int);

    if (scanopt(cmd.forcecalc_method, "barnes"))
    {
        printf("  -- Barnes method for force calculation --\n");
    } else {
        if (scanopt(cmd.forcecalc_method, "normal")) {
            printf("  -- normal method for force calculation --\n");
        } else {
          if (scanopt(cmd.forcecalc_method, "direct")) {
              printf("  -- direct method for force calculation --\n");
          } else {
                if (scanopt(cmd.forcecalc_method, "motion_without_leader")) {
                    printf("Motion without leader -- force calculation\n");
                } else
//                    printf("  -- No force calculations will be done -- \n");
                    printf("  -- Barnes method for force calculation (default) -- \n");
            }
        }
    }

    if (scanopt(cmd.force_models, "external_potential"))
        printf("  -- External potential included --\n");
    if (scanopt(cmd.force_models, "scalar_field_potential_no_eps")
			&& scanopt(cmd.force_models, "scalar_field_potential_eps") 
			&& scanopt(cmd.force_models, "gravcalc_2g")) {
        printf("This multiple option is not valid!\n");
        exit(1);
    }
    if (scanopt(cmd.force_models, "scalar_field_potential_no_eps") 
			&& scanopt(cmd.force_models, "scalar_field_potential_eps")) {
        printf("This multiple option is not valid!\n");
        exit(1);
    }
    if (scanopt(cmd.force_models, "scalar_field_potential_no_eps") 
			&& scanopt(cmd.force_models, "gravcalc_2g")) {
        printf("This multiple option is not valid!\n");
        exit(1);
    }
    if (scanopt(cmd.force_models, "scalar_field_potential_eps")
			&& scanopt(cmd.force_models, "gravcalc_2g")) {
        printf("This multiple option is not valid!\n");
        exit(1);
    }
    if (scanopt(cmd.force_models, "gravity_plummer"))
        printf("  -- Gravity (Plummer-smoothing) included --\n");
    if (scanopt(cmd.force_models, "gravcalc_2g"))
        printf("  -- Gravitation calculation using G=2 --\n");
    if (scanopt(cmd.force_models, "scalar_field_potential_eps"))
        printf("  -- Scalar field potential included (smoothed version) --\n");
    if (scanopt(cmd.force_models, "scalar_field_potential_no_eps"))
        printf("  -- Scalar field potential included (no smoothed version) --\n");

    force_models_input=GetParam("force_models");
    if(strnull(force_models_input))
        printf("  -- Gravity (Plummer-smoothing) included by default --\n");

    printf("  -- %s --\n", gd.model_comment);

    printf("\n%8s%8s%8s", "nbody", "dtime", "eps");
#if !defined(QUICKSCAN)
    printf("%8s", "theta");
#endif
    printf("%8s%8s%8s\n", "usequad", "dtout", "tstop");
    printf("%8d%8.2f%8.4f", cmd.nbody, gd.dtime, cmd.eps);
#if !defined(QUICKSCAN)
    printf("%8.2f", cmd.theta);
#endif
    printf("%8s%8.2f%8.4f\n", cmd.usequad ? "true" : "false", gd.dtout, cmd.tstop);
    if (! strnull(cmd.options))
        printf("\n\toptions: %s\n", cmd.options);

	if (! strnull(cmd.statefile)) {
		savestate(savestatetmp);
		sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
		printf("system: %s\n",buf);
		system(buf);
	}

    if (scanopt(cmd.force_models, "external_potential"))
        printf("\nExternal potential parameters:\n eps_pot= %6.3f, sigma_pot= %6.3f\n",
            cmd.eps_pot,cmd.sigma_pot);
}


void forcereport(void)
{
//    printf("\n\t%8s%8s%8s%8s%12s%12s%10s\n",
//           "rsize", "tdepth", "ftree",
//           "actmax", "nbbtot", "nbctot", "CPUfc");
//    printf("      %10.1f%8d%8.3f%8d%12d%12d%10.3f\n",
//           gd.rsize, gd.tdepth, (cmd.nbody + gd.ncell - 1) / ((real) gd.ncell),
//           gd.actmax, gd.nbbcalc, gd.nbccalc, gd.cpuforce);
}


void output(void)
{
    real cmabs, amabs, teff;
    char   buf[200];
    double cpustart, cpudt;
	int flaginfo=0;

    cpustart = cputime();                       
	teff = gd.tnow + gd.dtime/8.0;

    if (teff >= gd.toutinfo) {
		flaginfo=1;
		diagnostics();                              
		ABSV(cmabs, cmvel);                         
#if defined(THREEDIM)
		ABSV(amabs, amvec);                         
#else
		amabs=rabs(amvec);
#endif
		printf("\n\t%8s%8s%8s%8s%12s%12s%10s\n",
			"rsize", "tdepth", "ftree",
			"actmax", "nbbtot", "nbctot", "CPUfc");
		printf("      %10.1f%8d%8.3f%8d%12d%12d%10.3f\n",
			gd.rsize, gd.tdepth, (cmd.nbody + gd.ncell - 1) / ((real) gd.ncell),
			gd.actmax, gd.nbbcalc, gd.nbccalc, gd.cpuforce);

		printf("\n%14s%14s%14s%12s%16s\n",
			"time", "|K+U|", "K", "-U", "-K/U");
		printf("%14.3f%14.7f%14.7f%12.7f%16.7f\n",
			gd.tnow, ABS(etot[0]), etot[1], -etot[2], -etot[1]/etot[2]);
		printf("\n\t\t\t\t%16s%16s%10s\n",
			"|Vcom|", "|Jtot|", "CPUtot");
		printf("\t\t\t\t%16.7f%16.7f%10.3f\n",
			cmabs, amabs, cputime()-gd.cpuinit);

		PrintSummary(stdout);
		gd.toutinfo += gd.dtoutinfo;
	}

//	teff = gd.tnow + gd.dtime/8.0;  
    if (teff >= gd.tout) {
//		PrintSummary(stdout);
//        gd.tout += gd.dtout;
//	}

//    if (! strnull(cmd.snapoutfile) && teff >= gd.tout) { 
    if (! strnull(cmd.snapoutfile) ) { 
        switch(outfilefmt_int) {
            case SNAP_ASCII_FMT: 
                printf("\n\tsnap-ascii format output"); outputdata(); break;
            case NULL_FMT: 
                printf("\n\tsnap-ascii format output"); outputdata(); break;
            case SNAP_PV_FMT: 
                printf("\n\tsnap-pv format output"); outputpvdata(); break;
            case SNAP_BIN_FMT: 
                printf("\n\tsnap-bin format output"); outputbindata(); break;
            default: 
                printf("\n\toutput: Unknown output format...");
                printf("\n\tprinting in default format..."); outputdata(); break;
        }
//        gd.tout += gd.dtout;

//        if (! strnull(cmd.statefile)) {
//            savestate(savestatetmp);
//            sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
//            printf("system: %s",buf);
//            system(buf);
//        }
    }
        gd.tout += gd.dtout;
	}

	if (gd.nstep%cmd.stepState == 0) {
        if (! strnull(cmd.statefile)) {
			cpudt = cputime()-gd.cpuinit;
			gd.cputotal += cpudt;
//			Global_to_Header();
            savestate(savestatetmp);
			gd.cputotal -= cpudt;
            sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
            printf("system: %s\n",buf);
            system(buf);
        }
	}

	gd.cputotout += cputime()-cpustart;
	if (flaginfo)
		printf("Individual and Accum. output cpu time : %g %g\n",
			cputime()-cpustart, gd.cputotout);
}

local void PrintSummary(FILE *fp)
{
    stream outstr_energy;
    real cmabs, amabs, mtot;
	bodyptr p;

    ABSV(cmabs, cmvel);                         
#if defined(THREEDIM)
    ABSV(amabs, amvec);                         
#else
	amabs=rabs(amvec);
#endif
	mtot=0.;
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        mtot += Mass(p);                        

    outstr_energy = stropen(outfile_energy, "a");
	fprintf(outstr_energy,
		"%5d %8.4f %g %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f",
		gd.nstep, gd.tnow, etot[0], etot[1], etot[2], -etot[1]/etot[2],
		2.0*etot[1]+etot[2], rsqr(mtot)/rabs(etot[2]), cmabs, amabs);
	fprintf(outstr_energy,"\n");
    fclose(outstr_energy);
}

#undef outfile_energy

void outputdata(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;

	real tmpv[3];
	int k;

    sprintf(namebuf, cmd.snapoutfile, gd.nstep);           
    outstr = stropen(namebuf, "w!");
    out_int(outstr, cmd.nbody);                     
    out_int(outstr, NDIM);                      
    out_real(outstr, gd.tnow);                     
    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        out_real(outstr, Mass(p));              
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        out_vector(outstr, Pos(p));
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        out_vector(outstr, Vel(p));
    if (scanopt(cmd.options, "out-phi"))            
        for (p = bodytab; p < bodytab+cmd.nbody; p++)
            out_real(outstr, Phi(p));           
    if (scanopt(cmd.options, "out-acc"))            
        for (p = bodytab; p < bodytab+cmd.nbody; p++)
            out_vector(outstr, Acc(p));  
    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, gd.tnow);
}


local void outputbindata(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf, cmd.snapoutfile, gd.nstep);
    outstr = stropen(namebuf, "w!");         
    out_int_bin(outstr, cmd.nbody);
    out_int_bin(outstr, NDIM);                      
    out_real_bin(outstr, gd.tnow);                     
    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        out_real_bin(outstr, Mass(p));              
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        out_vector_bin(outstr, Pos(p));             
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        out_vector_bin(outstr, Vel(p));             
    if (scanopt(cmd.options, "out-phi"))           
        for (p = bodytab; p < bodytab+cmd.nbody; p++)
            out_real_bin(outstr, Phi(p));           
    if (scanopt(cmd.options, "out-acc"))
        for (p = bodytab; p < bodytab+cmd.nbody; p++)
            out_vector_bin(outstr, Acc(p));
    fclose(outstr);
    printf("\n\tdata output to file %s at time %f\n", namebuf, gd.tnow);
}

local void outputpvdata(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;

	real tmpv[3];
	int k;

    sprintf(namebuf, cmd.snapoutfile, gd.nstep);           
    outstr = stropen(namebuf, "w!");         
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {
		out_int_mar(outstr, (p-bodytab+1));

        out_vector_mar(outstr, Pos(p));
        out_vector_mar(outstr, Vel(p));

        if (scanopt(cmd.options, "out-phi"))
            out_real_mar(outstr, Phi(p));
        if (scanopt(cmd.options, "out-acc"))
            out_vector_mar(outstr, Acc(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);                         
    printf("\tpos-vel data output to file %s at time %f\n\n",namebuf,gd.tnow);
}


#define report "rep%03d"
#define rprt   "rpt%03d"
#define noderprt "ndr%03d"

local void treereport(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;
    int nodcount;

    sprintf(namebuf, report, gd.nstep);           
    if (stat(namebuf, &buf) != 0)              
        outstr = stropen(namebuf, "w");        
    else                                       
        outstr = stropen(namebuf, "a");        
    nodcount = 0;
    out_int(outstr, gd.nstep_grav);                     
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {
        nodcount++;
        out_int_mar(outstr, nodcount);         
        out_real_mar(outstr, Mass(p));         
        out_vector_mar(outstr, Pos(p));        
        out_vector_mar(outstr, Vel(p));
        out_vector_mar(outstr, Acc(p));        
        out_real(outstr, Phi(p));
    }
    fclose(outstr); 
    printf("\n\tdata output to file %s at time %f\n", namebuf, gd.tnow);
}


local void reporttree(string comment)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
	nodeptr p;
    int nodcount;

    sprintf(namebuf, rprt, gd.nstep);           
    if (stat(namebuf, &buf) != 0)            
        outstr = stropen(namebuf, "w");      
    else                                     
        outstr = stropen(namebuf, "a");      

	fprintf(outstr, "\n%s\n","**************************************");
    fprintf(outstr, "\n%s%s\n", "begin:",comment); 
    fprintf(outstr,"\n%s%8d\n","nstep_grav=",gd.nstep_grav);
	fprintf(outstr,"%s%8d%8s%8d%8s%8d\n",
		"nbody=",cmd.nbody,"ncell=",gd.ncell,"tdepth=",gd.tdepth);
    fprintf(outstr, "\n%4s%8s%8s%8s%10s%10s\n",
           "nidx","ntyp","updt","mass", "pos", "vel");
    nodcount = 0;

	p = (nodeptr) gd.root;
	while (p != NULL) {
		if ( Type(p) == BODY ) {
			nodcount++;
			out_int_mar(outstr, nodcount);              
			out_int_mar(outstr, Type(p));
			out_bool_mar(outstr,Update(p));
			out_real_mar(outstr, Mass(p));              
			out_vector_mar(outstr, Pos(p));             
			out_vector(outstr, Vel(p));			
			p = Next(p);
		} else {
			nodcount++;
			out_int_mar(outstr, nodcount);              
			out_int_mar(outstr, Type(p));
			out_bool_mar(outstr,Update(p));
			out_real_mar(outstr, Mass(p));              
			out_vector_mar(outstr, Pos(p));             
			out_real(outstr, Rcrit2(p));           
			p = More(p);
		}
	} 

    fprintf(outstr, "\n%s%s\n", "end:",comment);                 
	fprintf(outstr, "\n%s\n","**************************************");
    fclose(outstr);                             
}


local void nodereport(nodeptr p)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;

    sprintf(namebuf, noderprt, gd.nstep);        
    if (stat(namebuf, &buf) != 0)             
        outstr = stropen(namebuf, "w");       
    else                                      
        outstr = stropen(namebuf, "a");       
    out_int_mar(outstr, gd.nstep_grav);
    out_int_mar(outstr, Type(p));
    out_real_mar(outstr, Mass(p));            
    out_vector_mar(outstr, Pos(p));           
    if (Type(p) == BODY) {
        out_vector_mar(outstr, Vel(p));
        out_real(outstr, Phi(p));
    }else {
        out_real(outstr, Rcrit2(p));          
    }
    fclose(outstr);                           
}


local void diagnostics(void)
{
    register bodyptr p;
    real velsq;
    vector tmpv;
    matrix tmpt;
#if !defined(THREEDIM)
	real tmps;
#endif

    mtot = 0.0;                                 
    etot[1] = etot[2] = 0.0;                    
    CLRM(keten);                                
    CLRM(peten);                                
#if defined(THREEDIM)
    CLRV(amvec);                                
#else
	amvec=0.0;
#endif    
	CLRV(cmpos);                                
    CLRV(cmvel);                                
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        mtot += Mass(p);                        
        DOTVP(velsq, Vel(p), Vel(p));           
        etot[1] += 0.5 * Mass(p) * velsq;       
        etot[2] += 0.5 * Mass(p) * Phi(p);      
        MULVS(tmpv, Vel(p), 0.5 * Mass(p));     
        OUTVP(tmpt, tmpv, Vel(p));
        ADDM(keten, keten, tmpt);
        MULVS(tmpv, Pos(p), Mass(p));           
        OUTVP(tmpt, tmpv, Acc(p));
        ADDM(peten, peten, tmpt);
#if defined(THREEDIM)
        CROSSVP(tmpv, Vel(p), Pos(p));          
        MULVS(tmpv, tmpv, Mass(p));
        ADDV(amvec, amvec, tmpv);
#endif
#if defined(TWODIM)
        CROSSVP(tmps, Vel(p), Pos(p));          
        amvec += tmps * Mass(p);
#endif
#if defined(ONEDIM)
		amvec = 0.0;
#endif
        MULVS(tmpv, Pos(p), Mass(p));           
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           
        ADDV(cmvel, cmvel, tmpv);
    }
    etot[0] = etot[1] + etot[2];                
    DIVVS(cmpos, cmpos, mtot);                  
    DIVVS(cmvel, cmvel, mtot);

}

local void savestate(string pattern)
{
    char namebuf[256];
    stream str;
    int nchars;

    sprintf(namebuf, pattern, gd.nstep & 1);       
    str = stropen(namebuf, "w!");

    nchars = strlen(getargv0()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getargv0(), nchars * sizeof(char), str);

    nchars = strlen(getversion()) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(getversion(), nchars * sizeof(char), str);

    safewrite(&gd.dtime, sizeof(real), str);
#if !defined(QUICKSCAN)
    safewrite(&cmd.theta, sizeof(real), str);
#endif
    safewrite(&cmd.usequad, sizeof(bool), str);
    safewrite(&cmd.eps, sizeof(real), str);

    nchars = strlen(cmd.options) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.options, nchars * sizeof(char), str);

    safewrite(&cmd.tstop, sizeof(real), str);
    safewrite(&gd.dtout, sizeof(real), str);
    safewrite(&gd.tnow, sizeof(real), str);
    safewrite(&gd.tout, sizeof(real), str);
    safewrite(&gd.nstep, sizeof(int), str);
    safewrite(&gd.rsize, sizeof(real), str);
    safewrite(&cmd.nbody, sizeof(int), str);
    safewrite(bodytab, cmd.nbody * sizeof(body), str);

/* 
 * Items added in this new version
 */
    nchars = strlen(cmd.forcecalc_method) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.forcecalc_method, nchars * sizeof(char), str);

    nchars = strlen(cmd.force_models) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.force_models, nchars * sizeof(char), str);

    nchars = strlen(cmd.icfile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.icfile, nchars * sizeof(char), str);

    nchars = strlen(cmd.icfilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.icfilefmt, nchars * sizeof(char), str);

    nchars = strlen(cmd.snapoutfile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.snapoutfile, nchars * sizeof(char), str);

    nchars = strlen(cmd.snapoutfilefmt) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.snapoutfilefmt, nchars * sizeof(char), str);

    nchars = strlen(cmd.statefile) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.statefile, nchars * sizeof(char), str);

    safewrite(&cmd.dm_lambda, sizeof(real), str);
    safewrite(&cmd.dm_alpha, sizeof(real), str);
//    safewrite(&cmd.dm_inv_avgphi, sizeof(real), str);
    safewrite(&cmd.G, sizeof(real), str);
    safewrite(&cmd.eps_pot, sizeof(real), str);
    safewrite(&cmd.sigma_pot, sizeof(real), str);
    safewrite(&cmd.x_pot, sizeof(real), str);
    safewrite(&cmd.y_pot, sizeof(real), str);
    safewrite(&cmd.z_pot, sizeof(real), str);
    
    
    fclose(str);
}


void restorestate(string file)
{
    stream str;
    int nchars;
    string program, version;

    str = stropen(file, "r");

    saferead(&nchars, sizeof(int), str);				// Debe coincidir con la trayectoria
    program = (string) allocate(nchars * sizeof(char));	// con la cual fue corrida la simulacion
    saferead(program, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    version = (string) allocate(nchars * sizeof(char));
    saferead(version, nchars * sizeof(char), str);

    if (! streq(program, getargv0()) ||         
          ! streq(version, getversion()))
        printf("warning: state file may be outdated\n\n");

    saferead(&gd.dtime, sizeof(real), str);
#if !defined(QUICKSCAN)
    saferead(&cmd.theta, sizeof(real), str);
#endif
    saferead(&cmd.usequad, sizeof(bool), str);
    saferead(&cmd.eps, sizeof(real), str);

    saferead(&nchars, sizeof(int), str);
    cmd.options = (string) allocate(nchars * sizeof(char));
    saferead(cmd.options, nchars * sizeof(char), str);

    saferead(&cmd.tstop, sizeof(real), str);
    saferead(&gd.dtout, sizeof(real), str);
    saferead(&gd.tnow, sizeof(real), str);
    saferead(&gd.tout, sizeof(real), str);
    saferead(&gd.nstep, sizeof(int), str);
    saferead(&gd.rsize, sizeof(real), str);
    saferead(&cmd.nbody, sizeof(int), str);
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
    saferead(bodytab, cmd.nbody * sizeof(body), str);

/* 
 * Items added in this new version
 */
    saferead(&nchars, sizeof(int), str);
    cmd.forcecalc_method = (string) allocate(nchars * sizeof(char));
    saferead(cmd.forcecalc_method, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd.force_models = (string) allocate(nchars * sizeof(char));
    saferead(cmd.force_models, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd.icfile = (string) allocate(nchars * sizeof(char));
    saferead(cmd.icfile, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd.icfilefmt = (string) allocate(nchars * sizeof(char));
    saferead(cmd.icfilefmt, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd.snapoutfile = (string) allocate(nchars * sizeof(char));
    saferead(cmd.snapoutfile, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd.snapoutfilefmt = (string) allocate(nchars * sizeof(char));
    saferead(cmd.snapoutfilefmt, nchars * sizeof(char), str);

    saferead(&nchars, sizeof(int), str);
    cmd.statefile = (string) allocate(nchars * sizeof(char));
    saferead(cmd.statefile, nchars * sizeof(char), str);

    saferead(&cmd.dm_lambda, sizeof(real), str);
    saferead(&cmd.dm_alpha, sizeof(real), str);
//    saferead(&cmd.dm_inv_avgphi, sizeof(real), str);
    saferead(&cmd.G, sizeof(real), str);
    saferead(&cmd.eps_pot, sizeof(real), str);
    saferead(&cmd.sigma_pot, sizeof(real), str);
    saferead(&cmd.x_pot, sizeof(real), str);
    saferead(&cmd.y_pot, sizeof(real), str);
    saferead(&cmd.z_pot, sizeof(real), str);
    
    fclose(str);
}

#define out_diagnostics		"diag_out"
#define FMT  "%14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n"

local void output_diagnostics(void)
{
    real cmabs, amabs, teff;
	bodyptr p,q;

    char namebuf[256];
    struct stat buf;
    stream outstr;

    sprintf(namebuf, out_diagnostics, gd.nstep);        
    if (stat(namebuf, &buf) != 0)   {            
        outstr = stropen(namebuf, "w");         
    fprintf(outstr,"\n    %8s%13s%12s%8s%8s%8s%8s%8s\n",
           "time", "     |T+U|", "T", "-U", "-T/U", "|Vcom|", 
	"|Jtot|", "CPUtot");
    } else                                        
        outstr = stropen(namebuf, "a");         

    diagnostics();                              
    ABSV(cmabs, cmvel);                         
#if defined(THREEDIM)
    ABSV(amabs, amvec);                         
#else
	amabs=rabs(amvec);
#endif

	p = bodytab;
	q = bodytab+1;
    fprintf(outstr, FMT,
    gd.tnow, ABS(etot[0]), etot[1], -etot[2], -etot[1]/etot[2],
    cmabs, amabs, amvec[2], 
	Pos(p)[0], Pos(p)[1], Pos(p)[2], Pos(q)[0], Pos(q)[1], Pos(q)[2] );
}
#undef FMT


#define out_energy_diagnostics		"energy_diagnostics"
#define FMT  "%14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n"

local void output_energy_statistics(void)
{
    real cmabs, amabs, teff;
	bodyptr p,q;

    char namebuf[256];
    struct stat buf;
    stream outstr;

    sprintf(namebuf, out_energy_diagnostics, gd.nstep);        
    if (stat(namebuf, &buf) != 0)   {            
        outstr = stropen(namebuf, "w");         
    fprintf(outstr,"\n    %8s%13s%12s%16s%16s%16s%16s\n",
           "time", "     |T+U|", "T", "-U", "-T/U", "|Vcom|", 
	"|Jtot|");
    } else                                        
        outstr = stropen(namebuf, "a");         

    diagnostics();                              
    ABSV(cmabs, cmvel);                         
#if defined(THREEDIM)
    ABSV(amabs, amvec);                         
#else
	amabs=rabs(amvec);
#endif

	p = bodytab;
	q = bodytab+1;
    fprintf(outstr, FMT,
    gd.tnow, ABS(etot[0]), etot[1], -etot[2], -etot[1]/etot[2],
    cmabs, amabs );
}

#undef FMT
#undef out_energy_diagnostics

local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"snap-ascii") == 0) *outfmt_int = SNAP_ASCII_FMT;
    if (strnull(outfmt_str)) *outfmt_int = NULL_FMT;
    if (strcmp(outfmt_str,"snap-pv") == 0) *outfmt_int = SNAP_PV_FMT;
    if (strcmp(outfmt_str,"snap-bin") == 0) *outfmt_int = SNAP_BIN_FMT;
}

#undef SNAP_ASCII_FMT
#undef NULL_FMT
#undef SNAP_PV_FMT
#undef SNAP_BIN_FMT

local void outputdata_body(string filepar,int ni,int nf)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
	bodyptr p;

    sprintf(namebuf, filepar, gd.nstep);           
    if (stat(namebuf, &buf) != 0)               
        outstr = stropen(namebuf, "w");         
    else                                        
        outstr = stropen(namebuf, "a");
    for (p = bodytab+ni; p < bodytab+nf+1; p++) {
		out_int_mar(outstr, (p-bodytab+1));
        out_vector_mar(outstr, Pos(p));
        out_vector(outstr, Vel(p));
	}
    fclose(outstr);
}

#define stopfilename	"stop"
#define stopsavestate	"stop-state"

void checkstop(void)
{
    char   buf[200];
	FILE *fd;
	double cpudt;

	if ((fd=fopen(stopfilename,"r"))) {
		fclose(fd);
		gd.stopflag = 1;
		sprintf(buf,"rm -f %s", stopfilename);
		system(buf);
        printf("\nsaving a stop run state...\n\n");
		cpudt = cputime()-gd.cpuinit;
//		gd.cputotal += cpudt;
//		Global_to_Header();
        savestate(stopsavestate);
//		gd.cputotal -= cpudt;
	}
}

#undef stopfilename
#undef stopsavestate

void EndRun(void)
{
    char   buf[200];
	FILE *fd;
	double cpudt;

    if (! strnull(cmd.statefile)) {
        printf("\nsaving final run state...\n\n");
		cpudt = cputime()-gd.cpuinit;
//		gd.cputotal += cpudt;
//		Global_to_Header();
        savestate(savestatetmp);
//		gd.cputotal -= cpudt;
        sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
        system(buf);
    }

	fclose(gd.outlog);
//	gd.cputotal += cputime()-gd.cpuinit;
//	printf("\nFinal CPU time : %lf\n\n",gd.cputotal);
	printf("\nFinal CPU time : %lf\n\n",cputime()-gd.cpuinit);
}

#undef savestatetmp

