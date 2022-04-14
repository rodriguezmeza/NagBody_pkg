/*==============================================================================
	MODULE: nbody_n2_io.c				[nbody_n2]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	Modules, routines and external headers:	stdinc.h, mathfns.h, vectmath,
						vectmath.h, getparam.h, types.h, stat.h
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 2007; November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.																		!
==============================================================================*/

#include "globaldefs.h"


local void outputdata(void);                    
local void diagnostics(void);

local void outputpvdata(void);                   
local void outputbindata(void);                   
local void inputdata_ascii(void);                    
local void inputdata_bin(void);                    

local void savestate(string);
local void PrintRestoreParameters(int);

local void PrintSummary(void);

local real mtot;                                
local real etot[3];                             
local matrix keten;                            
local matrix peten;                             
local vector cmpos;                             
local vector cmvel;                             
local vector amvec;                             
local real cmabs;
local real amabs;
local real sumVir;
local real rVir;

local void outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;

#define SNAP_FMT	0
#define NULL_FMT	1
#define PV_FMT		2
#define SNAP_FMT_BIN	3

void inputdata(void)
{
    if (strcmp(cmd.icfilefmt,"snap-bin")==0)
        inputdata_bin();
    else
        inputdata_ascii();
}

local void inputdata_ascii(void)
{
    stream instr;
    int ndim;
    bodyptr p;
    
    strcpy(gd.model_comment, "Input data file (snap-ascii)");
    
    instr = stropen(cmd.icfile, "r");
    in_int(instr, &cmd.nbody);
    if (cmd.nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd.nbody);
    in_int(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real(instr, &gd.tnow);
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
    
    DO_BODY(p, bodytab, bodytab+cmd.nbody)
    in_real(instr, &Mass(p));
    DO_BODY(p, bodytab, bodytab+cmd.nbody)
    in_vector(instr, Pos(p));
    DO_BODY(p, bodytab, bodytab+cmd.nbody)
    in_vector(instr, Vel(p));
    fclose(instr);
    if (scanopt(cmd.options, "reset-time"))
        gd.tnow = 0.0;
    DO_BODY(p, bodytab, bodytab+cmd.nbody)
    Type(p) = BODY;
}

local void inputdata_bin(void)
{
    stream instr;
    int ndim;
    bodyptr p;

    strcpy(gd.model_comment, "Input data file (snap-bin)");

    instr = stropen(cmd.icfile, "r");
    in_int_bin(instr, &cmd.nbody);
    if (cmd.nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd.nbody);
    in_int_bin(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    in_real_bin(instr, &gd.tnow);
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        in_real_bin(instr, &Mass(p));
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        in_vector_bin(instr, Pos(p));
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        in_vector_bin(instr, Vel(p));
    fclose(instr);
    if (scanopt(cmd.options, "reset-time"))
        gd.tnow = 0.0;
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        Type(p) = BODY;
}

#define outfile_energy	"energy.dat"
#define savestatetmp	"savestate-tmp"

void StartOutput(void)
{
    stream outstr_energy;
	char buf[200];

	if (strnull(cmd.restorefile)) {
		sprintf(buf,"%s",outfile_energy);
		if(!(outstr_energy=fopen(buf,gd.mode)))
		{
			error("StartOutput: can't open file `%s`\n",buf);
		}
		fprintf(outstr_energy,"%1s%4s%8s%8s%8s%8s%8s%8s%8s%9s%9s",
			"#","nstep","time","TotE","KinE","PotE","-K/Pot",
			"Virial","rVirial", "|Vcom|", "|Jtot|");
		fprintf(outstr_energy,"\n");
		fprintf(outstr_energy,
			"%1s%4s%8s%8s%8s%8s%8s%8s%8s%8s%8s",
			"#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>","<10>");
		fprintf(outstr_energy,"\n");
		fclose(outstr_energy);
	}

    fprintf(stdout,"  \t -- %s --\n", gd.model_comment);
	fprintf(stdout,"\n  \t running ... direct force calculation method ...\n");

    outfilefmt_string_to_int(cmd.snapoutfilefmt, &outfilefmt_int);
    fprintf(stdout,"\n%8s%8s%8s", "nbody", "dtime", "eps");
    fprintf(stdout,"%8s%8s\n", "dtout", "tstop");
    fprintf(stdout,"%8d%8.2f%8.4f", cmd.nbody, gd.dtime, cmd.eps);
    fprintf(stdout,"%8.2f%8.4f\n", gd.dtout, cmd.tstop);
    if (! strnull(cmd.options))
        fprintf(stdout,"\n\toptions: %s\n", cmd.options);

	if (! strnull(cmd.statefile)) {
		savestate(savestatetmp);
		sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
		printf("system: %s\n",buf);
		system(buf);
	}
}

void output(void)
{
    real teff;
    char   buf[200];

	teff = gd.tnow + gd.dtime/8.0;

    if (teff >= gd.toutinfo) {
		diagnostics();
		printf("\n%8s%8s%8s%8s%8s%8s%8s%10s%8s%8s\n", "time", "|K+U|", "K",
				"-U", "-K/U", "|Vcom|", "|Jtot|","nbbtot", "CPUfc", "CPUtot");
		printf("%8.3f%8.5f%8.5f%8.5f%8.5f%8.5f%8.5f%10d%8.3f%8.3f\n",
			gd.tnow, ABS(etot[0]), etot[1], -etot[2], -etot[1]/etot[2],
			cmabs, amabs, gd.nbbcalc, gd.cpuforce, cputime()-gd.cpuinit);
		PrintSummary();
		gd.toutinfo += gd.dtoutinfo;
	}

    if (! strnull(cmd.snapoutfile) && teff >= gd.tout) { 
        switch(outfilefmt_int) {
            case SNAP_FMT: 
                printf("\n\tsnap format output"); outputdata(); break;
            case NULL_FMT: 
                printf("\n\tsnap format output"); outputdata(); break;
            case PV_FMT: 
                printf("\n\tpv format output"); outputpvdata(); break;
            case SNAP_FMT_BIN: 
                printf("\n\tsnap-bin format output"); outputbindata(); break;
            default: 
                printf("\n\toutput: Unknown output format...");
                printf("\n\tprinting in default format..."); outputdata();break;
        }
        gd.tout += gd.dtout;
    }

	if (gd.nstep%cmd.stepState == 0) {
        if (! strnull(cmd.statefile)) {
            savestate(savestatetmp);
            sprintf(buf,"cp %s %s",savestatetmp,cmd.statefile);
            printf("\nsystem: %s\n",buf);
            system(buf);
        }
	}
}

local void outputdata(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf, cmd.snapoutfile, gd.nstep);
    outstr = stropen(namebuf, "w!");
    out_int(outstr, cmd.nbody);
    out_int(outstr, NDIM);
    out_real(outstr, gd.tnow);
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_real(outstr, Mass(p));
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_vector(outstr, Pos(p));
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_vector(outstr, Vel(p));
    if (scanopt(cmd.options, "out-phi"))
		DO_BODY(p, bodytab, bodytab+cmd.nbody)
			out_real(outstr, Phi(p));
    if (scanopt(cmd.options, "out-acc"))
		DO_BODY(p, bodytab, bodytab+cmd.nbody)
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
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_real_bin(outstr, Mass(p));
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_vector_bin(outstr, Pos(p));
	DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_vector_bin(outstr, Vel(p));
    if (scanopt(cmd.options, "out-phi"))
		DO_BODY(p, bodytab, bodytab+cmd.nbody)
            out_real_bin(outstr, Phi(p));
    if (scanopt(cmd.options, "out-acc"))
		DO_BODY(p, bodytab, bodytab+cmd.nbody)
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

    sprintf(namebuf, cmd.snapoutfile, gd.nstep);
    outstr = stropen(namebuf, "w!");
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        out_vector_mar(outstr, Pos(p));             
        out_vector_mar(outstr, Vel(p));
        if (scanopt(cmd.options, "out-phi"))
            out_real_mar(outstr, Phi(p));
        if (scanopt(cmd.options, "out-acc"))
            out_vector_mar(outstr, Acc(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);                             
    printf("\n\tdata output to file %s at time %f\n", namebuf, gd.tnow);
}

local void diagnostics(void)
{
    register bodyptr p;
    real velsq;
    vector tmpv;
    matrix tmpt;

    mtot = 0.0;                                 
    etot[1] = etot[2] = 0.0;                    
    CLRM(keten);                                
    CLRM(peten);                                
    CLRV(amvec);                                
    CLRV(cmpos);                                
    CLRV(cmvel);                                
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
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
        CROSSVP(tmpv, Vel(p), Pos(p));          
        MULVS(tmpv, tmpv, Mass(p));
        ADDV(amvec, amvec, tmpv);
        MULVS(tmpv, Pos(p), Mass(p));           
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           
        ADDV(cmvel, cmvel, tmpv);
    }
    etot[0] = etot[1] + etot[2];                
    DIVVS(cmpos, cmpos, mtot);                  
    DIVVS(cmvel, cmvel, mtot);

    ABSV(cmabs, cmvel);
    ABSV(amabs, amvec);
}

#define ENERGYFMT	\
"%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f"

local void PrintSummary(void)
{
    stream outstr_energy;

    outstr_energy = stropen(outfile_energy, "a");
	fprintf(outstr_energy, ENERGYFMT,
		gd.nstep, gd.tnow, etot[0], etot[1], etot[2], -etot[1]/etot[2],
           2.0*etot[1]+etot[2], mtot*mtot/rabs(etot[2]),cmabs, amabs);
	fprintf(outstr_energy,"\n");
    fclose(outstr_energy);
}

#undef ENERGYFMT
#undef outfile_energy

void savestate(string pattern)
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
    safewrite(&cmd.eps, sizeof(real), str);
    nchars = strlen(cmd.options) + 1;
    safewrite(&nchars, sizeof(int), str);
    safewrite(cmd.options, nchars * sizeof(char), str);
    safewrite(&cmd.tstop, sizeof(real), str);
    safewrite(&gd.dtout, sizeof(real), str);
    safewrite(&gd.dtoutinfo, sizeof(real), str);
    safewrite(&gd.tnow, sizeof(real), str);
    safewrite(&gd.tout, sizeof(real), str);
    safewrite(&gd.toutinfo, sizeof(real), str);
    safewrite(&gd.nstep, sizeof(int), str);
    safewrite(&cmd.stepState, sizeof(int), str);
    safewrite(&cmd.nbody, sizeof(int), str);
    safewrite(bodytab, cmd.nbody * sizeof(body), str);
    fclose(str);
}

void restorestate(string file)
{
    stream str;
    int nchars;
    string program, version;

    str = stropen(file, "r");
    saferead(&nchars, sizeof(int), str);
    program = (string) allocate(nchars * sizeof(char));
    saferead(program, nchars * sizeof(char), str);
    saferead(&nchars, sizeof(int), str);
    version = (string) allocate(nchars * sizeof(char));
    saferead(version, nchars * sizeof(char), str);
    if (! streq(program, getargv0()) ||        
          ! streq(version, getversion()))
        printf("warning: state file may be outdated\n\n");
    saferead(&gd.dtime, sizeof(real), str);
    saferead(&cmd.eps, sizeof(real), str);
    saferead(&nchars, sizeof(int), str);
    cmd.options = (string) allocate(nchars * sizeof(char));
    saferead(cmd.options, nchars * sizeof(char), str);
    saferead(&cmd.tstop, sizeof(real), str);
    saferead(&gd.dtout, sizeof(real), str);
    saferead(&gd.dtoutinfo, sizeof(real), str);
    saferead(&gd.tnow, sizeof(real), str);
    saferead(&gd.tout, sizeof(real), str);
    saferead(&gd.toutinfo, sizeof(real), str);
    saferead(&gd.nstep, sizeof(int), str);
    saferead(&cmd.stepState, sizeof(int), str);
    saferead(&cmd.nbody, sizeof(int), str);
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
    saferead(bodytab, cmd.nbody * sizeof(body), str);
    fclose(str);
}

local void PrintRestoreParameters(int mode)
{
	fprintf(gd.outlog,"\n\n[%d] Restore parameters: \n",mode);
	fprintf(gd.outlog,"dtime = %g\n",gd.dtime);
	fprintf(gd.outlog,"eps = %g\n",cmd.eps);
	fprintf(gd.outlog,"options = %s\n",cmd.options);
	fprintf(gd.outlog,"tstop = %g\n",cmd.tstop);
	fprintf(gd.outlog,"dtout = %g\n",gd.dtout);
	fprintf(gd.outlog,"dtoutinfo = %g\n",gd.dtoutinfo);
	fprintf(gd.outlog,"tnow = %g\n",gd.tnow);
	fprintf(gd.outlog,"tout = %g\n",gd.tout);
	fprintf(gd.outlog,"toutinfo = %g\n",gd.toutinfo);
	fprintf(gd.outlog,"nstep = %d\n",gd.nstep);
	fprintf(gd.outlog,"stepState = %d\n",cmd.stepState);
	fprintf(gd.outlog,"nbody = %d\n",cmd.nbody);
	fflush(gd.outlog);
}

local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"snap-ascii") == 0) *outfmt_int = 0;
    if (strnull(outfmt_str)) *outfmt_int = 1;
    if (strcmp(outfmt_str,"snap-pv") == 0) *outfmt_int = 2;
    if (strcmp(outfmt_str,"snap-bin") == 0) *outfmt_int = 3;
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
        savestate(stopsavestate);
	}
}

#undef stopfilename
#undef stopsavestate

void EndRun(void)
{
    char   buf[200];
	FILE *fd;

    if (! strnull(cmd.statefile)) {
        printf("\nsaving final run state...\n");
        savestate(savestatetmp);
        sprintf(buf,"cp %s %s", savestatetmp, cmd.statefile);
        system(buf);
    }

	fclose(gd.outlog);
	printf("\nFinal CPU time : %lf\n\n", cputime() - gd.cpuinit);
}

#undef savestatetmp

