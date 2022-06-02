/*==============================================================================
	MODULE: models_io.c		[model]
	Written by: Mario A. Rodriguez-Meza
	Starting date:	January, 2005
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:	stdinc.h, mathfns.h, vectmath
						vectmath.h, getparam.h
						types.h, stat.h, inout.h
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 23, 2007; may 17, 2022
	Copyright: (c) 2005-2022 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

//
//  Estructuras y definiciones para escribir en formato Gadget
//

local void output_printsnap(void);
local void outputgadgetdata(void);

typedef struct io_header_1 {
    int npart[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    int npartTotal[6];
    int flag_cooling;
    int num_files;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
} local_io_header_1;

typedef struct {
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;
    
    float Rho, U, Temp, Ne;
} local_particle_data, *local_particle_data_ptr;

local local_particle_data_ptr local_P;
local local_io_header_1 local_header1;

local void outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;

//
//  FIN :::: Estructuras y definiciones para escribir en formato Gadget
//



void StartOutput(void)
{
    fprintf(gd.outlog,"  \t -- %s --\n", gd.model_comment);
    printf("\n%s: ", "nbody");
    printf("%d\n", cmd.nbody);
    if (! strnull(cmd.options))                     
        printf("options: %s\n", cmd.options);
}

void outputdata_body(string filepar,int ni,int nf)
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

int hType;
int dType;
int bType;

void InputData(void)
{
	int ifile, tnbodies, nbodies[gd.nfiles], ndim[gd.nfiles], k;
	bodyptr btab[gd.nfiles], p, q, btabtmp;
	bool exist_snap;
	int step=0;

    strcpy(gd.model_comment, "Input data files");

	gd.headerfmt = "snap-blj-ascii";

	for (ifile=0; ifile<gd.nfiles; ifile++) {
        fprintf(gd.outlog,"\n\nReading file %s, snap format %s (%d/%d) ...",
            gd.filenames[ifile], gd.filenamefmts[ifile], ifile+1, gd.nfiles);

		inputdata(gd.filenames[ifile], gd.filenamefmts[ifile], gd.headerfmt,
			step,&nbodies[ifile], &ndim[ifile], &gd.tnow, &exist_snap, 
			&hdr,cmd.options, gd.model_comment);

// bodytab is allocated in nagbody_io
// and should be the same as btab[ifile]::
// but the last one must be kept and the former freed
		btab[ifile] = (bodyptr) allocate(nbodies[ifile] * sizeof(body));
		q = btab[ifile];
		DO_BODY(p, bodytab, bodytab+nbodies[ifile]) {
			Mass(q) = Mass(p);
			Type(q) = Type(p);
            Id(q) = Id(p);
			for (k=0; k<NDIM; k++) {
				Pos(q)[k] = Pos(p)[k];
				Vel(q)[k] = Vel(p)[k];
			}
			q++;
		}
		free(bodytab);
		fprintf(gd.outlog," done\n");
	}

	tnbodies = 0;
	for (ifile=0; ifile<gd.nfiles; ifile++)
		tnbodies += nbodies[ifile];

    fprintf(gd.outlog,"\ntotal bodies : %d\n", tnbodies);

	cmd.nbody=tnbodies;
    bodytab = (bodyptr) allocate(tnbodies * sizeof(body));

    real hMass, dMass, bMass, tMass;
    vector rcmH, rcmD, rcmB;
    vector vcmH, vcmD, vcmB;
    CLRV(rcmH);
    CLRV(vcmH);
    CLRV(rcmD);
    CLRV(vcmD);
    CLRV(rcmB);
    CLRV(vcmB);
    hType=0;
    dType=0;
    bType=0;
    hMass=0.;
    dMass=0.;
    bMass=0.;
	q = bodytab;
	for (ifile=0; ifile<gd.nfiles; ifile++) {
		DO_BODY(p, btab[ifile], btab[ifile]+nbodies[ifile]) {
			Mass(q) = Mass(p);
			Type(q) = Type(p);
            if (Type(p)==1) {
                hType++;
                hMass += Mass(p);
            }
            if (Type(p)==2) {
                dType++;
                dMass += Mass(p);
            }
            if (Type(p)==3) {
                bType++;
                bMass += Mass(p);
            }
			for (k=0; k<NDIM; k++) {
				Pos(q)[k] = Pos(p)[k];
				Vel(q)[k] = Vel(p)[k];
			}
			q++;
		}
	}
    for (ifile=0; ifile<gd.nfiles; ifile++) {
        DO_BODY(p, btab[ifile], btab[ifile]+nbodies[ifile]) {
            if (Type(p)==1) {
                ADDMULVS(rcmH, Pos(p), Mass(p) / hMass);
                ADDMULVS(vcmH, Vel(p), Mass(p) / hMass);
            }
            if (Type(p)==2) {
                ADDMULVS(rcmD, Pos(p), Mass(p) / dMass);
                ADDMULVS(vcmD, Vel(p), Mass(p) / dMass);
            }
            if (Type(p)==3) {
                ADDMULVS(rcmB, Pos(p), Mass(p) / bMass);
                ADDMULVS(vcmB, Vel(p), Mass(p) / bMass);
            }
        }
    }
    if (scanopt(cmd.options, "fix-CM-HDB")) {
        fprintf(stdout,"\nFixing CM-HDB...");
        DO_BODY(q, bodytab, bodytab+cmd.nbody) {
            if (Type(q)==1) {
                SUBV(Pos(q), Pos(q), rcmH);
                SUBV(Vel(q), Vel(q), vcmH);
            }
            if (Type(q)==2) {
                SUBV(Pos(q), Pos(q), rcmD);
                SUBV(Vel(q), Vel(q), vcmD);
            }
            if (Type(q)==3) {
                SUBV(Pos(q), Pos(q), rcmB);
                SUBV(Vel(q), Vel(q), vcmB);
            }
        }
    }

    gd.tnow=0.;
    tMass = hMass+dMass+bMass;
// Fix center of mass

    vector rcm, vcm;
    CLRV(rcm);
    CLRV(vcm);
    rcm[0] = rcmH[0]*hMass/tMass + rcmD[0]*dMass/tMass + rcmB[0]*bMass/tMass;
    rcm[1] = rcmH[1]*hMass/tMass + rcmD[1]*dMass/tMass + rcmB[1]*bMass/tMass;
    rcm[2] = rcmH[2]*hMass/tMass + rcmD[2]*dMass/tMass + rcmB[2]*bMass/tMass;
    vcm[0] = vcmH[0]*hMass/tMass + vcmD[0]*dMass/tMass + vcmB[0]*bMass/tMass;
    vcm[1] = vcmH[1]*hMass/tMass + vcmD[1]*dMass/tMass + vcmB[1]*bMass/tMass;
    vcm[2] = vcmH[2]*hMass/tMass + vcmD[2]*dMass/tMass + vcmB[2]*bMass/tMass;
    if (scanopt(cmd.options, "fix-CM")) {
        fprintf(stdout,"\nFixing CM...");
        DO_BODY(q, bodytab, bodytab+cmd.nbody) {
            SUBV(Pos(q), Pos(q), rcm);
            SUBV(Vel(q), Vel(q), vcm);
        }
    }

    fprintf(stdout,"\nNumber and mass of DM type = %d %g\n",hType, hMass);
    fprintf(stdout,"Number and mass of Disk type = %d %g\n",dType, dMass);
    fprintf(stdout,"Number and mass of Bulge type = %d %g\n",bType, bMass);
    fprintf(stdout,"DM CM Pos and Vel : %g %g %g %g %g %g\n",
            rcmH[0],rcmH[1],rcmH[2], vcmH[0], vcmH[1], vcmH[2]);
    fprintf(stdout,"Disk CM Pos and Vel : %g %g %g %g %g %g\n",
        rcmD[0],rcmD[1],rcmD[2], vcmD[0], vcmD[1], vcmD[2]);
    fprintf(stdout,"Bulge CM Pos and Vel : %g %g %g %g %g %g\n",
        rcmB[0],rcmB[1],rcmB[2], vcmB[0], vcmB[1], vcmB[2]);
    fprintf(stdout,"Total CM Pos and Vel : %g %g %g %g %g %g\n",
        rcm[0],rcm[1],rcm[2], vcm[0], vcm[1], vcm[2]);
}

// Inputdata_gadget driver for SPECIAL 
// Particle data structure to manipulate I/O  N > 10^6 purpose...
void InputData_long_original(void) // esta tiene problemas ...
{
	int ifile, tnbodies, nbodies[gd.nfiles], ndim[gd.nfiles], k;
	bodyptr_long btab[gd.nfiles], p, q, btabtmp;
	bool exist_snap;
	int step=0;

    strcpy(gd.model_comment, "Input data files [io]");

	gd.headerfmt = "snap-blj-ascii";

	for (ifile=0; ifile<gd.nfiles; ifile++) {
		fprintf(gd.outlog,"\nReading file %s, snap format %s (%d/%d) ...",
			gd.filenames[ifile], gd.filenamefmts[ifile], ifile+1, gd.nfiles);
		inputdata(gd.filenames[ifile], gd.filenamefmts[ifile], gd.headerfmt,
			step,&nbodies[ifile], &ndim[ifile], &gd.tnow, &exist_snap, 
			&hdr,cmd.options, gd.model_comment);
		btab[ifile] = (bodyptr_long) allocate(nbodies[ifile] * sizeof(body));
		q = btab[ifile];
		DO_BODY(p, bodytab_long, bodytab_long+nbodies[ifile]) {
			Mass_long(q) = Mass_long(p);
			Type_long(q) = Type_long(p);
			for (k=0; k<NDIM; k++) {
				Pos_long(q)[k] = Pos_long(p)[k];
				Vel_long(q)[k] = Vel_long(p)[k];
			}
			q++;
		}
		free(bodytab_long);
		fprintf(gd.outlog," done\n");
	}

	tnbodies = 0;
	for (ifile=0; ifile<gd.nfiles; ifile++)
		tnbodies += nbodies[ifile];

	fprintf(gd.outlog,"\ntotal bodies : %d\n", tnbodies);

	cmd.nbody=tnbodies;
    bodytab_long = (bodyptr_long) allocate(tnbodies * sizeof(body_long));

	q = bodytab_long;
	for (ifile=0; ifile<gd.nfiles; ifile++) {
		DO_BODY(p, btab[ifile], btab[ifile]+nbodies[ifile]) {
			Mass_long(q) = Mass_long(p);
			Type_long(q) = Type_long(p);
			for (k=0; k<NDIM; k++) {
				Pos_long(q)[k] = Pos_long(p)[k];
				Vel_long(q)[k] = Vel_long(p)[k];
			}
			q++;
		}
	}
	gd.tnow=0.;
}

//
void InputData_long(void)
{
	int ifile, tnbodies, nbodies[gd.nfiles], ndim[gd.nfiles], k;
	bodyptr_long btab[gd.nfiles], p, q, btabtmp;
	bool exist_snap;
	int step=0;

    strcpy(gd.model_comment, "Input data files [io]");

	gd.headerfmt = "snap-blj-ascii";

ifile=0;
		fprintf(gd.outlog,"\nReading file %s, snap format %s (%d/%d) ...",
			gd.filenames[ifile], gd.filenamefmts[ifile], ifile+1, gd.nfiles);
		inputdata(gd.filenames[ifile], gd.filenamefmts[ifile], gd.headerfmt,
			step,&nbodies[ifile], &ndim[ifile], &gd.tnow, &exist_snap, 
			&hdr,cmd.options, gd.model_comment);
	cmd.nbody=nbodies[ifile];

	gd.tnow=0.;
}
//
//

#define SNAP_FMT_PRINTSNAP                  30
#define SNAP_FMT_SNAP_GADGET                31

void output(void)
{
    outfilefmt_string_to_int(cmd.outfmt, &outfilefmt_int);
    if (! strnull(cmd.out) ) {
        switch(outfilefmt_int) {
            case SNAP_FMT_SNAP_GADGET:
                printf("\n\tsnap-gadget-bin format output...\n"); outputgadgetdata(); break;
            case SNAP_FMT_PRINTSNAP:
                printf("\n\tsnap-pv-printfmt format output...\n"); output_printsnap(); break;
            default:
                outputdata(cmd.out, cmd.outfmt, gd.headerfmt, gd.nstep,
                           cmd.nbody, gd.tnow, &hdr, cmd.options);
        }
    } else
		error("\n\noutput: You should give an output file name\n\n");
}

void EndRun(void)
{
	fclose(gd.outlog);
	printf("\nFinal CPU time : %lf\n\n", cputime()-gd.cpuinit);
}

local void output_printsnap(void)
{
    stream outstr;
    bodyptr p;
    char fname[200];

    sprintf(fname, cmd.out, gd.nstep);

    outstr = stropen(fname, "w!");
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {
        out_vector_mar(outstr, Pos(p));
        out_vector_mar(outstr, Vel(p));
        out_int_mar(outstr, (p-bodytab+1));
        out_int_mar(outstr, Type(p));
        if (scanopt(cmd.options, "out-phi"))
            out_real_mar(outstr, Phi(p));
        if (scanopt(cmd.options, "out-acc"))
            out_vector_mar(outstr, Acc(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);
    printf("\tpos-vel-id-type data output to file %s at time %f\n\n",fname,gd.tnow);
}

//
// RUTINAS PARA ESCRIBIR FORMATO GADGET (VERSION STARSCREAM)
//
//
local void outputgadgetdata(void)
{
    bodyptr p;
    char fname[200];
    
    FILE *fp1, *fp2;
    int i, j, k, dummy, ntot_withmasses;
    int t,n,off,pc,pc_new,pc_sph;
    int files = 1;
    int Ids[cmd.nbody];
    char buf[200];
#define SKIP2 fwrite(&dummy, sizeof(dummy), 1, fp1);
    
    int      NumPart, Ngas;
    double   Time, Redshift;

    NumPart = cmd.nbody;
    fprintf(stdout,"\nDifference total number of bodies vs gadget number of types : %d\n",
            NumPart-hType-dType-bType);
    // Set everything to zero and overwrite it later if needed.
    local_header1.npart[0] = 0;
    local_header1.npart[1] = 0;
    local_header1.npart[2] = 0;
    local_header1.npart[3] = 0;
    local_header1.npart[4] = 0;
    local_header1.npart[5] = 0;
    local_header1.npartTotal[0] = 0;
    local_header1.npartTotal[1] = 0;
    local_header1.npartTotal[2] = 0;
    local_header1.npartTotal[3] = 0;
    local_header1.npartTotal[4] = 0;
    local_header1.npartTotal[5] = 0;
    local_header1.mass[0] = 0.0;
    local_header1.mass[1] = 0.0;
    local_header1.mass[2] = 0.0;
    local_header1.mass[3] = 0.0;
    local_header1.mass[4] = 0.0;
    local_header1.mass[5] = 0.0;

    local_header1.npart[1] = cmd.nbody;
    local_header1.npart[2] = 0;
    local_header1.npart[3] = 0;
    local_header1.npartTotal[1] = cmd.nbody;
    local_header1.npartTotal[2] = 0;
    local_header1.npartTotal[3] = 0;

// Gadget types::
    if (scanopt(cmd.options, "gadget-types")) {
    local_header1.npart[1] = hType;
    local_header1.npartTotal[1] = hType;
    local_header1.npart[2] = dType;
    local_header1.npartTotal[2] = dType;
    local_header1.npart[3] = bType;
    local_header1.npartTotal[3] = bType;
    }
//

    if (scanopt(cmd.options, "DM-type")) {
        local_header1.npart[1] = cmd.nbody;
        local_header1.npartTotal[1] = cmd.nbody;
    }
    if (scanopt(cmd.options, "Disk-type")) {
        local_header1.npart[2] = cmd.nbody;
        local_header1.npartTotal[2] = cmd.nbody;
        local_header1.npart[1] = 0;
        local_header1.npartTotal[1] = 0;
    }
    if (scanopt(cmd.options, "Bulge-type")) {
        local_header1.npart[3] = cmd.nbody;
        local_header1.npartTotal[3] = cmd.nbody;
        local_header1.npart[1] = 0;
        local_header1.npartTotal[1] = 0;
    }

    fprintf(stdout,"\nTypes: %d %d %d",hType,dType,bType);

    local_header1.time = 0.0;
    local_header1.redshift = 0.0;
    local_header1.flag_sfr = 0.0;
    local_header1.flag_feedback = 0.0;
    local_header1.flag_cooling = 0.0;
    local_header1.num_files = 1;
    local_header1.BoxSize = 0.0;
    local_header1.Omega0 = 0.0;
    local_header1.OmegaLambda = 0.0;
    local_header1.HubbleParam = 1.0;

    if (!(local_P=malloc(NumPart*sizeof(local_particle_data)))) {
        fprintf(stderr,"Unable to create particle data structure in memory.");
        exit(0);
    }
    local_P--;
    
// Transfer the particle data from the Starscream galaxy data structure to the Gadget
// position_data structure.
// Debug:: is it or 1? It is one.
//    j = 0;
    j = 1;
//Gadget Format Halo
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        local_P[j].Pos[0] = Pos(p)[0];
        local_P[j].Pos[1] = Pos(p)[1];
        local_P[j].Pos[2] = Pos(p)[2];
        local_P[j].Vel[0] = Vel(p)[0];
        local_P[j].Vel[1] = Vel(p)[1];
        local_P[j].Vel[2] = Vel(p)[2];
        local_P[j].Mass = Mass(p);
        local_P[j].Type = Type(p);
        if (scanopt(cmd.options, "DM-type"))
            local_P[j].Type = 1;
        else
            if (scanopt(cmd.options, "Disk-type"))
                local_P[j].Type = 2;
            else
                if (scanopt(cmd.options, "Bulge-type"))
                    local_P[j].Type = 3;
        ++j;
    }

// Recheck types::
    int rhType=0, rdType=0, rbType=0;
    for (j=1; j<=NumPart; j++) {
        if (local_P[j].Type==1)
            rhType++;
        if (local_P[j].Type==2)
            rdType++;
        if (local_P[j].Type==3)
            rbType++;
    }
    fprintf(stdout,"\nTypes: rhType, rdType and rbType: %d %d %d",rhType,rdType,rbType);
//

    sprintf(fname, cmd.out, gd.nstep);

    for(i=0, pc=1; i<files; i++, pc=pc_new){
        if(files>1)
            sprintf(buf,"%s.%d",fname,i);
        else
            sprintf(buf,"%s",fname);
        
        if(!(fp1=fopen(buf,"w"))){
            fprintf(stderr,"can't open file `%s`\n",buf);
            exit(0);
        }
        fflush(stdout);
        
        dummy = sizeof(local_header1);
        fwrite(&dummy, sizeof(dummy), 1, fp1);
        fwrite(&local_header1, sizeof(local_header1), 1, fp1);
        fwrite(&dummy, sizeof(dummy), 1, fp1);
        for(k=0, ntot_withmasses=0; k<6; k++){
            if(local_header1.mass[k]==0)
                ntot_withmasses+= local_header1.npart[k];
        }
        dummy = 3*sizeof(float)*NumPart;
        SKIP2;
        for(k=0,pc_new=0;k<6;k++){
            for(n=0;n<local_header1.npart[k];n++){
                fwrite(&local_P[pc_new].Pos[0], sizeof(float), 3, fp1);
                pc_new++;
            }
        }
        SKIP2;
        SKIP2;
        for(k=0,pc_new=pc;k<6;k++){
            for(n=0;n<local_header1.npart[k];n++){
                fwrite(&local_P[pc_new].Vel[0], sizeof(float), 3, fp1);
                pc_new++;
            }
        }
        SKIP2;

        dummy = sizeof(int)*NumPart;
        SKIP2;
        
        for(k=0,pc_new=pc;k<6;k++){
            for(n=0;n<local_header1.npart[k];n++){
                Ids[pc_new] = pc_new;
                fwrite(&Ids[pc_new], sizeof(int), 1, fp1);
                pc_new++;
            }
        }
        SKIP2;
        
        if(ntot_withmasses>0) {
            dummy = sizeof(float)*NumPart;
            SKIP2;
        }
        for(k=0, pc_new=pc; k<6; k++){
            for(n=0;n<local_header1.npart[k];n++){
                if(local_header1.mass[k]==0)
                    fwrite(&local_P[pc_new].Mass, sizeof(float), 1, fp1);
                else
                    local_P[pc_new].Mass= local_header1.mass[k];
                pc_new++;
            }
        }
        if(ntot_withmasses>0) {
            SKIP2;
        }
        //Solo si hay gas o SPH
        if(local_header1.npart[0]>0){
            dummy = sizeof(float)*local_header1.npart[0];
            SKIP2;
            for(n=0, pc_sph=pc; n<local_header1.npart[0];n++){
                fwrite(&local_P[pc_sph].U, sizeof(float), 1, fp1);
                pc_sph++;
            }
            SKIP2;
        }
    }

    fprintf(stdout,"\n\tdata output written to file %s \n",cmd.out);
    local_P++; free(local_P);
    fclose(fp1);
}

local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"snap-gadget") == 0)  *outfmt_int = SNAP_FMT_SNAP_GADGET;
    if (strcmp(outfmt_str,"gadget") == 0)       *outfmt_int = SNAP_FMT_SNAP_GADGET;
    if (strcmp(outfmt_str,"snap-pv-printfmt") == 0) *outfmt_int = SNAP_FMT_PRINTSNAP;
}

#undef SNAP_FMT_SNAP_GADGET
#undef SNAP_FMT_PRINTSNAP

//
// FIN :::: RUTINAS PARA ESCRIBIR FORMATO GADGET (VERSION STARSCREAM)
//
//
