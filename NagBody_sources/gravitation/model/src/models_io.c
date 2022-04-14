/*==============================================================================
	MODULE: models_io.c		[model]
	Written by: M.A. Rodriguez-Meza
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

	Major revisions: July 23, 2007
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
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
    fprintf(stdout,"  \t -- %s --\n", gd.model_comment);
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

void InputData(void)
{
	int ifile, tnbodies, nbodies[gd.nfiles], ndim[gd.nfiles], k;
	bodyptr btab[gd.nfiles], p, q, btabtmp;
	bool exist_snap;
	int step=0;

    strcpy(gd.model_comment, "Input data files");

	gd.headerfmt = "snap-blj-ascii";

	for (ifile=0; ifile<gd.nfiles; ifile++) {
		fprintf(gd.outlog,"\nReading file %s, snap format %s (%d/%d) ...",
			gd.filenames[ifile], gd.filenamefmts[ifile], ifile+1, gd.nfiles);
		inputdata(gd.filenames[ifile], gd.filenamefmts[ifile], gd.headerfmt,
			step,&nbodies[ifile], &ndim[ifile], &gd.tnow, &exist_snap, 
			&hdr,cmd.options, gd.model_comment);
		btab[ifile] = (bodyptr) allocate(nbodies[ifile] * sizeof(body));
		q = btab[ifile];
		DO_BODY(p, bodytab, bodytab+nbodies[ifile]) {
			Mass(q) = Mass(p);
			Type(q) = Type(p);
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

	q = bodytab;
	for (ifile=0; ifile<gd.nfiles; ifile++) {
		DO_BODY(p, btab[ifile], btab[ifile]+nbodies[ifile]) {
			Mass(q) = Mass(p);
			Type(q) = Type(p);
			for (k=0; k<NDIM; k++) {
				Pos(q)[k] = Pos(p)[k];
				Vel(q)[k] = Vel(p)[k];
			}
			q++;
		}
	}
	gd.tnow=0.;
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

#define SNAP_FMT_GADGET	4

void output(void)
{
    outfilefmt_string_to_int(cmd.outfmt, &outfilefmt_int);
    if (! strnull(cmd.out) ) {
        switch(outfilefmt_int) {
            case SNAP_FMT_GADGET:
                printf("\n\tsnap-gadget-bin format output...\n"); outputgadgetdata(); break;
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
    //    int Ids[gal->num_part[0]+gal->num_part[1]+gal->num_part[2]];
    int Ids[cmd.nbody];
    char buf[200];
#define SKIP2 fwrite(&dummy, sizeof(dummy), 1, fp1);
    
    int      NumPart, Ngas;
    double   Time, Redshift;

//    printf("\n In outputgadgetdata... \n");

    //    NumPart = gal->num_part[0]+gal->num_part[1]+gal->num_part[2];
    NumPart = cmd.nbody;
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
    
    // Set the header values to some defaults.
    //    header1.npart[1] = gal->num_part[1];
    //    header1.npart[2] = gal->num_part[0];
    //    header1.npart[3] = gal->num_part[2];
    //    header1.npartTotal[1] = gal->num_part[1];
    //    header1.npartTotal[2] = gal->num_part[0];
    //    header1.npartTotal[3] = gal->num_part[2];
    local_header1.npart[1] = cmd.nbody;
    local_header1.npart[2] = 0;
    local_header1.npart[3] = 0;
    local_header1.npartTotal[1] = cmd.nbody;
    local_header1.npartTotal[2] = 0;
    local_header1.npartTotal[3] = 0;
    //
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

//    printf("\n In outputgadgetdata... \n");

    if (!(local_P=malloc(NumPart*sizeof(local_particle_data)))) {
        fprintf(stderr,"Unable to create particle data structure in memory.");
        exit(0);
    }
    local_P--;
    
    // Transfer the particle data from the Starscream galaxy data structure to the Gadget
    // position_data structure.
    j = 0;
    //Gadget Format Halo
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        local_P[j].Pos[0] = Pos(p)[0];
        local_P[j].Pos[1] = Pos(p)[1];
        local_P[j].Pos[2] = Pos(p)[2];
        local_P[j].Vel[0] = Vel(p)[0];
        local_P[j].Vel[1] = Vel(p)[1];
        local_P[j].Vel[2] = Vel(p)[2];
        local_P[j].Mass = Mass(p);
        local_P[j].Type = 1;
        ++j;
    }
    /*        while(j < gal->num_part[0]+gal->num_part[1]+gal->num_part[2]) {
     for (i = gal->num_part[0]; i < gal->num_part[0]+gal->num_part[1]; ++i) {
     P[j].Pos[0] = gal->x[i];
     P[j].Pos[1] = gal->y[i];
     P[j].Pos[2] = gal->z[i];
     P[j].Vel[0] = gal->vel_x[i];
     P[j].Vel[1] = gal->vel_y[i];
     P[j].Vel[2] = gal->vel_z[i];
     //P[j].Mass = gal->mass[i]/unit_mass;
     P[j].Mass = gal->mass[i];
     P[j].Type = 1;
     ++j;
     } */
    //Gadget Format Disk
    /*        for (i = 0; i < gal->num_part[0]; ++i) {
     P[j].Pos[0] = gal->x[i];
     P[j].Pos[1] = gal->y[i];
     P[j].Pos[2] = gal->z[i];
     P[j].Vel[0] = gal->vel_x[i];
     P[j].Vel[1] = gal->vel_y[i];
     P[j].Vel[2] = gal->vel_z[i];
     //P[j].Mass = gal->mass[i]/unit_mass;
     P[j].Mass = gal->mass[i];
     P[j].Type = 2;
     ++j;
     }
     //Gadget Format Bulge
     for (i = gal->num_part[0]+gal->num_part[1]; i < gal->num_part[0]+gal->num_part[1]+gal->num_part[2]; ++i) {
     P[j].Pos[0] = gal->x[i];
     P[j].Pos[1] = gal->y[i];
     P[j].Pos[2] = gal->z[i];
     P[j].Vel[0] = gal->vel_x[i];
     P[j].Vel[1] = gal->vel_y[i];
     P[j].Vel[2] = gal->vel_z[i];
     //P[j].Mass = gal->mass[i]/unit_mass;
     P[j].Mass = gal->mass[i];
     P[j].Type = 3;
     ++j;
     } */
    //    }

    sprintf(fname, cmd.out, gd.nstep);
    
    //    fprintf(stderr,"Writing initial conditions... \n");
//    fprintf(stdout,"Writing initial conditions... \n");
    
    for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
        if(files>1)
            sprintf(buf,"%s.%d",fname,i);
        else
            sprintf(buf,"%s",fname);
        
        if(!(fp1=fopen(buf,"w")))
        {
            fprintf(stderr,"can't open file `%s`\n",buf);
            exit(0);
        }
        fflush(stdout);
        
        dummy = sizeof(local_header1);
        fwrite(&dummy, sizeof(dummy), 1, fp1);
        fwrite(&local_header1, sizeof(local_header1), 1, fp1);
        fwrite(&dummy, sizeof(dummy), 1, fp1);
        for(k=0, ntot_withmasses=0; k<6; k++)
        {
            if(local_header1.mass[k]==0)
                ntot_withmasses+= local_header1.npart[k];
            
        }
        dummy = 3*sizeof(float)*NumPart;
        SKIP2;
        for(k=0,pc_new=0;k<6;k++)
        {
            for(n=0;n<local_header1.npart[k];n++)
            {
                fwrite(&local_P[pc_new].Pos[0], sizeof(float), 3, fp1);
                pc_new++;
            }
        }
        SKIP2;
        SKIP2;
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<local_header1.npart[k];n++)
            {
                fwrite(&local_P[pc_new].Vel[0], sizeof(float), 3, fp1);
                pc_new++;
            }
        }
        SKIP2;
        
        dummy = sizeof(int)*NumPart;
        SKIP2;
        
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<local_header1.npart[k];n++)
            {
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
        for(k=0, pc_new=pc; k<6; k++)
        {
            for(n=0;n<local_header1.npart[k];n++)
            {
                local_P[pc_new].Type=k;
                
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
        if(local_header1.npart[0]>0)
        {
            dummy = sizeof(float)*local_header1.npart[0];
            SKIP2;
            for(n=0, pc_sph=pc; n<local_header1.npart[0];n++)
            {
                fwrite(&local_P[pc_sph].U, sizeof(float), 1, fp1);
                pc_sph++;
            }
            SKIP2;
            //Solo si se quiere escibir densidad
            /*	  SKIP2;
             printf("Writing particle density data with buffer size = %d \n",dummy);
             for(n=0, pc_sph=pc; n<header1.npart[0];n++)
             {
             fwrite(&P[pc_sph].Rho, sizeof(float), 1, fp1);
             pc_sph++;
             }
             SKIP2;
             
             if(header1.flag_cooling)
             {
             SKIP2;
             printf("Writing particle smoothing data with buffer size = %d \n",dummy);
             for(n=0, pc_sph=pc; n<header1.npart[0];n++)
             {
             fwrite(&P[pc_sph].Ne, sizeof(float), 1, fp1);
             pc_sph++;
             }
             SKIP2;
             }*/
        }
        
    }

    fprintf(stdout,"\tdata output written to file %s \n",cmd.out);
    local_P++; free(local_P);
    fclose(fp1);
}

local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
//    if (strcmp(outfmt_str,"snap-ascii") == 0) *outfmt_int = 0;
//    if (strnull(outfmt_str)) *outfmt_int = 1;
//    if (strcmp(outfmt_str,"snap-pv") == 0) *outfmt_int = 2;
//    if (strcmp(outfmt_str,"snap-bin") == 0) *outfmt_int = 3;
    if (strcmp(outfmt_str,"snap-gadget") == 0) *outfmt_int = 4;
}

#undef SNAP_FMT_GADGET

//
// FIN :::: RUTINAS PARA ESCRIBIR FORMATO GADGET (VERSION STARSCREAM)
//
//
