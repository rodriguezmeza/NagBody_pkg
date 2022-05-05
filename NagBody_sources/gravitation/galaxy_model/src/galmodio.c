/* =============================================================================
	MODULE: galmodio.c			[galaxy_models]
	Written by: M.A. Rodriguez-Meza
	Starting date:	May 2006
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: June 11, 2008;
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


local void outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;
#define SNAP_FMT			0
#define NULL_FMT			1
#define PV_FMT				2
#define SNAP_FMT_BIN		3
#define GADGET_FMT_ASCII	4
#define HERNQUIST_FMT		5
#define SNAP_BDHG_FMT		6

local void fprint_header(FILE *);
local void outputdata_galmod(void);
local void outputdata_pv(void);
local void outputdata_hernquist(void);
local void outputdata_bdhg(void);

void StartOutput(void)
{
    outfilefmt_string_to_int(cmd.filenamefmt, &outfilefmt_int);
}

void output(void)
{
        switch(outfilefmt_int) {
            case SNAP_FMT: 
                printf("\n\tsnap format output"); outputdata_galmod(); break;
            case NULL_FMT: 
                printf("\n\tsnap bdhg format output"); outputdata_bdhg(); break;
            case PV_FMT: 
                printf("\n\tpv format output"); outputdata_pv(); break;
            case HERNQUIST_FMT:
                printf("\n\thernquist format output"); 
				outputdata_hernquist(); break;
            case SNAP_BDHG_FMT:
                printf("\n\tsnap bdhg format output"); 
				outputdata_bdhg(); break;
            default:
                printf("\n\toutput: Unknown output format...");
                printf("\n\tprinting in default bdhg format..."); 
				outputdata_bdhg(); break;
        }
}

#undef SNAP_FMT
#undef NULL_FMT
#undef PV_FMT
#undef SNAP_FMT_BIN
#undef GADGET_FMT_ASCII
#undef HERNQUIST_FMT
#undef SNAP_BDHG_FMT

void EndRun(void)
{
	FreeArrays();
	fclose(gd.outlog);
	printf("\n\nTotal running cpu time: %gm\n\n",cputime()-gd.cpuinit);
}

local void outputdata_galmod(void)
{
    stream outstr;
    int p;
	real tnow;

	printf("\nWriting galaxy's data in snap format...\n");
 
    if(!(outstr=fopen(cmd.filename,"w"))) {
        fprintf(stdout,"outputdata : error opening file '%s' \n",cmd.filename);
        exit(0);
    }

	tnow=0.0;

    out_int(outstr, gd.nbodies);
    out_int(outstr, NDIM);
    out_real(outstr, tnow);
    for (p =1; p <= gd.nbodies; p++)
		fprintf(outstr,"%g\n",pmass[p]);
    for (p = 1; p <= gd.nbodies; p++)
		fprintf(outstr,"%g %g %g\n",x[p],y[p],z[p]);
    for (p = 1; p <= gd.nbodies; p++)
		fprintf(outstr,"%g %g %g\n",vx[p],vy[p],vz[p]);
    fclose(outstr);
	printf("done.\n");
}

local void outputdata_pv(void)
{
    stream outstr;
    int p;

	printf("\nWriting galaxy's data in pv format...\n");
 
    if(!(outstr=fopen(cmd.filename,"w"))) {
        fprintf(stdout,"outputdata : error opening file '%s' \n",cmd.filename);
        exit(0);
    }

    for (p = 1; p <= gd.nbodies; p++)
		fprintf(outstr,"%g %g %g %g %g %g\n",x[p],y[p],z[p],vx[p],vy[p],vz[p]);
    fclose(outstr);
	printf("done.\n");
}

local void outputdata_bdhg(void)
{
    stream outstr;
	real tnow;
	int ndim,i;
	real t,xt,yt,zt,vxt,vyt,vzt,one,pi,costhet1,sinthet1,cosphi1,
		sinphi1,costhet2,sinthet2,cosphi2,sinphi2;

	printf("\nWriting galaxy's data in snap bdhg format...\n");
 
    if(!(outstr=fopen(cmd.filename,"w"))) {
        fprintf(stdout,"outputdata : error opening file '%s' \n",cmd.filename);
        exit(0);
    }

	ndim=3;
	tnow=0.0;

	if (cmd.addmods)
		out_int(outstr, 2*gd.nbodies);
	else
		out_int(outstr, gd.nbodies);

    out_int(outstr, ndim);
    out_real(outstr, tnow);

	for (i=1; i<=cmd.ndgas; i++)
		fprintf(outstr,"%g\n",pmass[i]);
 
	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++)
			fprintf(outstr,"%g\n",pmass[i]);
 
	for (i=cmd.ndgas+1; i<=gd.nbodies; i++)
		fprintf(outstr,"%g\n",pmass[i]);
 
	if (cmd.addmods)
		for (i=cmd.ndgas+1; i<=gd.nbodies; i++)
			fprintf(outstr,"%g\n",pmass[i]);
 
	if (!cmd.addmods) {
		gd.xmod1=0.0;
		gd.ymod1=0.0;
		gd.zmod1=0.0;
		gd.vxmod1=0.0;
		gd.vymod1=0.0;
		gd.vzmod1=0.0;
		cmd.thetmod1=0.0;
		cmd.phimod1=0.0;
		costhet1=1.0;
		sinthet1=0.0;
		cosphi1=1.0;
		sinphi1=0.0;
	} else {
		cmd.thetmod1=2.0*PI*cmd.thetmod1/360.0;
		cmd.phimod1=2.0*PI*cmd.phimod1/360.0;
		cmd.thetmod2=2.0*PI*cmd.thetmod2/360.0;
		cmd.phimod2=2.0*PI*cmd.phimod2/360.0;
		costhet1=rcos(cmd.thetmod1);
		sinthet1=rsin(cmd.thetmod1);
		cosphi1=rcos(cmd.phimod1);
		sinphi1=rsin(cmd.phimod1);
		costhet2=rcos(cmd.thetmod2);
		sinthet2=rsin(cmd.thetmod2);
		cosphi2=rcos(cmd.phimod2);
		sinphi2=rsin(cmd.phimod2);
	}

	for (i=1; i<=cmd.ndgas; i++) {
		xt= x[i]*cosphi1+y[i]*sinphi1*costhet1+z[i]*sinphi1*sinthet1;
		yt= -x[i]*sinphi1+y[i]*cosphi1*costhet1+z[i]*cosphi1*sinthet1;
		zt= -y[i]*sinthet1+z[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", xt+gd.xmod1,yt+gd.ymod1,zt+gd.zmod1);
	}
 
	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++) {
			xt= x[i]*cosphi2+y[i]*sinphi2*costhet2+z[i]*sinphi2*sinthet2;
			yt= -x[i]*sinphi2+y[i]*cosphi2*costhet2+z[i]*cosphi2*sinthet2;
			zt= -y[i]*sinthet2+z[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", xt+gd.xmod2,yt+gd.ymod2,zt+gd.zmod2);
		}
 
	for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
		xt= x[i]*cosphi1+y[i]*sinphi1*costhet1+z[i]*sinphi1*sinthet1;
		yt= -x[i]*sinphi1+y[i]*cosphi1*costhet1+z[i]*cosphi1*sinthet1;
		zt= -y[i]*sinthet1+z[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", xt+gd.xmod1,yt+gd.ymod1,zt+gd.zmod1);
	}
 
	if (cmd.addmods)
		for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
			xt= x[i]*cosphi2+y[i]*sinphi2*costhet2+z[i]*sinphi2*sinthet2;
			yt= -x[i]*sinphi2+y[i]*cosphi2*costhet2+z[i]*cosphi2*sinthet2;
			zt= -y[i]*sinthet2+z[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", xt+gd.xmod2,yt+gd.ymod2,zt+gd.zmod2);
		}

	for (i=1; i<=cmd.ndgas; i++) {
		vxt= vx[i]*cosphi1+vy[i]*sinphi1*costhet1+vz[i]*sinphi1*sinthet1;
		vyt= -vx[i]*sinphi1+vy[i]*cosphi1*costhet1+vz[i]*cosphi1*sinthet1;
		vzt= -vy[i]*sinthet1+vz[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", vxt+gd.vxmod1,vyt+gd.vymod1,vzt+gd.vzmod1);
	}
 
	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++) {
			vxt= vx[i]*cosphi2+vy[i]*sinphi2*costhet2+vz[i]*sinphi2*sinthet2;
			vyt= -vx[i]*sinphi2+vy[i]*cosphi2*costhet2+vz[i]*cosphi2*sinthet2;
			vzt= -vy[i]*sinthet2+vz[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", 
					vxt+gd.vxmod2,vyt+gd.vymod2,vzt+gd.vzmod2);
		}
 
	for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
		vxt= vx[i]*cosphi1+vy[i]*sinphi1*costhet1+vz[i]*sinphi1*sinthet1;
		vyt= -vx[i]*sinphi1+vy[i]*cosphi1*costhet1+vz[i]*cosphi1*sinthet1;
		vzt= -vy[i]*sinthet1+vz[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", vxt+gd.vxmod1,vyt+gd.vymod1,vzt+gd.vzmod1);
	}
 
	if (cmd.addmods)
		for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
			vxt= vx[i]*cosphi2+vy[i]*sinphi2*costhet2+vz[i]*sinphi2*sinthet2;
			vyt= -vx[i]*sinphi2+vy[i]*cosphi2*costhet2+vz[i]*cosphi2*sinthet2;
			vzt= -vy[i]*sinthet2+vz[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", 
					vxt+gd.vxmod2,vyt+gd.vymod2,vzt+gd.vzmod2);
		}

	if (cmd.outpteps) { 

		for (i=cmd.ndgas+1; i<=gd.ndisk; i++)
			fprintf(outstr,"%g\n", cmd.epsdisk);
 
		for (i=1; i<=cmd.nbulge; i++)
			fprintf(outstr,"%g\n", cmd.epsbulge);
 
		for (i=1; i<=cmd.nhalo; i++)
			fprintf(outstr,"%g\n", cmd.epshalo);
 
		for (i=1; i<=cmd.nsat; i++)
			fprintf(outstr,"%g\n", cmd.epssat);
  
		if (cmd.addmods) {

			for (i=cmd.ndgas+1; i<=gd.ndisk; i++)
				fprintf(outstr,"%g\n", cmd.epsdisk);
 
			for (i=1; i<=cmd.nbulge; i++)
				fprintf(outstr,"%g\n", cmd.epsbulge);
 
			for (i=1; i<=cmd.nhalo; i++)
				fprintf(outstr,"%g\n", cmd.epshalo);
		}
	}

	if (cmd.usegas)
		for (i=1; i<=cmd.ndgas; i++)
			fprintf(outstr,"%g\n", cmd.gastemp);

	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++)
			fprintf(outstr,"%g\n", cmd.gastemp);

    fclose(outstr);
	printf("done.\n");
}

void outputdata_hernquist(void)
{
	stream outstr;
	int ndim,ns,i;
	real t,xt,yt,zt,vxt,vyt,vzt,one,pi,costhet1,sinthet1,cosphi1,
		sinphi1,costhet2,sinthet2,cosphi2,sinphi2;

	printf("\nOut_Bods...\n");
 
	for (i=1; i<=cmd.ndgas; i++)
		if(!cmd.selfggas) pmass[i]=gd.xgasmass/cmd.ndgas;
 
	ndim=3;
	t=0.;
	ns=0;

    if(!(outstr=fopen(cmd.filename,"w"))) {
        fprintf(stdout,"galmod_io : 001 : error opening file '%s' \n",
				cmd.filename);
        exit(0);
    }
 
	if (cmd.ndgas==0) {
		if (!cmd.addmods)
			fprintf(outstr,"%d %d %d\n",gd.nbodies,cmd.ndgas,cmd.ndgas);
		else
			fprintf(outstr,"%d %d %d\n",2*gd.nbodies,cmd.ndgas,cmd.ndgas);

		fprintf(outstr,"%d\n",ndim);
		fprintf(outstr,"%g\n",t);
	} else {
		if (!cmd.addmods)
			fprintf(outstr,"%d %d %d\n",gd.nbodies,cmd.ndgas,0);
		else
			fprintf(outstr,"%d %d %d\n",2*gd.nbodies,2*cmd.ndgas,0);

		fprintf(outstr,"%d\n",ndim);
		fprintf(outstr,"%g\n",t);
	}

	fprintf(stdout,"\n%%%%%%%%%%\n");
	if (!cmd.addmods)
		fprintf(stdout,"%d %d %g %d", gd.nbodies,ndim,t,cmd.ndgas);
	else
		fprintf(stdout,"%d %d %g %d",2*gd.nbodies,ndim,t,2*cmd.ndgas);
	fprintf(stdout,"\n%%%%%%%%%%\n");

	for (i=1; i<=cmd.ndgas; i++)
		fprintf(outstr,"%g\n",pmass[i]);
 
	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++)
			fprintf(outstr,"%g\n",pmass[i]);
 
	for (i=cmd.ndgas+1; i<=gd.nbodies; i++)
		fprintf(outstr,"%g\n",pmass[i]);
 
	if (cmd.addmods)
		for (i=cmd.ndgas+1; i<=gd.nbodies; i++)
			fprintf(outstr,"%g\n",pmass[i]);
 
	if (!cmd.addmods) {
		gd.xmod1=0.0;
		gd.ymod1=0.0;
		gd.zmod1=0.0;
		gd.vxmod1=0.0;
		gd.vymod1=0.0;
		gd.vzmod1=0.0;
		cmd.thetmod1=0.0;
		cmd.phimod1=0.0;
		costhet1=1.0;
		sinthet1=0.0;
		cosphi1=1.0;
		sinphi1=0.0;
	} else {
		cmd.thetmod1=2.0*PI*cmd.thetmod1/360.0;
		cmd.phimod1=2.0*PI*cmd.phimod1/360.0;
		cmd.thetmod2=2.0*PI*cmd.thetmod2/360.0;
		cmd.phimod2=2.0*PI*cmd.phimod2/360.0;
		costhet1=rcos(cmd.thetmod1);
		sinthet1=rsin(cmd.thetmod1);
		cosphi1=rcos(cmd.phimod1);
		sinphi1=rsin(cmd.phimod1);
		costhet2=rcos(cmd.thetmod2);
		sinthet2=rsin(cmd.thetmod2);
		cosphi2=rcos(cmd.phimod2);
		sinphi2=rsin(cmd.phimod2);
	}

	for (i=1; i<=cmd.ndgas; i++) {
		xt= x[i]*cosphi1+y[i]*sinphi1*costhet1+z[i]*sinphi1*sinthet1;
		yt= -x[i]*sinphi1+y[i]*cosphi1*costhet1+z[i]*cosphi1*sinthet1;
		zt= -y[i]*sinthet1+z[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", xt+gd.xmod1,yt+gd.ymod1,zt+gd.zmod1);
	}
 
	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++) {
			xt= x[i]*cosphi2+y[i]*sinphi2*costhet2+z[i]*sinphi2*sinthet2;
			yt= -x[i]*sinphi2+y[i]*cosphi2*costhet2+z[i]*cosphi2*sinthet2;
			zt= -y[i]*sinthet2+z[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", xt+gd.xmod2,yt+gd.ymod2,zt+gd.zmod2);
		}
 
	for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
		xt= x[i]*cosphi1+y[i]*sinphi1*costhet1+z[i]*sinphi1*sinthet1;
		yt= -x[i]*sinphi1+y[i]*cosphi1*costhet1+z[i]*cosphi1*sinthet1;
		zt= -y[i]*sinthet1+z[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", xt+gd.xmod1,yt+gd.ymod1,zt+gd.zmod1);
	}
 
	if (cmd.addmods)
		for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
			xt= x[i]*cosphi2+y[i]*sinphi2*costhet2+z[i]*sinphi2*sinthet2;
			yt= -x[i]*sinphi2+y[i]*cosphi2*costhet2+z[i]*cosphi2*sinthet2;
			zt= -y[i]*sinthet2+z[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", xt+gd.xmod2,yt+gd.ymod2,zt+gd.zmod2);
		}

	for (i=1; i<=cmd.ndgas; i++) {
		vxt= vx[i]*cosphi1+vy[i]*sinphi1*costhet1+vz[i]*sinphi1*sinthet1;
		vyt= -vx[i]*sinphi1+vy[i]*cosphi1*costhet1+vz[i]*cosphi1*sinthet1;
		vzt= -vy[i]*sinthet1+vz[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", vxt+gd.vxmod1,vyt+gd.vymod1,vzt+gd.vzmod1);
	}
 
	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++) {
			vxt= vx[i]*cosphi2+vy[i]*sinphi2*costhet2+vz[i]*sinphi2*sinthet2;
			vyt= -vx[i]*sinphi2+vy[i]*cosphi2*costhet2+vz[i]*cosphi2*sinthet2;
			vzt= -vy[i]*sinthet2+vz[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", 
					vxt+gd.vxmod2,vyt+gd.vymod2,vzt+gd.vzmod2);
		}
 
	for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
		vxt= vx[i]*cosphi1+vy[i]*sinphi1*costhet1+vz[i]*sinphi1*sinthet1;
		vyt= -vx[i]*sinphi1+vy[i]*cosphi1*costhet1+vz[i]*cosphi1*sinthet1;
		vzt= -vy[i]*sinthet1+vz[i]*costhet1;
		fprintf(outstr,"%g %g %g\n", vxt+gd.vxmod1,vyt+gd.vymod1,vzt+gd.vzmod1);
	}
 
	if (cmd.addmods)
		for (i=cmd.ndgas+1; i<=gd.nbodies; i++) {
			vxt= vx[i]*cosphi2+vy[i]*sinphi2*costhet2+vz[i]*sinphi2*sinthet2;
			vyt= -vx[i]*sinphi2+vy[i]*cosphi2*costhet2+vz[i]*cosphi2*sinthet2;
			vzt= -vy[i]*sinthet2+vz[i]*costhet2;
			fprintf(outstr,"%g %g %g\n", 
					vxt+gd.vxmod2,vyt+gd.vymod2,vzt+gd.vzmod2);
		}

	if (cmd.outpteps) { 

		for (i=cmd.ndgas+1; i<=gd.ndisk; i++)
			fprintf(outstr,"%g\n", cmd.epsdisk);
 
		for (i=1; i<=cmd.nbulge; i++)
			fprintf(outstr,"%g\n", cmd.epsbulge);
 
		for (i=1; i<=cmd.nhalo; i++)
			fprintf(outstr,"%g\n", cmd.epshalo);
 
		for (i=1; i<=cmd.nsat; i++)
			fprintf(outstr,"%g\n", cmd.epssat);
  
		if (cmd.addmods) {

			for (i=cmd.ndgas+1; i<=gd.ndisk; i++)
				fprintf(outstr,"%g\n", cmd.epsdisk);
 
			for (i=1; i<=cmd.nbulge; i++)
				fprintf(outstr,"%g\n", cmd.epsbulge);
 
			for (i=1; i<=cmd.nhalo; i++)
				fprintf(outstr,"%g\n", cmd.epshalo);
		}
	}

	if (cmd.usegas)
		for (i=1; i<=cmd.ndgas; i++)
			fprintf(outstr,"%g\n", cmd.gastemp);

	if (cmd.addmods)
		for (i=1; i<=cmd.ndgas; i++)
			fprintf(outstr,"%g\n", cmd.gastemp);
 
	fclose(outstr);
}

local void fprint_header(FILE *fd)
{
  fprintf(fd,"%d\n", header1.npart[0]);
  fprintf(fd,"%d\n", header1.npart[1]);
  fprintf(fd,"%d\n", header1.npart[2]);
  fprintf(fd,"%d\n", header1.npart[3]);
  fprintf(fd,"%d\n", header1.npart[4]);
  fprintf(fd,"%d\n", header1.npart[5]);
  fprintf(fd,"%e\n", header1.mass[0]);
  fprintf(fd,"%e\n", header1.mass[1]);
  fprintf(fd,"%e\n", header1.mass[2]);
  fprintf(fd,"%e\n", header1.mass[3]);
  fprintf(fd,"%e\n", header1.mass[4]);
  fprintf(fd,"%e\n", header1.mass[5]);
  fprintf(fd,"%e\n", header1.time);
  fprintf(fd,"%e\n", header1.redshift);
  fprintf(fd,"%d\n", header1.flag_sfr);
  fprintf(fd,"%d\n", header1.flag_feedback);
  fprintf(fd,"%d\n", header1.npartTotal[0]);
  fprintf(fd,"%d\n", header1.npartTotal[1]);
  fprintf(fd,"%d\n", header1.npartTotal[2]);
  fprintf(fd,"%d\n", header1.npartTotal[3]);
  fprintf(fd,"%d\n", header1.npartTotal[4]);
  fprintf(fd,"%d\n", header1.npartTotal[5]);
  fprintf(fd,"%d\n", header1.flag_cooling);
  fprintf(fd,"%d\n", header1.num_files);
  fprintf(fd,"%e\n", header1.BoxSize);
  fprintf(fd,"%e\n", header1.Omega0);
  fprintf(fd,"%e\n", header1.OmegaLambda);
  fprintf(fd,"%e\n", header1.HubbleParam);
}

local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"snap-ascii") == 0) *outfmt_int = 0;
    if (strnull(outfmt_str)) *outfmt_int = 1;
    if (strcmp(outfmt_str,"snap-pv") == 0) *outfmt_int = 2;
    if (strcmp(outfmt_str,"snap-bin") == 0) *outfmt_int = 3;
    if (strcmp(outfmt_str,"gadget-ascii") == 0) *outfmt_int = 4;
    if (strcmp(outfmt_str,"hernquist") == 0) *outfmt_int = 5;
    if (strcmp(outfmt_str,"snap-bdhg") == 0) *outfmt_int = 6;
}

