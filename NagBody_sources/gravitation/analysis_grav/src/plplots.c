/*==============================================================================
	MODULE: plplots.c				[analysis_grav]
	Written by: Mario A. Rodriguez-Meza
	Starting date:	May 2006
	Purpose:
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

local void animation_xy(string, int);
local void animation_xy_type(string, int);
local void animation_xy_trajectory(string, int);
local void animation_rhotheta(string, int);
local void animation_rho(string, int);
//local void animation_vc(string, int);
local void animation_vcr(string, int);
local void animation_general_fx(string, int);

local void setlegends(int, char **, int, bool *, bool *, bool *, int *, 
	real *, int *, int *, int *, int *, int *);

local void copyrighttext(bool);

local int defpenwidth=1;

void snap_anim_plplot_trajectory(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;
	double cpustart;

	cpustart = cputime();

	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);

	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);

	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

    plinit();

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
			if (exist_snap) {
				++snapcount;
				PrintSnap_Spherical(snapname, snapnametmp, gd.nbody, cmd.options);
				animation_xy_trajectory(snapname, snapcount);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval);
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
			++step;
	}

    plend();
}

// AUN HAY PROBLEMAS EN EL CONSUMO DE MEMORIA...
void snap_anim_plplot(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	double cpustart;

	cpustart = cputime();

	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (scanopt(cmd.options, "save")) {
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					PrintSnap_Spherical(snapname, snapnametmp, gd.nbody, cmd.options);
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					plsetopt("-dev","psc");
					plsfnam(namebuf);
					plsetopt("-ori",cmd.pl_ori);
					if (!strnull(cmd.pl_a))
						plsetopt("-a",cmd.pl_a);
					plinit();
					if (scanopt(cmd.options, "type"))
						animation_xy_type(snapname, snapcount);
					else
						animation_xy(snapname, snapcount);
					plend();
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					PrintSnap_Spherical(snapname, snapnametmp, gd.nbody, cmd.options);
					if (scanopt(cmd.options, "type"))
						animation_xy_type(snapname, snapcount);
					else
						animation_xy(snapname, snapcount);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}

void slab_anim_plplot(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	double cpustart;

	cpustart = cputime();

	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (scanopt(cmd.options, "save")) {
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					PrintSnap_Slab(snapname, snapnametmp, gd.nbody, cmd.options);
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					plsetopt("-dev","psc");
					plsfnam(namebuf);
					plsetopt("-ori",cmd.pl_ori);
					if (!strnull(cmd.pl_a))
						plsetopt("-a",cmd.pl_a);
					plinit();
					if (scanopt(cmd.options, "type"))
						animation_xy_type(snapname, snapcount);
					else
						animation_xy(snapname, snapcount);
					plend();
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					PrintSnap_Slab(snapname, snapnametmp, gd.nbody, cmd.options);
					if (scanopt(cmd.options, "type"))
						animation_xy_type(snapname, snapcount);
					else
						animation_xy(snapname, snapcount);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}

void bodies_anim_plplot(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	double cpustart;
	char mode[4];

	cpustart = cputime();

//printf("\nAqui voy (1)\n");

	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (scanopt(cmd.options, "save")) {
		strcpy(mode,"w!");
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					if (!strnull(cmd.bodiesID))
						PrintBodies(snapname, snapnametmp, gd.nbody, 
							gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
					else
						PrintBodiesSets_Spherical(snapname, snapnametmp, gd.nbody,
							gd.nbodiesSets, &gd.bodyIDMin[0], &gd.bodyIDMax[0], 
							 cmd.RMax, cmd.options);
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					plsetopt("-dev","psc");
					plsfnam(namebuf);
					plsetopt("-ori",cmd.pl_ori);
					if (!strnull(cmd.pl_a))
						plsetopt("-a",cmd.pl_a);
					plinit();
					if (scanopt(cmd.options, "type"))
						animation_xy_type(snapname, snapcount);
					else
						animation_xy(snapname, snapcount);
					plend();
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
//printf("\nAqui voy (2)\n");
//printf("\nAqui voy (3)\n");
		strcpy(mode,"w!");
//printf("\nAqui voy (3)\n");
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
//printf("\nAqui voy (3)\n");
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow,
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
//printf("\nAqui voy (3)\n");
				if (exist_snap) {
					++snapcount;
					if (!strnull(cmd.bodiesID)) {
//printf("\nAqui voy (2)\n");
						PrintBodies(snapname, snapnametmp, gd.nbody,
                                    gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
//                                    gd.nbodiesID, &gd.bodyID[0], cmd.options, "w!");
					} else
						PrintBodiesSets_Spherical(snapname, snapnametmp, gd.nbody,
							gd.nbodiesSets, &gd.bodyIDMin[0], &gd.bodyIDMax[0], 
							 cmd.RMax, cmd.options);
					if (scanopt(cmd.options, "type"))
						animation_xy_type(snapname, snapcount);
					else
						animation_xy(snapname, snapcount);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
//printf("\nAqui voy (3)\n");
}

void bodies_anim_plplot_trajectory(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	double cpustart;
	char mode[4];

	cpustart = cputime();

	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (scanopt(cmd.options, "save")) {
		strcpy(mode,"w!");
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					plsetopt("-dev","psc");
					plsfnam(namebuf);
					plsetopt("-ori",cmd.pl_ori);
					if (!strnull(cmd.pl_a))
						plsetopt("-a",cmd.pl_a);
					plinit();
					animation_xy_trajectory(snapname, snapcount);
					plend();
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		strcpy(mode,"w!");
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
						animation_xy_trajectory(snapname, snapcount);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}

void vcr_anim_plplot(string fname, string fnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;

	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (scanopt(cmd.options, "save")) {
		gd.VcR_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					EvalVcR(fname, fnametmp);
					if (gd.VcR_flag) {
						sprintf(namebuf1, cmd.out, snapcount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
						plinit();
						animation_vcr(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
						plend();
					}
					free_memory();
					free(bodytab);
					gd.VcR_flag=0;
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		gd.VcR_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					EvalVcR(fname, fnametmp);
					if (gd.VcR_flag) {
						animation_vcr(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
					}
					free_memory();
					free(bodytab);
					gd.VcR_flag=0;
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}

void rhotheta_anim_plplot(string fname, string fnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	double cpustart;

	cpustart = cputime();

	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	plssub(1, 4);

	if (scanopt(cmd.options, "save")) {
		gd.RhoTheta_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					EvalRhoTheta(fname, fnametmp);
					if (gd.RhoTheta_flag) {
						sprintf(namebuf1, cmd.out, snapcount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
						plssub(1, 4);
						plinit();
						animation_rhotheta(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
						plend();
					}
					gd.RhoTheta_flag=0;
					free_memory();
					free(bodytab);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		gd.RhoTheta_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					EvalRhoTheta(fname, fnametmp);
					if (gd.RhoTheta_flag) {
						animation_rhotheta(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
					}
					gd.RhoTheta_flag=0;
					free_memory();
					free(bodytab);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}

void general_fx_anim_plplot(string fname, string snapnametmp)
{
    struct stat buf;
    stream instr;

	int step=0, moreSteps=1, filecount=0;
	bool exist_file;
    char namebuf[200], namebuf1[200];

	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (scanopt(cmd.options, "save")) {
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				sprintf(namebuf, cmd.in, step);
				if (stat(namebuf, &buf) != 0)
					exist_file = FALSE;
				else
					exist_file = TRUE;
				if (exist_file) {
					++filecount;
						sprintf(namebuf1, cmd.out, filecount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
						if (!strnull(cmd.pl_a))
							plsetopt("-a",cmd.pl_a);
						plinit();
						animation_general_fx(cmd.in, step);
						free(npltd.xval); free(npltd.yval);
						plend();
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				sprintf(namebuf, cmd.in, step);
				if (stat(namebuf, &buf) != 0)
					exist_file = FALSE;
				else
					exist_file = TRUE;
				if (exist_file) {
					++filecount;
					animation_general_fx(cmd.in, step);
						free(npltd.xval); free(npltd.yval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}

// ------------------ RHO ANIM PLPLOT -----------------------------------------
void rho_anim_plplot(string fname, string fnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	double cpustart;

	cpustart = cputime();

	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);
	if (!strnull(cmd.pl_geo))
		plsetopt("-geo",cmd.pl_geo);
	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);
	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);
	if (!strnull(cmd.pl_bg))
		plsetopt("-bg",cmd.pl_bg);
	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);
	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (scanopt(cmd.options, "save")) {
		gd.Rho_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					EvalRhoDist(fname, fnametmp);
					if (gd.Rho_flag) {
						sprintf(namebuf1, cmd.out, snapcount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
//						plssub(1, 4);
						plinit();
						animation_rho(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
						plend();
					}
					gd.Rho_flag=0;
					free_memory();
					free(bodytab);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		gd.Rho_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					EvalRhoDist(fname, fnametmp);
					if (gd.Rho_flag) {
						animation_rho(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
					}
					gd.Rho_flag=0;
					free_memory();
					free(bodytab);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}
// -----------------------------------------------------------------------------


// ------------------ VC ANIM PLPLOT -----------------------------------------
// -----------------------------------------------------------------------------

local void animation_xy(string snapname, int snapcount)
{
  int i, col1, col2, npoint;
  real *xnc[1], *ync[1];
  int axestype;
//  int axeswidth;
  int axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
//  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor;
//  int labelfontweight;
//  int symbolcolor, symbolweight; //, symboltype;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  real viewportchr;

	printf("\n\nX-Y Animation ...\n\n");

	col1=gd.column1;
	col2=gd.column2;

	readin_pl_snap(snapname, col1, col2, &npoint);

	printf("\nnpoint = %d\n",npoint);

	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];

	if (gd.x_autoscale) {
			gd.xmin = xnc[0][0];
			gd.xmax = xnc[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.xmin>xnc[0][i]) gd.xmin=xnc[0][i];
				if (gd.xmax<xnc[0][i]) gd.xmax=xnc[0][i];
			}
	}

	if (gd.y_autoscale) {
			gd.ymin = ync[0][0];
			gd.ymax = ync[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.ymin>ync[0][i]) gd.ymin=ync[0][i];
				if (gd.ymax<ync[0][i]) gd.ymax=ync[0][i];
			}
	}

	printf("\nscaling done ...\n");

	axestype = 1;
//	if (scanopt(cmd.options, "save"))
//		axeswidth = 2;
//	else
//		axeswidth = 1;
	axescolor = 1;

	pllsty(axestype);
	plwid(cmd.axeswidth);
    plcol0(axescolor);

	pladv(1);
	plclear();
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, viewportchr);		// Set the scaled size of the labels
	plvsta();
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
//	nlsize=1.0;

    plsxax(nlxdigmax, 0);
    plsyax(nlydigmax, 0);
	plcol0(nlcolor);
	plschr(0.0,cmd.nlsize);

	frame=TRUE;
	xaxis=FALSE;
	yaxis=FALSE;
	frame_xopt="bcnst";
	frame_xtick=0.0;
	frame_nxsub=0;
	frame_yopt="bcnstv";
	frame_ytick=0.0;
	frame_nysub=0;
	xorigin=0.;
	yorigin=0.;

	if (frame)
		plbox(frame_xopt,frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);
	if (xaxis)
		plaxes(xorigin,yorigin,frame_xopt,frame_xtick,frame_nxsub,"",frame_ytick,frame_nysub);
	if (yaxis)
		plaxes(xorigin,yorigin,"",frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);

	plwid(defpenwidth); 

	labelcolor=2;
//	if (scanopt(cmd.options, "save"))
//		labelfontweight=2;
//	else
//		labelfontweight=1;

    plcol0(labelcolor);
	plschr(0.0, cmd.labelfontsize);
	plwid(cmd.labelfontweight);
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

//	symbolcolor=3;
//	if (scanopt(cmd.options, "save"))
//		symbolweight=2;
//	else
//		symbolweight=1;
//	symboltype=4;

	if (cmd.withsymbols) {
		plcol0(cmd.symbolcolor);
		plssym(0.0, cmd.symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint, xnc[0], ync[0], cmd.symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(cmd.symbolcolor);
		plssym(0.0, 0.0);
		plwid(cmd.symbolweight);
		plpoin(npoint, xnc[0], ync[0], 1);
		plwid(defpenwidth); 
	}

    sprintf(timebuf, "%g%", gd.tnow);
	timecolor=3;
	timesize=1.0;
	if (scanopt(cmd.options, "save"))
		timeweight=2;
	else
		timeweight=1;
	timeside="t";
	timedisp=1.0;
	timepos=0.75;
	timejust=0;
	plcol0(timecolor);
	plschr(0.0, timesize);
	plwid(timeweight);
	plmtex(timeside,timedisp,timepos,timejust,timebuf);
	plwid(defpenwidth);

	copyrighttext(TRUE);

	printf("\n\n... done snap animation plotting! ...\n\n");
}

local void animation_xy_type(string snapname, int snapcount)
{
#if !defined(BODY1)
#define BODY1	1
#endif
#if !defined(BODY2)
#define BODY2	2
#endif
  int i, col1, col2, col3, npoint;
  real *xnc[1], *ync[1];
  int *ival[1];
  real *x1, *y1, *x2, *y2;
  int npoint1, npoint2, i1, i2;
  int axestype, axeswidth, axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor, labelfontweight;
  int symbolcolor1, symbolcolor2;
  int symbolweight; //, symboltype;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;

	printf("\n\nX-Y Animation type ...\n\n");

	col1=gd.column1;
	col2=gd.column2;
	col3=7;

	readin_pl_snap_type(snapname, col1, col2, col3, &npoint);

	printf("\nnpoint = %d\n",npoint);

	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0]; ival[0]=&npltd.Typeval[0];

	if (gd.x_autoscale) {
			gd.xmin = xnc[0][0];
			gd.xmax = xnc[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.xmin>xnc[0][i]) gd.xmin=xnc[0][i];
				if (gd.xmax<xnc[0][i]) gd.xmax=xnc[0][i];
			}
	}

	if (gd.y_autoscale) {
			gd.ymin = ync[0][0];
			gd.ymax = ync[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.ymin>ync[0][i]) gd.ymin=ync[0][i];
				if (gd.ymax<ync[0][i]) gd.ymax=ync[0][i];
			}
	}

	npoint1=npoint2=0;
	for (i = 0; i < npoint; i++) {
		if (ival[0][i]==BODY1) ++npoint1;
		if (ival[0][i]==BODY2) ++npoint2;
	}
	fprintf(stdout,"\nReserving memory for %d %d points both types\n",
		npoint1, npoint2);
	x1 = (real *) allocate(npoint1 * sizeof(real));
	y1 = (real *) allocate(npoint1 * sizeof(real));
	x2 = (real *) allocate(npoint2 * sizeof(real));
	y2 = (real *) allocate(npoint2 * sizeof(real));
	i1=i2=0;
	for (i = 0; i < npoint; i++) {
		if (ival[0][i]==BODY1) {
			x1[i1]=xnc[0][i];
			y1[i1]=ync[0][i];
			++i1;
		}
		if (ival[0][i]==BODY2) {
			x2[i2]=xnc[0][i];
			y2[i2]=ync[0][i];
			++i2;
		}
	}
	fprintf(stdout,"Total points assigned %d %d both types\n",
		i1, i2);
//#undef BODY1
//#undef BODY2

	printf("\nscaling done ...\n");

	axestype = 1;
	if (scanopt(cmd.options, "save"))
		axeswidth = 2;
	else
		axeswidth = 1;
	axescolor = 1;

	pllsty(axestype);
	plwid(axeswidth);
    plcol0(axescolor);

	pladv(1);
	plclear();
	plvsta();
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
	nlsize=1.0;

    plsxax(nlxdigmax, 0);
    plsyax(nlydigmax, 0);
	plcol0(nlcolor);
	plschr(0.0,nlsize);

	frame=TRUE;
	xaxis=FALSE;
	yaxis=FALSE;
	frame_xopt="bcnst";
	frame_xtick=0.0;
	frame_nxsub=0;
	frame_yopt="bcnstv";
	frame_ytick=0.0;
	frame_nysub=0;
	xorigin=0.;
	yorigin=0.;

	if (frame)
		plbox(frame_xopt,frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);
	if (xaxis)
		plaxes(xorigin,yorigin,frame_xopt,frame_xtick,frame_nxsub,"",frame_ytick,frame_nysub);
	if (yaxis)
		plaxes(xorigin,yorigin,"",frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);

	plwid(defpenwidth); 

	labelcolor=2;
	if (scanopt(cmd.options, "save"))
		labelfontweight=2;
	else
		labelfontweight=1;

    plcol0(labelcolor);
	plschr(0.0, cmd.labelfontsize);
	plwid(labelfontweight);
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	symbolcolor1=2;
	symbolcolor2=3;
	if (scanopt(cmd.options, "save"))
		symbolweight=2;
	else
		symbolweight=1;
//	symboltype=4;

/*
	if (cmd.withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, cmd.symbolsize);
		plwid(symbolweight);
		plpoin(npoint, xnc[0], ync[0], symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
	}
*/
	if (cmd.withsymbols) {
		plcol0(symbolcolor1);
		plssym(0.0, cmd.symbolsize);
		plwid(symbolweight);
		plpoin(npoint1, x1, y1, cmd.symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor1);
		plssym(0.0, 0.0);
		plpoin(npoint1, x1, y1, 1);
	}

	if (cmd.withsymbols) {
		plcol0(symbolcolor2);
		plssym(0.0, cmd.symbolsize);
		plwid(symbolweight);
		plpoin(npoint2, x2, y2, cmd.symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor2);
		plssym(0.0, 0.0);
		plpoin(npoint2, x2, y2, 1);
	}

    sprintf(timebuf, "%g%", gd.tnow);
	timecolor=3;
	timesize=1.0;
	if (scanopt(cmd.options, "save"))
		timeweight=2;
	else
		timeweight=1;
	timeside="t";
	timedisp=1.0;
	timepos=0.75;
	timejust=0;
	plcol0(timecolor);
	plschr(0.0, timesize);
	plwid(timeweight);
	plmtex(timeside,timedisp,timepos,timejust,timebuf);
	plwid(defpenwidth);

	copyrighttext(TRUE);

	printf("\n\n... done snap animation type plotting! ...\n\n");
}

local void animation_xy_trajectory(string snapname, int snapcount)
{
  int i, col1, col2, npoint;
  real *xnc[1], *ync[1];
  int axestype, axeswidth, axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor, labelfontweight;
  real labelfontsize;
  string xlabel, ylabel, plotlabel;
  bool withsymbols, withdots;
  int symbolcolor, symbolweight; //, symboltype;
  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool f_name;
  PLINT cur_strm, new_strm;
  char   namebuf[200], namebuf1[200];

	printf("\n\nX-Y Animation (Trajectories) ...\n\n");

//	col1=1;
//	col2=2;
	col1=gd.column1;
	col2=gd.column2;

	readin_pl_snap(snapname, col1, col2, &npoint);

	printf("\nnpoint = %d\n",npoint);

	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];

	if (gd.x_autoscale) {
			gd.xmin = xnc[0][0];
			gd.xmax = xnc[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.xmin>xnc[0][i]) gd.xmin=xnc[0][i];
				if (gd.xmax<xnc[0][i]) gd.xmax=xnc[0][i];
			}
	}

	if (gd.y_autoscale) {
			gd.ymin = ync[0][0];
			gd.ymax = ync[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.ymin>ync[0][i]) gd.ymin=ync[0][i];
				if (gd.ymax<ync[0][i]) gd.ymax=ync[0][i];
			}
	}

/*	
	if (snapcount==1) {
		gd.xmin=-0.5*gd.Box[0];
		gd.xmax=0.5*gd.Box[0];
		gd.ymin=-0.5*gd.Box[1];
		gd.ymax=0.5*gd.Box[1];
	}
*/
	printf("\nscaling done ...\n");

	if (snapcount==1) {

	axestype = 1;
	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);							// Line type (8 styles). See 'plstyl'
	plwid(axeswidth);							// Line width
    plcol0(axescolor);							// Set color of frame - viewport

	pladv(1);									// Set to zero to advance to the next subpage
	plvsta();									// Sets up a standard viewport, leaving apropriate margins
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);				// Sets up the world-coordinates windows

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
	nlsize=1.0;

    plsxax(nlxdigmax, 0);						// To set digmax on x-axis
    plsyax(nlydigmax, 0);						// To set digmax on y-axis
	plcol0(nlcolor);							// Set color for numeric labels. Also set color for axes...
	plschr(0.0,nlsize);

	frame=TRUE;
	xaxis=FALSE;
	yaxis=FALSE;
	frame_xopt="bcnst";
	frame_xtick=0.0;
	frame_nxsub=0;
	frame_yopt="bcnstv";
	frame_ytick=0.0;
	frame_nysub=0;
	xorigin=0.;
	yorigin=0.;

	if (frame)
		plbox(frame_xopt,frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);
	if (xaxis)
		plaxes(xorigin,yorigin,frame_xopt,frame_xtick,frame_nxsub,"",frame_ytick,frame_nysub);
	if (yaxis)
		plaxes(xorigin,yorigin,"",frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);

	plwid(defpenwidth); 

	labelcolor=2;
	labelfontsize=1.;
	labelfontweight=2;
	xlabel="x";
	ylabel="y";
	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, labelfontsize);
	plwid(labelfontweight);
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	withsymbols=FALSE;
	withdots=TRUE;
	symbolcolor=3;
	symbolsize=0.5;
	symbolweight=1;
//	symboltype=4;

	}

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, symbolsize);
		plwid(symbolweight);
		plpoin(npoint, xnc[0], ync[0], cmd.symboltype);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
	}

/*
    sprintf(timebuf, "%g%", tnow);
	timecolor=3;
	timesize=1.0;
	timeweight=1;
	timeside="t";
	timedisp=1.0;
	timepos=0.75;
	timejust=0;
	plcol0(timecolor);
	plschr(0.0, timesize);
	plwid(timeweight);
	plmtex(timeside,timedisp,timepos,timejust,timebuf);
	plwid(defpenwidth);
*/

// No se despliega la informacion del tiempo porque se encima el texto...

	copyrighttext(TRUE);

	if (scanopt(cmd.options, "save")) {
		sprintf(namebuf1, cmd.out, snapcount);
		sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
		printf("The current plot was saved in color Postscript under the name `%s'.\n", namebuf);
		plgstrm(&cur_strm);
		plmkstrm(&new_strm);
		plsfnam(namebuf);
		plsdev("psc");
		plcpstrm(cur_strm, 0);
		plreplot();
		plend1();
		plsstrm(cur_strm);
    }

	printf("\n\n... done snap trajectory plotting! ...\n\n");
}

local void animation_vcr(string filename, int snapcount)
{
  int i, col1, col2, npoint;
  real *xnc[2], *ync[2];
  int axestype, axeswidth, axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor, labelfontweight;
  real labelfontsize;
  string xlabel, ylabel, plotlabel;
  bool plotjoined, withsymbols, withdots; 
  int linecolor, linewidth, linetype;
  int symbolcolor, symbolweight, symboltype;
  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool printcopyright;

	printf("\n\nVcR Animation ...\n\n");

	col1=1;
	col2=2;

	readin_pl_snap(filename, col1, col2, &npoint);
	printf("\nnpoint = %d\n",npoint);
	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];	
	if (gd.x_autoscale) {
		gd.xmin = xnc[0][0];
		gd.xmax = xnc[0][0];
		for (i = 1; i < npoint; i++) {
			if (gd.xmin>xnc[0][i]) gd.xmin=xnc[0][i];
			if (gd.xmax<xnc[0][i]) gd.xmax=xnc[0][i];
		}
	}
	if (gd.y_autoscale) {
		gd.ymin = ync[0][0];
		gd.ymax = ync[0][0]; 
		for (i = 1; i < npoint; i++) {
			if (gd.ymin>ync[0][i]) gd.ymin=ync[0][i];
			if (gd.ymax<ync[0][i]) gd.ymax=ync[0][i];
		}
	}

	printf("\nscaling done ...\n");

	axestype = 1;
	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(axeswidth);
    plcol0(axescolor);

	pladv(1);
	plclear();
	plvsta();
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
	nlsize=1.0;

    plsxax(nlxdigmax, 0);
    plsyax(nlydigmax, 0);
	plcol0(nlcolor);
	plschr(0.0,nlsize);

	frame=TRUE;
	xaxis=FALSE;
	yaxis=FALSE;
	frame_xopt="bcnst";
	frame_xtick=0.0;
	frame_nxsub=0;
	frame_yopt="bcnstv";
	frame_ytick=0.0;
	frame_nysub=0;
	xorigin=0.;
	yorigin=0.;

	if (frame)
		plbox(frame_xopt,frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);
	if (xaxis)
		plaxes(xorigin,yorigin,frame_xopt,frame_xtick,frame_nxsub,"",frame_ytick,frame_nysub);
	if (yaxis)
		plaxes(xorigin,yorigin,"",frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);

	plwid(defpenwidth);

	labelcolor=2;
	labelfontsize=1.2;
	labelfontweight=2;
	xlabel="r";
	ylabel="Circular Velocity";
	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, labelfontsize);
	plwid(labelfontweight);
    pllab(xlabel, ylabel, plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	plotjoined=TRUE;
	withsymbols=FALSE;
	withdots=FALSE;

	linecolor=3;
	linewidth=3;
	linetype=1;

	symbolcolor=3;
	symbolsize=0.5;
	symbolweight=2;
	symboltype=4;

	if (plotjoined) {
		plcol0(linecolor);
		pllsty(linetype);
		plwid(linewidth);
		plline(npoint, xnc[0], ync[0]);
		plwid(defpenwidth); 
	}

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, symbolsize);
		plwid(symbolweight);
		plpoin(npoint, xnc[0], ync[0], symboltype);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
	}

    sprintf(timebuf, "%g%", gd.tnow);
	timecolor=3;
	timesize=1.0;
	timeweight=2;
	timeside="t";
	timedisp=1.0;
	timepos=0.75;
	timejust=0;
	plcol0(timecolor);
	plschr(0.0, timesize);
	plwid(timeweight);
	plmtex(timeside,timedisp,timepos,timejust,timebuf);
	plwid(defpenwidth);

	printcopyright=TRUE;

	copyrighttext(printcopyright);

	printf("\n\n... done vcr plotting! ...\n\n");
}

local void animation_rhotheta(string rhofilename, int snapcount)
{
  int i, col1, col2, npoint;
  real *xnc[2], *ync[2];
  int axestype, axeswidth, axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor, labelfontweight;
  real labelfontsize;
  string xlabel, ylabel, plotlabel;
  bool plotjoined, withsymbols, withdots; 
  int linecolor, linewidth, linetype;
  int symbolcolor, symbolweight; //, symboltype;
  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool printcopyright;
  int ifile;
  int nmaxfiles=4;		// Old value 3...

	printf("\n\nRHOTHETA Animation ...\n\n");

for (ifile=0; ifile<nmaxfiles; ifile++) {

	col1=(ifile+1)*2-1;
	col2=(ifile+1)*2;

	readin_pl_snap(rhofilename, col1, col2, &npoint);
	printf("\nnpoint = %d\n",npoint);
	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];	
	if (gd.x_autoscale) {
		gd.xmin = xnc[0][0];
		gd.xmax = xnc[0][0];
		for (i = 1; i < npoint; i++) {
			if (gd.xmin>xnc[0][i]) gd.xmin=xnc[0][i];
			if (gd.xmax<xnc[0][i]) gd.xmax=xnc[0][i];
		}
	}

	if (gd.y_autoscale) {
		gd.ymin = ync[0][0];
		gd.ymax = ync[0][0]; 
		for (i = 1; i < npoint; i++) {
			if (gd.ymin>ync[0][i]) gd.ymin=ync[0][i];
			if (gd.ymax<ync[0][i]) gd.ymax=ync[0][i];
		}
	} else {
		gd.ymin = gd.yrangeMin[ifile];
		gd.ymax = gd.yrangeMax[ifile];
	}

	printf("\nscaling done ...\n");

	axestype = 1;
	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(axeswidth);
    plcol0(axescolor);

	if (ifile==0)
		pladv(1);
	else
		pladv(0);
	plclear();
	plvsta();
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
	nlsize=1.0;

    plsxax(nlxdigmax, 0);
    plsyax(nlydigmax, 0);
	plcol0(nlcolor);
	plschr(0.0,nlsize);

	frame=TRUE;
	xaxis=FALSE;
	yaxis=FALSE;
	frame_xopt="bcnst";
	frame_xtick=0.0;
	frame_nxsub=0;
	frame_yopt="bcnstv";
	frame_ytick=0.0;
	frame_nysub=0;
	xorigin=0.;
	yorigin=0.;

	if (frame)
		plbox(frame_xopt,frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);
	if (xaxis)
		plaxes(xorigin,yorigin,frame_xopt,frame_xtick,frame_nxsub,"",frame_ytick,frame_nysub);
	if (yaxis)
		plaxes(xorigin,yorigin,"",frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);

	plwid(defpenwidth);

	labelcolor=2;
	labelfontsize=1.2;
	labelfontweight=2;
	if (ifile==0) xlabel="#gf";
	if (ifile==1) xlabel="#gf";
	if (ifile==2) xlabel="#gf";
	if (ifile==3) xlabel="#gf";
	if (ifile==4) xlabel="#gf";
	if (ifile==5) xlabel="#gf";
//	if (ifile==0) ylabel="Density Bulge";
	if (ifile==0) ylabel="Densidad Disco";
	if (ifile==1) ylabel="Densidad Halo";
//	if (ifile==3) ylabel="Acc Rad";
	if (ifile==2) ylabel="Acc Tan";
	if (ifile==3) ylabel="Vel Tan";
//	ylabel="Density profile";
	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, labelfontsize);
	plwid(labelfontweight);
    pllab(xlabel, ylabel, plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	plotjoined=TRUE;
	withsymbols=FALSE;
	withdots=FALSE;

	linecolor=3;
	linewidth=3;
	linetype=1;

	symbolcolor=3;
	symbolsize=0.5;
	symbolweight=2;
//	symboltype=4;

	if (plotjoined) {
		plcol0(linecolor);
		pllsty(linetype);
		plwid(linewidth);
		plline(npoint, xnc[0], ync[0]);
		plwid(defpenwidth); 
	}

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, symbolsize);
		plwid(symbolweight);
		plpoin(npoint, xnc[0], ync[0], cmd.symboltype);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
	}

    sprintf(timebuf, "%g%", gd.tnow);
	timecolor=3;
	timesize=1.0;
	timeweight=2;
	timeside="t";
	timedisp=1.0;
	timepos=0.75;
	timejust=0;
	plcol0(timecolor);
	plschr(0.0, timesize);
	plwid(timeweight);
	plmtex(timeside,timedisp,timepos,timejust,timebuf);
	plwid(defpenwidth);

	printcopyright=TRUE;

	copyrighttext(printcopyright);
}
	printf("\n\n... done rhotheta plotting! ...\n\n");
}


local void animation_rho(string rhofilename, int snapcount)
{
  int i, col1, col2, npoint;
  real *xnc[2], *ync[2];
  int axestype, axeswidth, axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor, labelfontweight;
  real labelfontsize;
  string xlabel, ylabel, plotlabel;
  bool plotjoined, withsymbols, withdots; 
  int linecolor, linewidth, linetype;
  int symbolcolor, symbolweight; //, symboltype;
  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool printcopyright;
  int ifile;
  int nmaxfiles=1;		// Old value 3...
  real tmp;

	printf("\n\nRHO Animation ...\n\n");

for (ifile=0; ifile<nmaxfiles; ifile++) {

	col1=(ifile+1)*2-1;
	col2=(ifile+1)*2;

	readin_pl_snap(rhofilename, col1, col2, &npoint);
	printf("\nnpoint = %d\n",npoint);
	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];	
	if (gd.x_autoscale) {
		gd.xmin = xnc[0][0];
		gd.xmax = xnc[0][0];
		for (i = 1; i < npoint; i++) {
			if (gd.xmin>xnc[0][i]) gd.xmin=xnc[0][i];
			if (gd.xmax<xnc[0][i]) gd.xmax=xnc[0][i];
		}
	}

	if (gd.y_autoscale) {
		gd.ymin = ync[0][0];
		gd.ymax = ync[0][0]; 
		for (i = 1; i < npoint; i++) {
			if (gd.ymin>ync[0][i]) gd.ymin=ync[0][i];
			if (gd.ymax<ync[0][i]) gd.ymax=ync[0][i];
		}
	} else {
		gd.ymin = gd.yrangeMin[ifile];
		gd.ymax = gd.yrangeMax[ifile];
	}

	for (i=0; i<npoint; i++) {								// Log plot ...
		tmp = rlog10(ync[0][i]);
		ync[0][i] = tmp;
		tmp = rlog10(xnc[0][i]);
		xnc[0][i] = tmp;
	}

	gd.xmin = rlog10(gd.xmin);
	gd.xmax = rlog10(gd.xmax);
	gd.ymin = rlog10(gd.ymin);
	gd.ymax = rlog10(gd.ymax);

	printf("\nscaling done ...\n");

	axestype = 1;
	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(axeswidth);
    plcol0(axescolor);

	if (ifile==0)
		pladv(1);
	else
		pladv(0);
	plclear();
	plvsta();
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
	nlsize=1.0;

    plsxax(nlxdigmax, 0);
    plsyax(nlydigmax, 0);
	plcol0(nlcolor);
	plschr(0.0,nlsize);

	frame=TRUE;
	xaxis=FALSE;
	yaxis=FALSE;
//	frame_xopt="bcnst";
	frame_xopt="bcnstl";								// Log plot ...
	frame_xtick=0.0;
	frame_nxsub=0;
//	frame_yopt="bcnstv";
	frame_yopt="bcnstvl";								// Log plot ...
	frame_ytick=0.0;
	frame_nysub=0;
	xorigin=0.;
	yorigin=0.;

	if (frame)
		plbox(frame_xopt,frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);
	if (xaxis)
		plaxes(xorigin,yorigin,frame_xopt,frame_xtick,frame_nxsub,"",frame_ytick,frame_nysub);
	if (yaxis)
		plaxes(xorigin,yorigin,"",frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);

	plwid(defpenwidth);

	labelcolor=2;
	labelfontsize=1.2;
	labelfontweight=2;
	if (ifile==0) xlabel="r";
	if (ifile==1) xlabel="#gf";
	if (ifile==2) xlabel="#gf";
	if (ifile==3) xlabel="#gf";
	if (ifile==4) xlabel="#gf";
	if (ifile==5) xlabel="#gf";
	if (ifile==0) ylabel="Density profile";
//	if (ifile==0) ylabel="Densidad Disco";
	if (ifile==1) ylabel="Densidad Halo";
//	if (ifile==3) ylabel="Acc Rad";
	if (ifile==2) ylabel="Acc Tan";
	if (ifile==3) ylabel="Vel Tan";
//	ylabel="Density profile";
	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, labelfontsize);
	plwid(labelfontweight);
    pllab(xlabel, ylabel, plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	plotjoined=TRUE;
	withsymbols=FALSE;
	withdots=FALSE;

	linecolor=3;
	linewidth=3;
	linetype=1;

	symbolcolor=3;
	symbolsize=0.5;
	symbolweight=2;
//	symboltype=4;

	if (plotjoined) {
		plcol0(linecolor);
		pllsty(linetype);
		plwid(linewidth);
		plline(npoint, xnc[0], ync[0]);
		plwid(defpenwidth); 
	}

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, symbolsize);
		plwid(symbolweight);
		plpoin(npoint, xnc[0], ync[0], cmd.symboltype);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
	}

    sprintf(timebuf, "%g%", gd.tnow);
	timecolor=3;
	timesize=1.0;
	timeweight=2;
	timeside="t";
	timedisp=1.0;
	timepos=0.75;
	timejust=0;
	plcol0(timecolor);
	plschr(0.0, timesize);
	plwid(timeweight);
	plmtex(timeside,timedisp,timepos,timejust,timebuf);
	plwid(defpenwidth);

	printcopyright=TRUE;

	copyrighttext(printcopyright);
}
	printf("\n\n... done rho plotting! ...\n\n");
}

local void animation_general_fx(string rdffilename, int snapcount)
{
    char namebuf[256];

  int i, col1, col2, npoint;
  real *xnc[2], *ync[2];
  int axestype, axeswidth, axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor, labelfontweight;
//  real labelfontsize;
//  string xlabel, ylabel, plotlabel;
  bool plotjoined;
//  bool withsymbols, withdots; 
  int linecolor, linewidth, linetype;
  int symbolcolor, symbolweight; //, symboltype;
//  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool printcopyright;

	printf("\n\nGeneral fx Animation ...\n\n");

//	col1=1;
//	col2=2;
	col1=gd.column1;
	col2=gd.column2;

    sprintf(namebuf, cmd.in, snapcount);

	readin_pl_snap(namebuf, col1, col2, &npoint);
//	readin_pl_snap(rdffilename, col1, col2, &npoint);
	printf("\nnpoint = %d\n",npoint);
	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];	
	if (gd.x_autoscale) {
		gd.xmin = xnc[0][0];
		gd.xmax = xnc[0][0];
		for (i = 1; i < npoint; i++) {
			if (gd.xmin>xnc[0][i]) gd.xmin=xnc[0][i];
			if (gd.xmax<xnc[0][i]) gd.xmax=xnc[0][i];
		}
	}
	if (gd.y_autoscale) {
		gd.ymin = ync[0][0];
		gd.ymax = ync[0][0]; 
		for (i = 1; i < npoint; i++) {
			if (gd.ymin>ync[0][i]) gd.ymin=ync[0][i];
			if (gd.ymax<ync[0][i]) gd.ymax=ync[0][i];
		}
	}

	printf("\nscaling done ...\n");

	axestype = 1;
	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(axeswidth);
    plcol0(axescolor);

	pladv(1);
	plclear();
	plvsta();
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
	nlsize=1.0;

    plsxax(nlxdigmax, 0);
    plsyax(nlydigmax, 0);
	plcol0(nlcolor);
	plschr(0.0,nlsize);

	frame=TRUE;
	xaxis=FALSE;
	yaxis=FALSE;
	frame_xopt="bcnst";
	frame_xtick=0.0;
	frame_nxsub=0;
	frame_yopt="bcnstv";
	frame_ytick=0.0;
	frame_nysub=0;
	xorigin=0.;
	yorigin=0.;

	if (frame)
		plbox(frame_xopt,frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);
	if (xaxis)
		plaxes(xorigin,yorigin,frame_xopt,frame_xtick,frame_nxsub,"",frame_ytick,frame_nysub);
	if (yaxis)
		plaxes(xorigin,yorigin,"",frame_xtick,frame_nxsub,frame_yopt,frame_ytick,frame_nysub);

	plwid(defpenwidth);

	labelcolor=2;
//	labelfontsize=1.2;
	labelfontweight=2;
//	xlabel="spatial coordinate bin";
//	ylabel="Profile";
//	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, cmd.labelfontsize);
	plwid(labelfontweight);
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	plotjoined=TRUE;
//	withsymbols=FALSE;
//	withdots=FALSE;

	linecolor=3;
	linewidth=3;
	linetype=1;

	symbolcolor=3;
//	symbolsize=0.5;
	symbolweight=2;
//	symboltype=4;

	if (plotjoined) {
		plcol0(linecolor);
		pllsty(linetype);
		plwid(linewidth);
		plline(npoint, xnc[0], ync[0]);
		plwid(defpenwidth); 
	}

	if (cmd.withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, cmd.symbolsize);
		plwid(symbolweight);
		plpoin(npoint, xnc[0], ync[0], cmd.symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
	}

//    sprintf(timebuf, "%g%", tnow);
    sprintf(timebuf, "%d%", snapcount);
	timecolor=3;
	timesize=1.0;
	timeweight=2;
	timeside="t";
	timedisp=1.0;
	timepos=0.75;
	timejust=0;
	plcol0(timecolor);
	plschr(0.0, timesize);
	plwid(timeweight);
	plmtex(timeside,timedisp,timepos,timejust,timebuf);
	plwid(defpenwidth);

	printcopyright=TRUE;

	copyrighttext(printcopyright);

	printf("\n\n... done general fx plotting! ...\n\n");
}

// DIBUJO DE LEGENDAS PARA LOS SIMBOLOS

#define TOPLEFT		0
#define TOPRIGHT	1
#define BOTTONLEFT	2
#define BOTTONRIGHT	3

local void setlegends(int nlegends, char **legendnames, int legendspos,
	bool *lplotjoined, bool *lwithsymbols, bool *lwithdots, int *lsymbolcolor, 
	real *lsymbolsize, int *lsymbolweight, int *lsymboltype, int *llinecolor,
	int *llinetype, int *llinewidth) {

  int i, j;
  real xu, yu, yfs, xstart, ystart, lsize, dyi, fl, fs, lwx, lwy;
  real xsym[1], ysym[1];

// Setting reference point
	xu = (gd.xmax-gd.xmin)/10.;
	yu = (gd.ymax-gd.ymin)/10.;
	yfs = 2.0*(gd.ymax-gd.ymin)/32.4;
	xstart = gd.xmin+1.0*xu;
	ystart = gd.ymax-0.6*yu;
	lsize = xu/2.0;
	fl = 0.65;
	fs =1.0;
	lwx = 5.0*xu;
	lwy = 5.0*yu;

    for (i = 0; i < nlegends; i++) {

	switch(legendspos) {

		case TOPLEFT:
			dyi = (real)(i)*yfs;
			if (lplotjoined[i]) {
				plcol0(llinecolor[i]);
				pllsty(llinetype[i]);
				plwid(llinewidth[i]);
				pljoin(xstart, ystart-dyi, xstart+lsize, ystart-dyi);
				plwid(defpenwidth); 
			}
			if (lwithsymbols[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, lsymbolsize[i]);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-dyi;
				plwid(lsymbolweight[i]);
				plpoin(1, xsym, ysym, lsymboltype[i]);
				plwid(defpenwidth); 
			}
			if (lwithdots[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, 0.0);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-dyi;
				plpoin(1, xsym, ysym, 1);
			}
			plptex(xstart+fs*xu, ystart-dyi, 0.1, 0.0, 0., legendnames[i]);
			break;

		case TOPRIGHT:
			dyi = (real)(i)*yfs;
			if (lplotjoined[i]) {
				plcol0(llinecolor[i]);
				pllsty(llinetype[i]);
				plwid(llinewidth[i]);
				pljoin(xstart+lwx,		 ystart-dyi, xstart+lwx+lsize, ystart-dyi);
				plwid(defpenwidth); 
			}
			if (lwithsymbols[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, lsymbolsize[i]);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-dyi;
				plwid(lsymbolweight[i]);
				plpoin(1, xsym, ysym, lsymboltype[i]);
				plwid(defpenwidth);
			}
			if (lwithdots[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, 0.0);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-dyi;
				plpoin(1, xsym, ysym, 1);
			}
			plptex(xstart+lwx+fs*xu, ystart-dyi, 0.1, 0.0, 0., legendnames[i]);
			break;


		case BOTTONLEFT:
			dyi = (real)(i)*yfs;
			if (lplotjoined[i]) {
				plcol0(llinecolor[i]);
				pllsty(llinetype[i]);
				plwid(llinewidth[i]);
				pljoin(xstart,		 ystart-lwy-dyi, xstart+lsize, ystart-lwy-dyi);
				plwid(defpenwidth);
			}
			if (lwithsymbols[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, lsymbolsize[i]);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plwid(lsymbolweight[i]);
				plpoin(1, xsym, ysym, lsymboltype[i]);
				plwid(defpenwidth);
			}
			if (lwithdots[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, 0.0);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plpoin(1, xsym, ysym, 1);
			}
			plptex(xstart+fs*xu, ystart-lwy-dyi, 0.1, 0.0, 0., legendnames[i]);
			break;

		case BOTTONRIGHT:
			dyi = (real)(i)*yfs;
			if (lplotjoined[i]) {
				plcol0(llinecolor[i]);
				pllsty(llinetype[i]);
				plwid(llinewidth[i]);
				pljoin(xstart+lwx,		 ystart-lwy-dyi, xstart+lwx+lsize, ystart-lwy-dyi);
				plwid(defpenwidth);
			}
			if (lwithsymbols[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, lsymbolsize[i]);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plwid(lsymbolweight[i]);
				plpoin(1, xsym, ysym, lsymboltype[i]);
				plwid(defpenwidth);
			}
			if (lwithdots[i]) {
				plcol0(lsymbolcolor[i]);
				plssym(0.0, 0.0);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plpoin(1, xsym, ysym, 1);
			}
			plptex(xstart+lwx+fs*xu, ystart-lwy-dyi, 0.1, 0.0, 0., legendnames[i]);
			break;
		}
    }
}

#undef TOPLEFT
#undef TOPRIGHT
#undef BOTTONLEFT
#undef BOTTONRIGHT

local void copyrighttext(bool cpr)
{
  string cprtext;
  int cprcolor, cprweight;
  real cprsize, cprdisp, cprpos, cprjust;
  string cprside;

	if (cpr) {
		plfontld(1);
		cprtext="#(274)2007 Mar";
		cprcolor=8;
		cprsize=0.5;
		cprweight=1;
		cprside="b";
		cprdisp=8;
		cprpos=0.9;
		cprjust=0;
		plcol0(cprcolor);
		plschr(0.0, cprsize);
		plwid(cprweight);
		plmtex(cprside,cprdisp,cprpos,cprjust,cprtext);
		plwid(defpenwidth);
		plfontld(0);
		plschr(0.0, 1.0);
	}
}

