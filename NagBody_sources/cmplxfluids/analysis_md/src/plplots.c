/*==============================================================================
	MODULE: plplots.c				[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date:	May 2006
	Purpose:
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes:
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: November 2008;
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"


local void animation_xy(string, int);
local void animation_xy_type(string, int);
local void animation_xy_type_4b(string, int);
#ifdef THREEDIM
local void animation_xyz_3d(string, int);
local void animation_xy_trajectory_3d(string, int);
#endif
local void animation_xy_trajectory(string, int);
local void animation_xy_trajectory_type(string, int);
local void animation_rdf(string, int);
local void animation_vel(string, int);
local void animation_rhoaxes(string, int);
local void animation_rhoaxes_type(string, int);
#ifdef THREEDIM
local void animation_nfrecaxes(string, string, string, int);
#else
local void animation_nfrecaxes(string, string, int);
#endif
local void animation_general_fx(string, int);

local void setlegends(int, char **, int, bool *, bool *, bool *, int *, 
	real *, int *, int *, int *, int *, int *);

local void copyrighttext(bool);
local void copyrighttext3d(bool);

local int defpenwidth=1;
														// CHECK 2D --- OK!!!
void snap_anim_plplot_trajectory(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	char mode[4];

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

	if (scanopt(cmd.options, "save-to-file")) {
		strcpy(mode,"a");
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					PrintSnap(namebuf, snapnametmp, gd.nbody, cmd.options, mode);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		return;
	}

	if (scanopt(cmd.options, "save-animation")) {
		strcpy(mode,"w!");
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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
		return;
	}

	if (scanopt(cmd.options, "save")) {
		strcpy(mode,"w!");
		sprintf(namebuf1, cmd.out, snapcount);
		sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
		plsetopt("-dev","psc");
		plsfnam(namebuf);
		plsetopt("-ori",cmd.pl_ori);
		if (!strnull(cmd.pl_a))
			plsetopt("-a",cmd.pl_a);
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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
	} else {
		strcpy(mode,"w!");
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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


														// CHECK 2D --- OK!!!
#ifdef THREEDIM
void snap_anim_plplot_trajectory_3d(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;
	char mode[4];
	strcpy(mode,"w!");

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
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
				&gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			if (exist_snap) {
				++snapcount;
				Header_to_Global();
				PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
				animation_xy_trajectory_3d(snapname, snapcount);
				free_memory();
				free(bodytab);
				free(npltd.xval); free(npltd.yval); free(npltd.zval);
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
			++step;
	}

    plend();
}
#endif
														// CHECK 2D --- OK!!!
#ifdef THREEDIM
void snap_anim_plplot_3d(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	char mode[4];
	strcpy(mode,"w!");

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
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					plsetopt("-dev","psc");
					plsfnam(namebuf);
					plsetopt("-ori",cmd.pl_ori);
					plinit();
					animation_xyz_3d(snapname, snapcount);
					plend();
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.zval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					animation_xyz_3d(snapname, snapcount);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.zval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}
#endif

// AUN HAY PROBLEMAS EN EL CONSUMO DE MEMORIA...
void snap_anim_plplot(string snapname, string snapnametmp)// CHECK 2D --- OK!!!
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	char mode[4];
	strcpy(mode,"w!");

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
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
//                    plsetopt("-dev","psc");
                    plsetopt("-dev","png");
					plsfnam(namebuf);
					plsetopt("-ori",cmd.pl_ori);
					if (!strnull(cmd.pl_a))
						plsetopt("-a",cmd.pl_a);
					plinit();
                    if (scanopt(cmd.options, "type")) {
						animation_xy_type(snapname, snapcount);
                    } else {
                        if (scanopt(cmd.options, "type4b")) {
                            animation_xy_type_4b(snapname, snapcount);
                        } else {
                            animation_xy(snapname, snapcount);
                        }
                    }
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
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
                    if (scanopt(cmd.options, "type")) {
						animation_xy_type(snapname, snapcount);
                    } else {
                        if (scanopt(cmd.options, "type4b")) {
                            animation_xy_type_4b(snapname, snapcount);
                        } else {
                            animation_xy(snapname, snapcount);
                        }
                    }
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

// SLAB ANIM ...

#ifdef THREEDIM
void slab_anim_plplot(string snapname, string snapnametmp)
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

    strcpy(mode,"w!");

	if (scanopt(cmd.options, "save")) {
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
                          &exist_snap, &hdr, cmd.options, gd.model_comment);
				fprintf(stdout,"\n\ninputdata CPU time: %g\n",cputime()-cpustart);
				if (exist_snap) {
					++snapcount;
					PrintSnap_Slab(snapname, snapnametmp, gd.nbody, cmd.options, mode);
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
					PrintSnap_Slab(snapname, snapnametmp, gd.nbody, cmd.options, mode);
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

void slab_anim_plplot_trajectory(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	char mode[4];
    
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
    
	if (scanopt(cmd.options, "save-to-file")) {
		strcpy(mode,"a");
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
                          &gd.nbody, &ndim, &gd.tnow, 
                          &exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					PrintSnap_Slab(namebuf, snapnametmp, gd.nbody, cmd.options, mode);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		return;
	}
    
	if (scanopt(cmd.options, "save-animation")) {
		strcpy(mode,"w!");
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
                          &gd.nbody, &ndim, &gd.tnow, 
                          &exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap_Slab(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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
		return;
	}
    
	if (scanopt(cmd.options, "save")) {
		strcpy(mode,"w!");
		sprintf(namebuf1, cmd.out, snapcount);
		sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
		plsetopt("-dev","psc");
		plsfnam(namebuf);
		plsetopt("-ori",cmd.pl_ori);
		if (!strnull(cmd.pl_a))
			plsetopt("-a",cmd.pl_a);
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
                          &gd.nbody, &ndim, &gd.tnow, 
                          &exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap_Slab(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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
	} else {
		strcpy(mode,"w!");
		plinit();
        
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
                          &gd.nbody, &ndim, &gd.tnow, 
                          &exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintSnap_Slab(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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
#endif

// END SLAB ANIM ...

void bodies_anim_plplot(string snapname, string snapnametmp)// CHECK 2D --- OK!!!
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	char mode[4];

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
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
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
		strcpy(mode,"w!");
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_type(snapname, snapcount);
					else
						animation_xy(snapname, snapcount);
					free_memory();
					free(bodytab);
//					free(xval); free(yval);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}
														// CHECK 2D --- OK!!!
void bodies_anim_plplot_trajectory(string snapname, string snapnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	char mode[4];

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

	if (scanopt(cmd.options, "save-to-file")) {
		strcpy(mode,"a");
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
//					PrintSnap(namebuf, snapnametmp, gd.nbody, cmd.options, mode);
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
					free_memory();
					free(bodytab);
					free(npltd.xval); free(npltd.yval); free(npltd.Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		return;
	}

	if (scanopt(cmd.options, "save-animation")) {
		strcpy(mode,"w!");
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
//					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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
		return;
	}

	if (scanopt(cmd.options, "save")) {
		strcpy(mode,"w!");
		sprintf(namebuf1, cmd.out, snapcount);
		sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
		plsetopt("-dev","psc");
		plsfnam(namebuf);
		plsetopt("-ori",cmd.pl_ori);
		if (!strnull(cmd.pl_a))
			plsetopt("-a",cmd.pl_a);
		plinit();
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
//					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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
	} else {
		strcpy(mode,"w!");
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
//					PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options, mode);
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options, mode);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
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


/*
	if (scanopt(cmd.options, "save")) {
		strcpy(mode,"w!");
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options);
					sprintf(namebuf1, cmd.out, snapcount);
					sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
					plsetopt("-dev","psc");
					plsfnam(namebuf);
					plsetopt("-ori",cmd.pl_ori);
					if (!strnull(cmd.pl_a))
						plsetopt("-a",cmd.pl_a);
					plinit();
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
						animation_xy_trajectory(snapname, snapcount);
					plend();
					free_memory();
					free(bodytab);
					free(xval); free(yval); free(Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();

		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					PrintBodies(snapname, snapnametmp, gd.nbody, 
						gd.nbodiesID, &gd.bodyID[0], cmd.options);
					if (scanopt(cmd.options, "type"))
						animation_xy_trajectory_type(snapname, snapcount);
					else
						animation_xy_trajectory(snapname, snapcount);
					free_memory();
					free(bodytab);
					free(xval); free(yval); free(Typeval);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
*/
}

void rdf_anim_plplot(string fname, string fnametmp)		// CHECK 2D --- OK!!!
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
		gd.Rdf_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalRdf(fname, fnametmp);
					if (gd.Rdf_flag) {
						sprintf(namebuf1, cmd.out, snapcount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
						plinit();
//						animation_rdf("rdf.dat", snapcount);
						animation_rdf(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
						plend();
					}
					free_memory();
					free(bodytab);
					gd.Rdf_flag=0;
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		gd.Rdf_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalRdf(fname, fnametmp);
					if (gd.Rdf_flag) {
//						animation_rdf("rdf.dat", snapcount);
						animation_rdf(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
					}
					gd.Rdf_flag=0;
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

void vel_anim_plplot(string fname, string fnametmp)		// CHECK 2D --- OK!!!
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
		gd.Vel_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalVelDist(fname, fnametmp);
					if (gd.Vel_flag) {
						sprintf(namebuf1, cmd.out, snapcount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
						plinit();
//						animation_vel("vel.dat", snapcount);
						animation_vel(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
						plend();
					}
					free_memory();
					free(bodytab);
					gd.Vel_flag=0;
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		gd.Vel_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalVelDist(fname, fnametmp);
					if (gd.Vel_flag) {
//						animation_vel("vel.dat", snapcount);
						animation_vel(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
					}
					free_memory();
					free(bodytab);
					gd.Vel_flag=0;
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
		plend();
	}
}

void rhoaxes_anim_plplot(string fname, string fnametmp)	// CHECK 2D --- OK!!!
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

	plssub(1, NDIM);

	if (scanopt(cmd.options, "save")) {
		gd.RhoAxes_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalRhoAxes(fname, fnametmp);
					if (gd.RhoAxes_flag) {
						sprintf(namebuf1, cmd.out, snapcount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
						plinit();
						if (scanopt(cmd.options, "type"))
							animation_rhoaxes_type(fname, snapcount);
						else
							animation_rhoaxes(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
						plend();
					}
					gd.RhoAxes_flag=0;
					free_memory();
					free(bodytab);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		gd.RhoAxes_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalRhoAxes(fname, fnametmp);
					if (gd.RhoAxes_flag) {
						if (scanopt(cmd.options, "type"))
							animation_rhoaxes_type(fname, snapcount);
						else
							animation_rhoaxes(fname, snapcount);
						free(npltd.xval); free(npltd.yval);
					}
					gd.RhoAxes_flag=0;
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

void nfrecaxes_anim_plplot(string fname, string fnametmp)// CHECK 2D --- OK!!!
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;
	char fnamex[200], fnamey[200];
	char fnamex1[200], fnamey1[200];
	char fnamex2[200], fnamey2[200];
	char fnamexd[200], fnameyd[200];
#if (NDIM==3)
	char fnamez[200];
	char fnamez1[200];
	char fnamez2[200];
	char fnamezd[200];

	sprintf(fnamez, "%s-z", fname);
	sprintf(fnamez1, "%s-z1", fname);
	sprintf(fnamez2, "%s-z2", fname);
	sprintf(fnamezd, "%s-zd", fname);
#endif

	sprintf(fnamex, "%s-x", fname);
	sprintf(fnamey, "%s-y", fname);
	sprintf(fnamex1, "%s-x1", fname);
	sprintf(fnamey1, "%s-y1", fname);
	sprintf(fnamex2, "%s-x2", fname);
	sprintf(fnamey2, "%s-y2", fname);
	sprintf(fnamexd, "%s-xd", fname);
	sprintf(fnameyd, "%s-yd", fname);


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

	plssub(1, NDIM);

	if (scanopt(cmd.options, "save")) {
		gd.NFrecAxes_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalNFrecAxes(fname, fnametmp);
					if (gd.NFrecAxes_flag) {
						sprintf(namebuf1, cmd.out, snapcount);
						sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
						plsetopt("-dev","psc");
						plsfnam(namebuf);
						plsetopt("-ori",cmd.pl_ori);
						plinit();
//						animation_nfrecaxes("rhoaxes.dat", "nfrecy.dat", "nfrecz.dat", snapcount);
//						if (scanopt(cmd.options, "nfrec"))
//							animation_nfrecaxes(fnamex,fnamey,fnamez,snapcount);
#ifdef THREEDIM
						if (scanopt(cmd.options, "nfrec-1"))
							animation_nfrecaxes(fnamex1,fnamey1,fnamez1,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-2"))
							animation_nfrecaxes(fnamex2,fnamey2,fnamez2,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-d"))
							animation_nfrecaxes(fnamexd,fnameyd,fnamezd,snapcount);
						else
							animation_nfrecaxes(fnamex,fnamey,fnamez,snapcount);
#else
						if (scanopt(cmd.options, "nfrec-1"))
							animation_nfrecaxes(fnamex1,fnamey1,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-2"))
							animation_nfrecaxes(fnamex2,fnamey2,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-d"))
							animation_nfrecaxes(fnamexd,fnameyd,snapcount);
						else
							animation_nfrecaxes(fnamex,fnamey,snapcount);
#endif
						free(npltd.xval); free(npltd.yval);
						plend();
					}
					gd.NFrecAxes_flag=0;
					free_memory();
					free(bodytab);
				}
			} else
				if (step >= cmd.fsnap) moreSteps = 0;
			++step;
		}
	} else {
		plinit();
		gd.NFrecAxes_flag=0;
		while (moreSteps) {
			if (step >= cmd.isnap && step <= cmd.fsnap) {
				inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
					&gd.nbody, &ndim, &gd.tnow, 
					&exist_snap, &hdr, cmd.options, gd.model_comment);
				if (exist_snap) {
					++snapcount;
					Header_to_Global();
					EvalNFrecAxes(fname, fnametmp);
					if (gd.NFrecAxes_flag) {
//						animation_nfrecaxes("nfrecx.dat", "nfrecy.dat", "nfrecz.dat", snapcount);
//						if (scanopt(cmd.options, "nfrec"))
//							animation_nfrecaxes(fnamex,fnamey,fnamez,snapcount);
#ifdef THREEDIM
						if (scanopt(cmd.options, "nfrec-1"))
							animation_nfrecaxes(fnamex1,fnamey1,fnamez1,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-2"))
							animation_nfrecaxes(fnamex2,fnamey2,fnamez2,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-d"))
							animation_nfrecaxes(fnamexd,fnameyd,fnamezd,snapcount);
						else
							animation_nfrecaxes(fnamex,fnamey,fnamez,snapcount);
#else
						if (scanopt(cmd.options, "nfrec-1"))
							animation_nfrecaxes(fnamex1,fnamey1,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-2"))
							animation_nfrecaxes(fnamex2,fnamey2,snapcount);
						else
						if (scanopt(cmd.options, "nfrec-d"))
							animation_nfrecaxes(fnamexd,fnameyd,snapcount);
						else
							animation_nfrecaxes(fnamex,fnamey,snapcount);
#endif
						free(npltd.xval); free(npltd.yval);
					}
					gd.NFrecAxes_flag=0;
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

void general_fx_anim_plplot(string fname, string snapnametmp)// CHECK 2D --- OK!!!
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

local void animation_xy(string snapname, int snapcount)	// CHECK 2D --- OK!!!
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
//  int symbolcolor;
//  int symbolweight; 
  int symboltype;
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
	symboltype=4;

	if (cmd.withsymbols) {
		plcol0(cmd.symbolcolor);
		plssym(0.0, cmd.symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint, xnc[0], ync[0], symboltype);
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

local void animation_xy_type(string snapname, int snapcount)// CHECK 2D --- OK!!!
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
  int symbolcolor1, symbolcolor2;
//  int symbolweight;
  int symboltype;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  real viewportchr;

	printf("\n\nX-Y Animation type ...\n\n");

	col1=gd.column1;
	col2=gd.column2;
#ifdef THREEDIM
	col3=8;
#else
	col3=6;
#endif

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

	symbolcolor1=2;
	symbolcolor2=3;
//	if (scanopt(cmd.options, "save"))
//		symbolweight=2;
//	else
//		symbolweight=1;
	symboltype=4;

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
		plwid(cmd.symbolweight);
		plpoin(npoint1, x1, y1, symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor1);
		plssym(0.0, 0.0);
		plwid(cmd.symbolweight);
		plpoin(npoint1, x1, y1, 1);
		plwid(defpenwidth); 
	}

	if (cmd.withsymbols) {
		plcol0(symbolcolor2);
		plssym(0.0, cmd.symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint2, x2, y2, symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor2);
		plssym(0.0, 0.0);
		plwid(cmd.symbolweight);
		plpoin(npoint2, x2, y2, 1);
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

	printf("\n\n... done snap animation type plotting! ...\n\n");
}


// ======================================================
// SIRD:
local void animation_xy_type_4b(string snapname, int snapcount)// CHECK 2D --- OK!!!
{
#if !defined(BODY1)
#define BODY1    1
#endif
#if !defined(BODY2)
#define BODY2    2
#endif
    
#if !defined(BODY3)
#define BODY3    3
#endif
#if !defined(BODY4)
#define BODY4    4
#endif

  int i, col1, col2, col3, npoint;
  real *xnc[1], *ync[1];
  int *ival[1];
  real *x1, *y1, *x2, *y2;
  int npoint1, npoint2, i1, i2;
  real *x3, *y3, *x4, *y4;
  int npoint3, npoint4, i3, i4;
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
  int symbolcolor1, symbolcolor2;
  int symbolcolor3, symbolcolor4;
//  int symbolweight;
  int symboltype;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  real viewportchr;

    printf("\n\nX-Y Animation type ...\n\n");

    col1=gd.column1;
    col2=gd.column2;
#ifdef THREEDIM
    col3=8;
#else
    col3=6;
#endif

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
    npoint3=npoint4=0;
    for (i = 0; i < npoint; i++) {
        if (ival[0][i]==BODY1) ++npoint1;
        if (ival[0][i]==BODY2) ++npoint2;
        if (ival[0][i]==BODY3) ++npoint3;
        if (ival[0][i]==BODY4) ++npoint4;
    }
    fprintf(stdout,"\nReserving memory for %d %d %d %d points all types\n",
        npoint1, npoint2, npoint3, npoint4);
    x1 = (real *) allocate(npoint1 * sizeof(real));
    y1 = (real *) allocate(npoint1 * sizeof(real));
    x2 = (real *) allocate(npoint2 * sizeof(real));
    y2 = (real *) allocate(npoint2 * sizeof(real));
    x3 = (real *) allocate(npoint3 * sizeof(real));
    y3 = (real *) allocate(npoint3 * sizeof(real));
    x4 = (real *) allocate(npoint4 * sizeof(real));
    y4 = (real *) allocate(npoint4 * sizeof(real));
    i1=i2=0;
    i3=i4=0;
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
        if (ival[0][i]==BODY3) {
            x3[i3]=xnc[0][i];
            y3[i3]=ync[0][i];
            ++i3;
        }
        if (ival[0][i]==BODY4) {
            x4[i4]=xnc[0][i];
            y4[i4]=ync[0][i];
            ++i4;
        }
    }
    fprintf(stdout,"Total points assigned %d %d %d %d all types\n",
        i1, i2, i3, i4);
//#undef BODY1
//#undef BODY2

    printf("\nscaling done ...\n");

    axestype = 1;
//    if (scanopt(cmd.options, "save"))
//        axeswidth = 2;
//    else
//        axeswidth = 1;
    axescolor = 1;

    pllsty(axestype);
    plwid(cmd.axeswidth);
    plcol0(axescolor);

    pladv(1);
    plclear();
    viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
    plschr(0.0, viewportchr);        // Set the scaled size of the labels
    plvsta();
    plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);

    nlxdigmax=2;
    nlydigmax=2;
    nlcolor=1;
//    nlsize=1.0;

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
//    if (scanopt(cmd.options, "save"))
//        labelfontweight=2;
//    else
//        labelfontweight=1;

    plcol0(labelcolor);
    plschr(0.0, cmd.labelfontsize);
    plwid(cmd.labelfontweight);
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);
    plwid(defpenwidth);

    printf("\nsetting viewport ... done ...\n");

    symbolcolor1=2;                     // Yellow
    symbolcolor2=1;                     // Red
//    symbolcolor3=4;
//    symbolcolor3=6;                     // Como gray
    symbolcolor3=3;                     // Green
//    symbolcolor4=5;
    symbolcolor4=6;                     // Como gray
//    if (scanopt(cmd.options, "save"))
//        symbolweight=2;
//    else
//        symbolweight=1;
    symboltype=4;

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
        plwid(cmd.symbolweight);
        plpoin(npoint1, x1, y1, symboltype);
        plwid(defpenwidth);
    }

    if (cmd.withdots) {
        plcol0(symbolcolor1);
        plssym(0.0, 0.0);
        plwid(cmd.symbolweight);
        plpoin(npoint1, x1, y1, 1);
        plwid(defpenwidth);
    }

    if (cmd.withsymbols) {
        plcol0(symbolcolor2);
        plssym(0.0, cmd.symbolsize);
        plwid(cmd.symbolweight);
        plpoin(npoint2, x2, y2, symboltype);
        plwid(defpenwidth);
    }

    if (cmd.withdots) {
        plcol0(symbolcolor2);
        plssym(0.0, 0.0);
        plwid(cmd.symbolweight);
        plpoin(npoint2, x2, y2, 1);
        plwid(defpenwidth);
    }

    if (cmd.withsymbols) {
        plcol0(symbolcolor3);
        plssym(0.0, cmd.symbolsize);
        plwid(cmd.symbolweight);
        plpoin(npoint3, x3, y3, symboltype);
        plwid(defpenwidth);
    }

    if (cmd.withdots) {
        plcol0(symbolcolor3);
        plssym(0.0, 0.0);
        plwid(cmd.symbolweight);
        plpoin(npoint3, x3, y3, 1);
        plwid(defpenwidth);
    }

    if (cmd.withsymbols) {
        plcol0(symbolcolor4);
        plssym(0.0, cmd.symbolsize);
        plwid(cmd.symbolweight);
        plpoin(npoint4, x4, y4, symboltype);
        plwid(defpenwidth);
    }

    if (cmd.withdots) {
        plcol0(symbolcolor4);
        plssym(0.0, 0.0);
        plwid(cmd.symbolweight);
        plpoin(npoint4, x4, y4, 1);
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

    printf("\n\n... done snap animation type (4b) plotting! ...\n\n");
}
// ======================================================

														// CHECK 2D --- OK!!!
local void animation_xy_trajectory_type(string snapname, int snapcount)
{
  int i, col1, col2, col3, npoint;
  real *xnc[1], *ync[1];
  int *ival[1];
  real *x1, *y1, *x2, *y2;
  int npoint1, npoint2, i1, i2;
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
//  real labelfontsize;
//  string xlabel, ylabel, plotlabel;
  bool withsymbols, withdots;
  int symbolcolor;
//  int symbolweight;
  int symboltype;
  int symbolcolor1, symbolcolor2;
  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool f_name;
  PLINT cur_strm, new_strm;
  char   namebuf[200], namebuf1[200];
  real viewportchr;

	printf("\n\nX-Y Animation (Trajectories-type) ...\n\n");

	col1=gd.column1;
	col2=gd.column2;
#ifdef THREEDIM
	col3=8;
#else
	col3=6;
#endif

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
	fprintf(stdout,"Total points assigned %d %d both types\n",i1, i2);

	printf("\nscaling done ...\n");

	if (snapcount==1) {

	axestype = 1;
//	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);							// Line type (8 styles). See 'plstyl'
	plwid(cmd.axeswidth);							// Line width
    plcol0(axescolor);							// Set color of frame - viewport

	pladv(1);									// Set to zero to advance to the next subpage
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, viewportchr);		// Set the scaled size of the labels
	plvsta();									// Sets up a standard viewport, leaving apropriate margins
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);				// Sets up the world-coordinates windows

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
//	nlsize=1.0;

    plsxax(nlxdigmax, 0);						// To set digmax on x-axis
    plsyax(nlydigmax, 0);						// To set digmax on y-axis
	plcol0(nlcolor);							// Set color for numeric labels. Also set color for axes...
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
//	labelfontsize=1.;
//	labelfontweight=2;
//	xlabel="x";
//	ylabel="y";
//	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, cmd.labelfontsize);
	plwid(cmd.labelfontweight);
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	withsymbols=FALSE;
	withdots=TRUE;

	symbolcolor1=2;
	symbolcolor2=3;

	symbolsize=0.5;
//	symbolweight=1;
	symboltype=4;

	}

	if (cmd.withsymbols) {
		plcol0(symbolcolor1);
		plssym(0.0, cmd.symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint1, x1, y1, symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor1);
		plssym(0.0, 0.0);
		plwid(cmd.symbolweight);
		plpoin(npoint1, x1, y1, 1);
		plwid(defpenwidth); 
	}

	if (cmd.withsymbols) {
		plcol0(symbolcolor2);
		plssym(0.0, cmd.symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint2, x2, y2, symboltype);
		plwid(defpenwidth); 
	}

	if (cmd.withdots) {
		plcol0(symbolcolor2);
		plssym(0.0, 0.0);
		plwid(cmd.symbolweight);
		plpoin(npoint2, x2, y2, 1);
		plwid(defpenwidth); 
	}

	copyrighttext(TRUE);

	if (scanopt(cmd.options, "save-animation")) {
		sprintf(namebuf1, cmd.out, snapcount);
		sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
		printf("The current plot (trajectory type) was saved in color Postscript (pdf) under the name `%s'.\n", namebuf);
		plgstrm(&cur_strm);
		plmkstrm(&new_strm);
		plsfnam(namebuf);
//        plsdev("psc");
//        plsdev("png");
//        plsdev("jpg");              // not available
//        plsdev("pngcairo");
//        plsdev("ps");
//        plsdev("epscairo");
//        plsdev("pdf");
        plsdev("psttfc");
//        plsdev("xfig");
//
// Plotting Options:
// < 1> aqt        AquaTerm - Mac OS X
// < 2> xwin       X-Window (Xlib)
// < 3> ps         PostScript File (monochrome)
// < 4> psc        PostScript File (color)
// < 5> xfig       Xfig file
// < 6> null       Null device
// < 7> mem        User-supplied memory device
// < 8> psttf      PostScript File (monochrome)
// < 9> psttfc     PostScript File (color)
// <10> svg        Scalable Vector Graphics (SVG 1.1)
// <11> pdf        Portable Document Format PDF
// <12> xcairo     Cairo X Windows Driver
// <13> pdfcairo   Cairo PDF Driver
// <14> epscairo   Cairo EPS Driver
// <15> pscairo    Cairo PS Driver
// <16> svgcairo   Cairo SVG Driver
// <17> pngcairo   Cairo PNG Driver
// <18> memcairo   Cairo memory driver
// <19> extcairo   Cairo external context driver

		plcpstrm(cur_strm, 0);
		plreplot();
		plend1();
		plsstrm(cur_strm);
    }

	printf("\n\n... done snap trajectory type plotting! ...\n\n");
}
														// CHECK 2D --- OK!!!
local void animation_xy_trajectory(string snapname, int snapcount)
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
//  real labelfontsize;
//  string xlabel, ylabel, plotlabel;
  bool withsymbols, withdots;
  int symbolcolor;
//  int symbolweight;
  int symboltype;
  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool f_name;
  PLINT cur_strm, new_strm;
  char   namebuf[200], namebuf1[200];
  real viewportchr;

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
//	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);							// Line type (8 styles). See 'plstyl'
	plwid(cmd.axeswidth);							// Line width
    plcol0(axescolor);							// Set color of frame - viewport

	pladv(1);									// Set to zero to advance to the next subpage
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, viewportchr);		// Set the scaled size of the labels
	plvsta();									// Sets up a standard viewport, leaving apropriate margins
	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);				// Sets up the world-coordinates windows

	nlxdigmax=2;
	nlydigmax=2;
	nlcolor=1;
//	nlsize=1.0;

    plsxax(nlxdigmax, 0);						// To set digmax on x-axis
    plsyax(nlydigmax, 0);						// To set digmax on y-axis
	plcol0(nlcolor);							// Set color for numeric labels. Also set color for axes...
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
//	labelfontsize=1.;
//	labelfontweight=2;
//	xlabel="x";
//	ylabel="y";
//	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, cmd.labelfontsize);
	plwid(cmd.labelfontweight);
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

	withsymbols=FALSE;
	withdots=TRUE;
	symbolcolor=3;
	symbolsize=0.5;
//	symbolweight=1;
	symboltype=4;

	}

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint, xnc[0], ync[0], symboltype);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plwid(cmd.symbolweight);
		plpoin(npoint, xnc[0], ync[0], 1);
		plwid(defpenwidth); 
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

	if (scanopt(cmd.options, "save-animation")) {
		sprintf(namebuf1, cmd.out, snapcount);
		sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
		printf("The current plot (trajectory) was saved in color Postscript (png) under the name `%s'.\n", namebuf);
		plgstrm(&cur_strm);
		plmkstrm(&new_strm);
		plsfnam(namebuf);
//		plsdev("psc");
        plsdev("png");
		plcpstrm(cur_strm, 0);
		plreplot();
		plend1();
		plsstrm(cur_strm);
    }

	printf("\n\n... done snap trajectory plotting! ...\n\n");
}
														// CHECK 2D --- OK!!!
#ifdef THREEDIM
local void animation_xy_trajectory_3d(string snapname, int snapcount)
{
  int i, col1, col2, col3, npoint;
  real *xnc[1], *ync[1], *znc[1];
  bool x_autoscale, y_autoscale;
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
  bool withsymbols, withdots;
  int symbolcolor, symbolweight, symboltype;
  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool f_name;
  PLINT cur_strm, new_strm;
  char   namebuf[200], namebuf1[200];
  real alt, az, lsx, lsy, bx[2], by[2], bz[2];
  real viewportchr;

	printf("\n\nX-Y-Z Animation (Trajectories 3D) ...\n\n");

	col1=1;
	col2=2;
	col3=3;

	readin_pl_snap_3d(snapname, col1, col2, col3, &npoint);

	printf("\nnpoint = %d\n",npoint);

	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0]; znc[0]=&npltd.zval[0];

	if (snapcount==1) {
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
		if (gd.z_autoscale) {
			gd.zmin = znc[0][0];
			gd.zmax = znc[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.zmin>znc[0][i]) gd.zmin=znc[0][i];
				if (gd.zmax<znc[0][i]) gd.zmax=znc[0][i];
			}
		}
	}

	printf("\nscaling done ...\n");

	if (snapcount==1) {

	axestype = 1;
	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(axeswidth);
    plcol0(axescolor);

	lsx=0.22*rabs(gd.xmax-gd.xmin);
	lsy=0.22*rabs(gd.ymax-gd.ymin);

	pladv(1);

// Movidos aqui para que queden dentro del area de vision...
	nlsize=1.0;
	labelfontsize=1;
	viewportchr=MAX(labelfontsize,nlsize);
	plschr(0.0, 1.1*viewportchr);		// Set the scaled size of the labels

	plvsta();
	plwind(gd.xmin-1.95*lsx, gd.xmax+1.5*lsx, gd.ymin-0.5*lsy, gd.ymax+3.5*lsy);

	alt=20.0; az=60.0;
	plw3d(rabs(gd.xmax-gd.xmin), rabs(gd.ymax-gd.ymin), rabs(gd.zmax-gd.zmin), 
		gd.xmin, gd.xmax, gd.ymin, gd.ymax, gd.zmin, gd.zmax, 
		alt, az);

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


/*
	plbox3("bnstu", "x", 0.0, 0,
	       "bnstu", "y", 0.0, 0,
	       "bcdmnstuv", "z", 0.0, 0);
*/
	plbox3("bnstu", cmd.xlabel, 0.0, 0,
	       "bnstu", cmd.ylabel, 0.0, 0,
	       "bmnstu", "z", 0.0, 0);

	plwid(defpenwidth);

	printf("\nsetting viewport ... done ...\n");

	withsymbols=FALSE;
	withdots=TRUE;
	symbolcolor=3;
	symbolsize=0.5;
	symbolweight=1;
	symboltype=4;

	}

// Completing the box
	plwid(2);
	bx[0]=gd.xmin;
	by[0]=gd.ymin;
	bz[0]=gd.zmin;
	bx[1]=gd.xmin;
	by[1]=gd.ymin;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmin;
	by[0]=gd.ymin;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymin;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmin;
	by[0]=gd.ymin;
	bz[0]=gd.zmax;
	bx[1]=gd.xmin;
	by[1]=gd.ymax;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmin;
	by[0]=gd.ymax;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymax;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmax;
	by[0]=gd.ymax;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymin;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmax;
	by[0]=gd.ymin;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymin;
	bz[1]=gd.zmin;
	plline3(2,&bx[0],&by[0],&bz[0]);
	plwid(defpenwidth);
//

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, symbolsize);
		plwid(symbolweight);
		plpoin3(npoint, xnc[0], ync[0], znc[0], symboltype);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin3(npoint, xnc[0], ync[0], znc[0], 1);
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

	printf("\n\n... done snap trajectory 3d plotting! ...\n\n");
}
#endif
														// CHECK 2D --- OK!!!
#ifdef THREEDIM
local void animation_xyz_3d(string snapname, int snapcount)
{
  int i, col1, col2, col3, npoint;
  real *xnc[1], *ync[1], *znc[1];
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
  bool withsymbols, withdots;
  int symbolweight, symboltype;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool f_name;
  PLINT cur_strm, new_strm;
  char   namebuf[200], namebuf1[200];
  real alt, az, lsx, lsy, bx[2], by[2], bz[2];
  real viewportchr;


	printf("\n\nX-Y-Z Animation (Trajectories 3D) ...\n\n");

	col1=1;
	col2=2;
	col3=3;

	readin_pl_snap_3d(snapname, col1, col2, col3, &npoint);

	printf("\nnpoint = %d\n",npoint);

	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0]; znc[0]=&npltd.zval[0];

	if (snapcount==1) {
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
		if (gd.z_autoscale) {
			gd.zmin = znc[0][0];
			gd.zmax = znc[0][0]; 
			for (i = 1; i < npoint; i++) {
				if (gd.zmin>znc[0][i]) gd.zmin=znc[0][i];
				if (gd.zmax<znc[0][i]) gd.zmax=znc[0][i];
			}
		}
	}

	printf("\nscaling done ...\n");

//	if (snapcount==1) {

	axestype = 1;
	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(axeswidth);
    plcol0(axescolor);

	lsx=0.22*rabs(gd.xmax-gd.xmin);
	lsy=0.22*rabs(gd.ymax-gd.ymin);

	pladv(1);
	plclear();

// Movidos aqui para que queden dentro del area de vision...
	nlsize=1.0;
	labelfontsize=1;
	viewportchr=MAX(labelfontsize,nlsize);
	plschr(0.0, 1.1*viewportchr);		// Set the scaled size of the labels

	plvsta();
	plwind(gd.xmin-1.95*lsx, gd.xmax+1.5*lsx, gd.ymin-0.5*lsy, gd.ymax+3.5*lsy);

	alt=20.0; az=60.0;
	plw3d(rabs(gd.xmax-gd.xmin), rabs(gd.ymax-gd.ymin), rabs(gd.zmax-gd.zmin), 
		gd.xmin, gd.xmax, gd.ymin, gd.ymax, gd.zmin, gd.zmax, 
		alt, az);

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

	plbox3("bnstu", cmd.xlabel, 0.0, 0,
	       "bnstu", cmd.ylabel, 0.0, 0,
	       "bmnstu", "z", 0.0, 0);

	plwid(defpenwidth);

	printf("\nsetting viewport ... done ...\n");

	withsymbols=TRUE;
	withdots=FALSE;
	symbolweight=2;
	symboltype=4;

//	}

	if (withsymbols) {
		plcol0(cmd.symbolcolor);
		plssym(0.0, cmd.symbolsize);
		plwid(symbolweight);
		plpoin3(npoint, xnc[0], ync[0], znc[0], symboltype);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(cmd.symbolcolor);
		plssym(0.0, 0.0);
		plpoin3(npoint, xnc[0], ync[0], znc[0], 1);
	}

// Completing the box
	plcol0(axescolor);
	plwid(2);
	bx[0]=gd.xmin;
	by[0]=gd.ymin;
	bz[0]=gd.zmin;
	bx[1]=gd.xmin;
	by[1]=gd.ymin;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmin;
	by[0]=gd.ymin;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymin;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmin;
	by[0]=gd.ymin;
	bz[0]=gd.zmax;
	bx[1]=gd.xmin;
	by[1]=gd.ymax;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmin;
	by[0]=gd.ymax;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymax;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmax;
	by[0]=gd.ymax;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymin;
	bz[1]=gd.zmax;
	plline3(2,&bx[0],&by[0],&bz[0]);
	bx[0]=gd.xmax;
	by[0]=gd.ymin;
	bz[0]=gd.zmax;
	bx[1]=gd.xmax;
	by[1]=gd.ymin;
	bz[1]=gd.zmin;
	plline3(2,&bx[0],&by[0],&bz[0]);
	plwid(defpenwidth);
//

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

	copyrighttext3d(TRUE);

	printf("\n\n... done snap xyz 3d plotting! ...\n\n");
}
#endif

local void animation_rdf(string rdffilename, int snapcount)// CHECK 2D --- OK!!!
{
  int i, col1, col2, npoint, npoint_exp;
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
  bool plotjoined, withsymbols, withdots, withsymbols_exp;
  int linecolor, linewidth, linetype;
  int symbolcolor, symbolweight, symboltype;
  real symbolsize;
  int symbolcolor_exp, symbolweight_exp, symboltype_exp;
  real symbolsize_exp;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  string cprtext;
  int cprcolor, cprweight;
  real cprsize, cprdisp, cprpos, cprjust;
  string cprside;
  bool legends;
  int nlegends, legendspos;
  char *legendnames[2];
  bool lplotjoined[2], lwithsymbols[2], lwithdots[2];
  int lsymbolcolor[2], lsymboltype[2], llinecolor[2], llinetype[2], 
	lsymbolweight[2], llinewidth[2];
  real lsymbolsize[2];

	printf("\n\nRDF Animation ...\n\n");

	col1=1;
	col2=2;

	readin_pl_snap(rdffilename, col1, col2, &npoint);
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

	if (scanopt(cmd.options, "experiments")) {
		readin_pl_snap("gdr.exp", col1, col2, &npoint_exp);
		printf("\nnpoint_exp = %d\n",npoint_exp);
		xnc[1] = &npltd.xval[0]; ync[1]=&npltd.yval[0];
		if (gd.x_autoscale) {
			for (i = 0; i < npoint_exp; i++) {
				if (gd.xmin>xnc[1][i]) gd.xmin=xnc[1][i];
				if (gd.xmax<xnc[1][i]) gd.xmax=xnc[1][i];
			}
		}
		if (gd.y_autoscale) {
			for (i = 0; i < npoint_exp; i++) {
				if (gd.ymin>ync[1][i]) gd.ymin=ync[1][i];
				if (gd.ymax<ync[1][i]) gd.ymax=ync[1][i];
			}
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
	ylabel="Radial Distribution Function";
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

	withsymbols_exp=TRUE;

	linecolor=3;
	linewidth=3;
	linetype=1;

	symbolcolor=3;
	symbolsize=0.5;
	symbolweight=1;
	symboltype=4;

	symbolcolor_exp=2;
	symbolsize_exp=1.0;
	symbolweight_exp=2;
	symboltype_exp=4;

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

	if (scanopt(cmd.options, "experiments")) {
		if (withsymbols_exp) {
			plcol0(symbolcolor_exp);
			plssym(0.0, symbolsize_exp);
			plwid(symbolweight_exp);
			plpoin(npoint_exp, xnc[1], ync[1], symboltype_exp);
			plwid(defpenwidth); 
		}
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

	if (scanopt(cmd.options, "experiments")) {
		legends=TRUE;
		nlegends=2;
		legendnames[0] = (string) malloc(30);
		legendnames[1] = (string) malloc(30);
		legendnames[0]="Theory";
		legendnames[1]="Experiment";
		legendspos=1;
		lplotjoined[0]=TRUE;
		lplotjoined[1]=FALSE;
		lwithsymbols[0]=FALSE;
		lwithsymbols[1]=TRUE;
		lwithdots[0]=FALSE;
		lwithdots[1]=FALSE;
		lsymbolcolor[0]=3;
		lsymbolcolor[1]=2;
		lsymboltype[0]=1;
		lsymboltype[1]=4;
		lsymbolsize[0]=symbolsize;
		lsymbolsize[1]=symbolsize_exp;
		lsymbolweight[0]=symbolweight;
		lsymbolweight[1]=symbolweight_exp;
		llinecolor[0]=linecolor;
		llinecolor[1]=2;
		llinetype[0]=linetype;
		llinetype[1]=1;	
		llinewidth[0]=linewidth;
		llinewidth[1]=2;

		if (legends)
			setlegends(nlegends, &legendnames[0], legendspos,
				&lplotjoined[0], &lwithsymbols[0], &lwithdots[0], &lsymbolcolor[0], 
				&lsymbolsize[0], &lsymbolweight[0], &lsymboltype[0], &llinecolor[0],
				&llinetype[0], &llinewidth[0]);
	}

	copyrighttext(TRUE);

	printf("\n\n... done rdf plotting! ...\n\n");
}

local void animation_vel(string rdffilename, int snapcount)	// CHECK 2D --- OK!!!
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

	printf("\n\nVEL Animation ...\n\n");

	col1=1;
	col2=2;

	readin_pl_snap(rdffilename, col1, col2, &npoint);
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
	xlabel="v";
	ylabel="Velocity Distribution";
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

	printf("\n\n... done vel plotting! ...\n\n");
}
														// CHECK 2D --- OK!!!
local void animation_rhoaxes(string rhofilename, int snapcount)
{
  int i, col1, col2, col3, col4, npoint;
  real *xnc[2], *ync[2];
  real *znc[2], *wnc[2];
  int axestype; 
//  int axeswidth;
  int axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor;
//  int labelfontweight;
//  real labelfontsize;
  string xlabel, ylabel, plotlabel;
//  bool plotjoined;
  bool withsymbols, withdots; 
  int linecolor; 
//  int linewidth;
  int linetype;
  int symbolcolor; 
//  int symbolweight;
  int symboltype;
//  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool printcopyright;
  int ifile;

	printf("\n\nRHOAXES Animation ...\n\n");

for (ifile=0; ifile<NDIM; ifile++) {

	col1=(2*ifile+1)*2-1;
	col2=(2*ifile+1)*2;
//	col3=(2*ifile+1)*2+1;
//	col4=(2*ifile+1)*2+2;

	readin_pl_snap(rhofilename, col1, col2, &npoint);
	printf("\nnpoint = %d\n",npoint);
	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];	
//	znc[0] = &zval[0]; wnc[0]=&wval[0];	
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
//		for (i = 0; i < npoint; i++) {
//			if (gd.ymin>znc[0][i]) gd.ymin=znc[0][i];
//			if (gd.ymax<znc[0][i]) gd.ymax=znc[0][i];
//		}
//		for (i = 0; i < npoint; i++) {
//			if (gd.ymin>wnc[0][i]) gd.ymin=wnc[0][i];
//			if (gd.ymax<wnc[0][i]) gd.ymax=wnc[0][i];
//		}
	}

	printf("\nscaling done ...\n");

	axestype = 1;
//	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(cmd.axeswidth);
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
//	labelfontsize=1.2;
//	labelfontweight=2;
	if (ifile==0) xlabel="x";
	if (ifile==1) xlabel="y";
#ifdef THREEDIM
	if (ifile==2) xlabel="z";
#endif
	ylabel="Density profile";
	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, cmd.labelfontsize);
	plwid(cmd.labelfontweight);
    pllab(xlabel, ylabel, plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

//	plotjoined=TRUE;
	withsymbols=FALSE;
	withdots=FALSE;

	linecolor=3;
//	linewidth=3;
	linetype=1;

	symbolcolor=3;
//	symbolsize=0.5;
//	symbolweight=2;
	symboltype=4;

	if (cmd.plotjoined) {
		plcol0(linecolor);
		pllsty(linetype);
		plwid(cmd.linewidth);
		plline(npoint, xnc[0], ync[0]);
//		plcol0(linecolor+1);
//		pllsty(linetype+1);
//		plline(npoint, xnc[0], znc[0]);
//		plcol0(linecolor+2);
//		pllsty(linetype+2);
//		plline(npoint, xnc[0], wnc[0]);
		plwid(defpenwidth); 
	}

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, cmd.symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint, xnc[0], ync[0], symboltype);
//		plpoin(npoint, xnc[0], znc[0], symboltype+1);
//		plpoin(npoint, xnc[0], wnc[0], symboltype+2);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
//		plcol0(symbolcolor+1);
//		plpoin(npoint, xnc[0], znc[0], 1);
//		plcol0(symbolcolor+2);
//		plpoin(npoint, xnc[0], wnc[0], 1);
	}

	if (ifile==0) {
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
	}

//	printcopyright=TRUE;

//	copyrighttext(printcopyright);
}
	printcopyright=TRUE;
	copyrighttext(printcopyright);
	printf("\n\n... done rhoaxes plotting! ...\n\n");
}
														// CHECK 2D --- OK!!!
local void animation_rhoaxes_type(string rhofilename, int snapcount)
{
  int i, col1, col2, col3, col4, npoint;
  real *xnc[2], *ync[2];
  real *znc[2], *wnc[2];
  int axestype; 
//  int axeswidth;
  int axescolor;
  int nlxdigmax, nlydigmax, nlcolor;
  real nlsize;
  bool frame, xaxis, yaxis;
  string frame_xopt, frame_yopt;
  real frame_xtick, frame_ytick;
  int frame_nxsub, frame_nysub;
  real xorigin, yorigin;
  int labelcolor;
//  int labelfontweight;
//  real labelfontsize;
  string xlabel, ylabel, plotlabel;
//  bool plotjoined;
  bool withsymbols, withdots; 
  int linecolor; 
//  int linewidth;
  int linetype;
  int symbolcolor; 
//  int symbolweight;
  int symboltype;
//  real symbolsize;
  char timebuf[20];
  int timecolor, timeweight;
  real timesize, timedisp, timepos, timejust;
  string timeside;
  bool printcopyright;
  int ifile;

	printf("\n\nRHOAXES-TYPE Animation ...\n\n");

for (ifile=0; ifile<NDIM; ifile++) {

	col1=(2*ifile+1)*2-1;
	col2=(2*ifile+1)*2;
	col3=(2*ifile+1)*2+1;
	col4=(2*ifile+1)*2+2;

	readin_pl_snap_4c(rhofilename, col1, col2, col3, col4, &npoint);
	printf("\nnpoint = %d\n",npoint);
	xnc[0] = &npltd.xval[0]; ync[0]=&npltd.yval[0];	
	znc[0] = &npltd.zval[0]; wnc[0]=&npltd.wval[0];	
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
		for (i = 0; i < npoint; i++) {
			if (gd.ymin>znc[0][i]) gd.ymin=znc[0][i];
			if (gd.ymax<znc[0][i]) gd.ymax=znc[0][i];
		}
		for (i = 0; i < npoint; i++) {
			if (gd.ymin>wnc[0][i]) gd.ymin=wnc[0][i];
			if (gd.ymax<wnc[0][i]) gd.ymax=wnc[0][i];
		}
	}

	printf("\nscaling done ...\n");

	axestype = 1;
//	axeswidth = 2;
	axescolor = 1;

	pllsty(axestype);
	plwid(cmd.axeswidth);
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
//	labelfontsize=1.2;
//	labelfontweight=2;
	if (ifile==0) xlabel="x";
	if (ifile==1) xlabel="y";
#ifdef THREEDIM
	if (ifile==2) xlabel="z";
#endif
	ylabel="Density profile";
	plotlabel="";

    plcol0(labelcolor);
	plschr(0.0, cmd.labelfontsize);
	plwid(cmd.labelfontweight);
    pllab(xlabel, ylabel, plotlabel);
	plwid(defpenwidth); 

	printf("\nsetting viewport ... done ...\n");

//	plotjoined=TRUE;
	withsymbols=FALSE;
	withdots=FALSE;

	linecolor=3;
//	linewidth=3;
	linetype=1;

	symbolcolor=3;
//	symbolsize=0.5;
//	symbolweight=2;
	symboltype=4;

	if (cmd.plotjoined) {
		plcol0(linecolor);
		pllsty(linetype);
		plwid(cmd.linewidth);
		plline(npoint, xnc[0], ync[0]);
		plcol0(linecolor+1);
		pllsty(linetype+1);
		plline(npoint, xnc[0], znc[0]);
		plcol0(linecolor+2);
		pllsty(linetype+2);
		plline(npoint, xnc[0], wnc[0]);
		plwid(defpenwidth); 
	}

	if (withsymbols) {
		plcol0(symbolcolor);
		plssym(0.0, cmd.symbolsize);
		plwid(cmd.symbolweight);
		plpoin(npoint, xnc[0], ync[0], symboltype);
		plpoin(npoint, xnc[0], znc[0], symboltype+1);
		plpoin(npoint, xnc[0], wnc[0], symboltype+2);
		plwid(defpenwidth); 
	}

	if (withdots) {
		plcol0(symbolcolor);
		plssym(0.0, 0.0);
		plpoin(npoint, xnc[0], ync[0], 1);
		plcol0(symbolcolor+1);
		plpoin(npoint, xnc[0], znc[0], 1);
		plcol0(symbolcolor+2);
		plpoin(npoint, xnc[0], wnc[0], 1);
	}

	if (ifile==0) {
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
	}

//	printcopyright=TRUE;

//	copyrighttext(printcopyright);
}
	printcopyright=TRUE;
	copyrighttext(printcopyright);
	printf("\n\n... done rhoaxes-type plotting! ...\n\n");
}

														// CHECK 2D --- OK!!!
#ifdef THREEDIM
local void animation_nfrecaxes(string nfxfilename, string nfyfilename, 
		string nfzfilename, int snapcount)
#else
local void animation_nfrecaxes(string nfxfilename, string nfyfilename,
		int snapcount)
#endif
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
  int ifile;
#ifdef THREEDIM
  char *filename[3];
#else
  char *filename[2];
#endif

	filename[0] = (char *) malloc(20);
	filename[1] = (char *) malloc(20);
	strcpy(filename[0], nfxfilename);
	strcpy(filename[1], nfyfilename);
#ifdef THREEDIM
	filename[2] = (char *) malloc(20);
	strcpy(filename[2], nfzfilename);
#endif

	printf("\n\nNFRECAXES Animation ...\n\n");

for (ifile=0; ifile<NDIM; ifile++) {

//	col1=(ifile+1)*2-1;
//	col2=(ifile+1)*2;
	col1=1;
	col2=2;

//	readin_pl_snap(rhofilename, col1, col2, &npoint);
	readin_pl_snap(filename[ifile], col1, col2, &npoint);
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
	if (ifile==0) xlabel="#gr(x)";
	if (ifile==1) xlabel="rhoy";
#ifdef THREEDIM
	if (ifile==2) xlabel="rhoz";
#endif
	ylabel="Frequency";
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

	if (ifile==0) {
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
	}
}
	printcopyright=TRUE;
	copyrighttext(printcopyright);
	printf("\n\n... done nfrecaxes plotting! ...\n\n");
}

														// CHECK 2D --- OK!!!
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
  int symbolcolor, symbolweight, symboltype;
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
	symboltype=4;

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
		plpoin(npoint, xnc[0], ync[0], symboltype);
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
														// CHECK 2D --- OK!!!
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

local void copyrighttext(bool cpr)						// CHECK 2D --- OK!!!
{
  string cprtext;
  int cprcolor, cprweight;
  real cprsize, cprdisp, cprpos, cprjust;
  string cprside;

	if (cpr) {
		plfontld(1);
		cprtext="#(274)2005-2009 Mar";
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

local void copyrighttext3d(bool cpr)						// CHECK 2D --- OK!!!
{
  string cprtext;
  int cprcolor, cprweight;
  real cprsize, cprdisp, cprpos, cprjust;
  string cprside;

	if (cpr) {
		plfontld(1);
		cprtext="#(274)2005-2009 Mar";
		cprcolor=8;
		cprsize=0.5;
		cprweight=1;
		cprside="b";
//		cprdisp=8;
//		cprpos=0.9;
		cprdisp=3;
		cprpos=0.8;
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

