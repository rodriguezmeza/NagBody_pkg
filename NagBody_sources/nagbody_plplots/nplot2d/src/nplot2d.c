/* =============================================================================
	MODULE: nplot2d.c				[nplot2d]
	Written by: M.A. Rodriguez-Meza
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

	Major revisions: November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/

#include "globaldefs.h"


static PLGraphicsIn gin;

/*
 * You can select a different set of symbols to use when plotting the
 * lines by changing the value of OFFSET. 
 */

//#define OFFSET  2

#define NORMAL		0
#define NORMALLOG	1
#define LOGNORMAL	2
#define LOGLOG		3

local int defpenwidth=1;

local void setploptions(void);

local void plot2d_normal(void);
local void plot2d_normallog(void);
local void plot2d_lognormal(void);
local void plot2d_loglog(void);
local void plot2d_graphicsarray(void);

local void plottext(void);
local void setlegends(int);

void Plot2DDriver(void)
{
    PLINT digmax;				// To manipulate numeric label representation
								// not active by now
	if (!(gd.graphicsarray_nrows==1 && gd.graphicsarray_ncols==1)) {
		printf("\n\nWe subdivide the page ...\n");
		plssub(gd.graphicsarray_ncols, gd.graphicsarray_nrows);
	}

	setploptions();

    plinit();					// Initializes plplot lib

    if (cmd.fontset) {
		plfontld(1);
		plfont(cmd.fontchr);	// Set character font (Normal, Roman, Italic, Script)
	} else
		plfontld(0);

	if (gd.graphicsarray_nrows==1 && gd.graphicsarray_ncols==1)
		switch (cmd.plottype) {
			case NORMAL:
				plot2d_normal(); break;
			case NORMALLOG:
				plot2d_normallog(); break;
			case LOGNORMAL:
				plot2d_lognormal(); break;
			case LOGLOG:
				plot2d_loglog(); break;
		}
	else
		plot2d_graphicsarray();

// Let's get some user input (locate mode)

    if (cmd.locatemode) {
	for (;;) {
	  if (! plGetCursor(&gin)) break;
	  if (gin.keysym == PLK_Escape) break;

	    pltext();
	    if (gin.keysym < 0xFF && isprint(gin.keysym)) 
		printf("subwin = %d, wx = %f,  wy = %f, dx = %f,  dy = %f,  c = '%c'\n",
		       gin.subwindow, gin.wX, gin.wY, gin.dX, gin.dY, gin.keysym);
	    else
		printf("subwin = %d, wx = %f,  wy = %f, dx = %f,  dy = %f,  c = 0x%02x\n",
		       gin.subwindow, gin.wX, gin.wY, gin.dX, gin.dY, gin.keysym);

	    plgra();
	}
    }

    plend();
}

#undef NORMAL
#undef NORMALLOG
#undef LOGNORMAL
#undef LOGLOG


//---------MENU COMMAND LINE PLPLOT OPTIONS------------------------------------
local void setploptions(void) {

// Option set to false == use default

	if (cmd.pl_showall)							// Si funciona
		plsetopt("-showall","");

	if (cmd.pl_h)								// No funciona (pero puede ser por que
		plsetopt("-h","");						// esta activado el mode_quiet)

	if (cmd.pl_v)								// No funciona (pero puede ser por que
		plsetopt("-v","");						// esta activado el mode_quiet, ver plargs.c)

	if (cmd.pl_verbose)							// Si funciona
		plsetopt("-verbose","");

	if (cmd.pl_debug)							// Si funciona
		plsetopt("-debug","");

	if (cmd.pl_hack)							// Si funciona 
		plsetopt("-hack","");

	if (!strnull(cmd.pl_dev))
		plsetopt("-dev",cmd.pl_dev);		// It can be used also 'plsdev(cmd.pl_dev)'

	if (!strnull(cmd.pl_o))
		plsetopt("-o",cmd.pl_o);

	if (!strnull(cmd.pl_display))
		plsetopt("-display",cmd.pl_display);

//	if (!strnull(pl_px))
//		plsetopt("-px",pl_px);

//	if (!strnull(pl_py))
//		plsetopt("-py",pl_py);

	if (!strnull(cmd.pl_geometry))
		plsetopt("-geometry",cmd.pl_geometry);

	if (!strnull(cmd.pl_wplt))
		plsetopt("-wplt",cmd.pl_wplt);

	if (!strnull(cmd.pl_mar))
		plsetopt("-mar",cmd.pl_mar);

	if (!strnull(cmd.pl_a))
		plsetopt("-a",cmd.pl_a);

	if (!strnull(cmd.pl_jx))
		plsetopt("-jx",cmd.pl_jx);

	if (!strnull(cmd.pl_jy))
		plsetopt("-jy",cmd.pl_jy);

	if (!strnull(cmd.pl_ori))
		plsetopt("-ori",cmd.pl_ori);

	if (cmd.pl_freeaspect)								// Si funciona 
		plsetopt("-freeaspect","");

	if (cmd.pl_portrait)								// Si funciona ... dudo
		plsetopt("-portrait","");

	if (!strnull(cmd.pl_width))							// Si funciona ... dudo
		plsetopt("-width",cmd.pl_width);

	if (!strnull(cmd.pl_bg))							// Set background color
		plsetopt("-bg",cmd.pl_bg);

	if (!strnull(cmd.pl_ncol0))
		plsetopt("-ncol0",cmd.pl_ncol0);

	if (!strnull(cmd.pl_ncol1))
		plsetopt("-ncol1",cmd.pl_ncol1);

	if (cmd.pl_fam)								// funciona?
		plsetopt("-fam","");

	if (!strnull(cmd.pl_fsiz))							// funciona?
		plsetopt("-fsiz",cmd.pl_fsiz);

	if (!strnull(cmd.pl_fbeg))							// funciona?
		plsetopt("-fbeg",cmd.pl_fbeg);

	if (!strnull(cmd.pl_finc))							// funciona?
		plsetopt("-finc",cmd.pl_finc);

	if (!strnull(cmd.pl_fflen))							// funciona?
		plsetopt("-fflen",cmd.pl_fflen);

	if (cmd.pl_nopixmap)								// funciona?
		plsetopt("-nopixmap","");

	if (cmd.pl_db)								// funciona?
		plsetopt("-db","");

	if (cmd.pl_np)								// funciona?
		plsetopt("-np","");

	if (!strnull(cmd.pl_bufmax))							// funciona?
		plsetopt("-bufmax",cmd.pl_bufmax);

//	if (!strnull(cmd.pl_server_name))							// funciona?
//		plsetopt("-server_name",cmd.pl_server_name);

//	if (!strnull(cmd.pl_plserver))								// funciona?
//		plsetopt("-plserver",cmd.pl_plserver);

//	if (!strnull(cmd.pl_plwindow))								// funciona?
//		plsetopt("-plwindow",cmd.pl_plwindow);

//	if (!strnull(cmd.pl_tcl_cmd))					// Depreciated - use -drvopt tcl_cmd= ...
//		plsetopt("-tcl_cmd",cmd.pl_tcl_cmd);

//	if (!strnull(cmd.pl_auto_path))								// funciona?
//		plsetopt("-auto_path",cmd.pl_auto_path);

//	if (!strnull(cmd.pl_tk_file))								// funciona?
//		plsetopt("-tk_file",cmd.pl_tk_file);

	if (!strnull(cmd.pl_dpi))								// funciona?
		plsetopt("-dpi",cmd.pl_dpi);

	if (!strnull(cmd.pl_compression))								// funciona?
		plsetopt("-compression",cmd.pl_compression);

	if (!strnull(cmd.pl_drvopt))								// funciona?
		plsetopt("-drvopt",cmd.pl_drvopt);

//-----------------------------------------------------------------------------

}

//-------------------------------- COMIENZA PLOT2D NORMAL ------------------------------------

local void plot2d_normal(void)
{
  int i, ifile;
  real *xnc[gd.nfiles], *ync[gd.nfiles], *yminnc[gd.nfiles], *ymaxnc[gd.nfiles];
  real *yminnctmp[gd.nfiles];
  int vnpoint[gd.nfiles];
  bool flag;
  real viewportchr;

	printf("\n\nNormal plotting ...\n\n");
// COMIENZA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	for (ifile=0; ifile<gd.nfiles; ifile++) {

		if (strnull(cmd.inputfile))
			TestData(&vnpoint[ifile]);
		else
			if (!strnull(cmd.witherrorbars) && gd.errorbars[ifile]==1)
				if (cmd.errorbarstype==1) {
					InputData_4c(gd.filenames[ifile], gd.vcol1[ifile], gd.vcol2[ifile],
						gd.vcol3[ifile], gd.vcol4[ifile], &vnpoint[ifile]);
				} else {
					InputData_3c(gd.filenames[ifile], gd.vcol1[ifile], gd.vcol2[ifile],
						gd.vcol3[ifile], &vnpoint[ifile]);
				}
			else
				InputData(gd.filenames[ifile], gd.vcol1[ifile], gd.vcol2[ifile], &vnpoint[ifile]);

		printf("\nvnpoint[%d] = %d\n",ifile,vnpoint[ifile]);

		gd.npoint=vnpoint[ifile];

		xnc[ifile] = &gd.xval[0]; ync[ifile]=&gd.yval[0];
		if (!strnull(cmd.witherrorbars) && gd.errorbars[ifile]==1) {
			if (cmd.errorbarstype==1) {
				yminnc[ifile]=&gd.yminval[0];
				ymaxnc[ifile]=&gd.ymaxval[0];
			} else {
				yminnctmp[ifile]=&gd.yminval[0];
				yminnc[ifile]=(real *)allocate(vnpoint[ifile]*sizeof(real));
				ymaxnc[ifile]=(real *)allocate(vnpoint[ifile]*sizeof(real));
			}
		}

		if (!strnull(cmd.witherrorbars) && gd.errorbars[ifile]==1) {
			if (cmd.errorbarstype==0) {
				for (i = 0; i < vnpoint[ifile]; i++) {
					yminnc[ifile][i] = ync[ifile][i] - yminnctmp[ifile][i];
					ymaxnc[ifile][i] = ync[ifile][i] + yminnctmp[ifile][i];
				}
			}
		}

		if (gd.x_autoscale) {
			if (ifile==0) {
				gd.xmin = xnc[ifile][0];
				gd.xmax = xnc[ifile][0]; 
				for (i = 1; i < vnpoint[ifile]; i++) {
					if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
					if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
				}
			} else
				for (i = 0; i < vnpoint[ifile]; i++) {
					if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
					if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
				}
		}

		if (gd.y_autoscale) {
			if (ifile==0) {
				gd.ymin = ync[ifile][0];
				gd.ymax = ync[ifile][0]; 
				for (i = 1; i < vnpoint[ifile]; i++) {
					if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
					if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
				}
			} else
				for (i = 0; i < vnpoint[ifile]; i++) {
					if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
					if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
				}
	
			if (!strnull(cmd.witherrorbars) && gd.errorbars[ifile]==1)
				for (i = 0; i < vnpoint[ifile]; i++) {
					if (gd.ymin>yminnc[ifile][i]) gd.ymin=yminnc[ifile][i];
					if (gd.ymax<ymaxnc[ifile][i]) gd.ymax=ymaxnc[ifile][i];
				}
		}

	}

printf("\nscaling done ...\n");

// TERMINA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	pllsty(cmd.axestype);								// Line type (8 styles). See 'plstyl'
	plwid(cmd.axeswidth);								// Line width
    plcol0(cmd.axescolor);								// Set color of frame - viewport

//    plenv(xmin, xmax, ymin, ymax, 0, 0);			// Set viewport. Axes scaled separately. Box drawn.

// The following commands are equivalent to the previous one
	pladv(0);							// Set to zero to advance to the next subpage

// Movidos aqui para que queden dentro del area de vision...
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, 1.1*viewportchr);		// Set the scaled size of the labels
	plvsta();						// Sets up a standard viewport, leaving apropriate margins

	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);		// Sets up the world-coordinates windows

    plsxax(cmd.nlxdigmax, 0);						// To set digmax on x-axis
    plsyax(cmd.nlydigmax, 0);						// To set digmax on y-axis
	plcol0(cmd.nlcolor);							// Set color for numeric labels. Also set color for axes...
	plschr(0.0,cmd.nlsize);

	if (strnull(cmd.framestyle)) {

	if (cmd.frame)
		if (!cmd.gridlines)
			plbox("bcnst", 0.0, 0, "bcnstv", 0.0, 0);		// Frame is drawn, axis, ticks, labels
		else
			plbox("bcgnst", 0.0, 0, "bcgnstv", 0.0, 0);		// Frame is drawn, axis, ticks, labels
	else
		plbox("bnst", 0.0, 0, "bnstv", 0.0, 0);				// Draw only x and y axis not a frame

	if (cmd.xaxis)
		plaxes(gd.xorigin, gd.yorigin, "anst", 0.0, 0, "", 0.0, 0);// Draw x-axis (y=yorigin)
	if (cmd.yaxis)
		plaxes(gd.xorigin, gd.yorigin, "", 0.0, 0, "anstv", 0.0, 0);// Draw y-axis (x=xorigin)

	} else {
		if (cmd.frame)
			plbox(gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
		if (cmd.xaxis)
			plaxes(gd.xorigin,gd.yorigin,gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,"",gd.frame_ytick,gd.frame_nysub);
		if (cmd.yaxis)
			plaxes(gd.xorigin,gd.yorigin,"",gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
	}

	plwid(defpenwidth);					// Return to default pen width 

    plcol0(cmd.labelcolor);					// Set labels color
	plschr(0.0, cmd.labelfontsize);			// Set the scaled size of the labels
	plwid(cmd.labelfontweight);				// Set font weight
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);	// text labels for b, l hand and t of the viewport
	plwid(defpenwidth);					// Return to default pen width 

	printf("\nsetting viewport ... done ...\n");

	for (ifile=0; ifile<gd.nfiles; ifile++) {
//		if (cmd.plotjoined) {
//	printf("\n\nAqui voy %d %d\n",gd.plotjoined[ifile],cmd.linetype+ifile);
		if (gd.plotjoined[ifile]) {
			plcol0(cmd.linecolor+ifile);						// Line color
//			pllsty(cmd.linetype+ifile);							// Line type (8 styles). See 'plstyl'
			pllsty(gd.linetype[ifile]);							// Line type (8 styles). See 'plstyl'
			plwid(cmd.linewidth);								// Line width
			plline(vnpoint[ifile], xnc[ifile], ync[ifile]);	// Draw line through data points
			plwid(defpenwidth);								// Return to default pen width 
		}
		if (!strnull(cmd.witherrorbars) && gd.errorbars[ifile]==1) {
			plwid(cmd.symbolweight);							// Symbol weight
			pllsty(1);							// Line type (8 styles). See 'plstyl'
			plerry(vnpoint[ifile], xnc[ifile], yminnc[ifile], ymaxnc[ifile]);
			plwid(defpenwidth);								// Return to default pen width 
		}
//		if (cmd.withsymbols) {
		if (gd.withsymbols[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plssym(0.0, cmd.symbolsize);
			plwid(cmd.symbolweight);							// Symbol weight
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
			plwid(defpenwidth);								// Return to default pen width 
		}
//		if (cmd.withdots) {
		if (gd.withdots[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plssym(0.0, 0.0);								// Set symbol to pixel
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
		}
	}

	if (!strnull(cmd.legends))
		setlegends(gd.nfiles);
	plottext();
	printf("\n\n... done plotting! ...\n\n");

}

//-------------------------------- TERMINA PLOT2D NORMAL ------------------------------------

//-------------------------------- COMIENZA PLOT2D NORMALLOG ------------------------------------

local void plot2d_normallog(void)
{
//  real *xn, *yn;
//  pointptr p;
  int i, ifile;
  real *xnc[gd.nfiles], *ync[gd.nfiles];
  real tmp;
  int vnpoint[gd.nfiles];
  real viewportchr;

	printf("\n\nNormal-Log plotting ...\n\n");
// COMIENZA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	for (ifile=0; ifile<gd.nfiles; ifile++) {

	if (strnull(cmd.inputfile))
		TestData(&vnpoint[ifile]);
	else
		InputData(gd.filenames[ifile],gd.vcol1[ifile], gd.vcol2[ifile], &vnpoint[ifile]);

//	InputData(filenames[ifile],vcol1[ifile], vcol2[ifile], &vnpoint[ifile]);

	printf("\nvnpoint[%d] = %d\n",ifile,vnpoint[ifile]);

	gd.npoint=vnpoint[ifile];
/*
	xn = (realptr) allocate(vnpoint[ifile]*sizeof(real));
	yn = (realptr) allocate(vnpoint[ifile]*sizeof(real));

	DO_POINT(p, pointtab, pointtab+vnpoint[ifile]) {
		i = p-pointtab;
		xn[i]=Xval(p);
		yn[i]=Yval(p);
	}
*/

	xnc[ifile] = &gd.xval[0]; ync[ifile]=&gd.yval[0];

//	DO_POINT(p, gd.pointtab, gd.pointtab+vnpoint[ifile]) {
	for (i=0; i<vnpoint[ifile]; ++i) {
//		i = p-gd.pointtab;
		tmp = rlog10(ync[ifile][i]);
		ync[ifile][i] = tmp;
	}

	if (gd.x_autoscale) {
		if (ifile==0) {
			gd.xmin = xnc[ifile][0];
			gd.xmax = xnc[ifile][0]; 
			for (i = 1; i < vnpoint[ifile]; i++) {
				if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
				if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
			}
		} else
			for (i = 0; i < vnpoint[ifile]; i++) {
				if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
				if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
			}
	}

	if (gd.y_autoscale) {
		if (ifile==0) {
			gd.ymin = ync[ifile][0];
			gd.ymax = ync[ifile][0]; 
			for (i = 1; i < vnpoint[ifile]; i++) {
				if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
				if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
			}
		} else
			for (i = 0; i < vnpoint[ifile]; i++) {
				if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
				if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
			}
	}
	
	}

	printf("\nscaling done ...\n");

// TERMINA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	pllsty(cmd.axestype);								// Line type (8 styles). See 'plstyl'
	plwid(cmd.axeswidth);								// Line width
    plcol0(cmd.axescolor);								// Set color of frame - viewport

	pladv(0);										// Set to zero to advance to the next subpage

// Movidos aqui para que queden dentro del area de vision...
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, 1.1*viewportchr);		// Set the scaled size of the labels
	plvsta();										// Sets up a standard viewport, leaving apropriate margins

	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);					// Sets up the world-coordinates windows

    plsxax(cmd.nlxdigmax, 0);						// To set digmax on x-axis
    plsyax(cmd.nlydigmax, 0);						// To set digmax on y-axis
	plcol0(cmd.nlcolor);							// Set color for numeric labels. Also set color for axes...
	plschr(0.0,cmd.nlsize);

	if (strnull(cmd.framestyle)) {

	if (cmd.frame)
		if (!cmd.gridlines)
			plbox("bcnst", 0.0, 0, "bcnstvl", 0.0, 0);	// Specify whether a frame is drawn, axis, ticks, labels
		else
			plbox("bcgnst", 0.0, 0, "bcgnstvl", 0.0, 0);	// Specify whether a frame is drawn, axis, ticks, labels
	else
		plbox("bnst", 0.0, 0, "bnstvl", 0.0, 0);		// Draw only x and y axis not a frame

	if (cmd.xaxis)
		plaxes(gd.xorigin, gd.yorigin, "anst", 0.0, 0, "", 0.0, 0);// Draw x-axis (y=yorigin)
	if (cmd.yaxis)
		plaxes(gd.xorigin, gd.yorigin, "", 0.0, 0, "anstv", 0.0, 0);// Draw y-axis (x=xorigin). 
															// In log scale can not be done when y=0
	} else {
		strcat(gd.frame_yopt, "l");
		if (cmd.frame)
			plbox(gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
		if (cmd.xaxis)
			plaxes(gd.xorigin,gd.yorigin,gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,"",gd.frame_ytick,gd.frame_nysub);
		if (cmd.yaxis)
			plaxes(gd.xorigin,gd.yorigin,"",gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
	}

	plwid(defpenwidth);							// Return to default pen width 

    plcol0(cmd.labelcolor);							// Set labels color
	plschr(0.0, cmd.labelfontsize);					// Set the scaled size of the labels
	plwid(cmd.labelfontweight);									// Set font weight
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);				// text labels for bottom, left hand and top of the viewport
	plwid(defpenwidth);							// Return to default pen width 

printf("\nsetting viewport ... done ...\n");

	for (ifile=0; ifile<gd.nfiles; ifile++) {
		if (gd.plotjoined[ifile]) {
			plcol0(cmd.linecolor+ifile);								// Line color
//			pllsty(cmd.linetype+ifile);							// Line type (8 styles). See 'plstyl'
			pllsty(gd.linetype[ifile]);							// Line type (8 styles). See 'plstyl'
			plwid(cmd.linewidth);							// Line width
			plline(vnpoint[ifile], xnc[ifile], ync[ifile]);		// Draw line through data points
			plwid(defpenwidth);							// Return to default pen width 
		}
//		if (cmd.withsymbols) {
		if (gd.withsymbols[ifile]) {
			plcol0(cmd.symbolcolor+ifile);								// symbol color
			plssym(0.0, cmd.symbolsize);
			plwid(cmd.symbolweight);							// Symbol weight
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
			plwid(defpenwidth);							// Return to default pen width 
		}
//		if (cmd.withdots) {
		if (gd.withdots[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plssym(0.0, 0.0);								// Set symbol to pixel
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
		}
	}
	
	if (!strnull(cmd.legends))
		setlegends(gd.nfiles);
	plottext();
	printf("\n\n... done normal-log plotting! ...\n\n");	
}

//-------------------------------- TERMINA PLOT2D NORMALLOG ------------------------------------

//-------------------------------- COMIENZA PLOT2D LOGNORMAL ------------------------------------

local void plot2d_lognormal(void)
{
//  real *xn, *yn;
//  pointptr p;
  int i, ifile;
  real *xnc[gd.nfiles], *ync[gd.nfiles];
  int vnpoint[gd.nfiles];
  real tmp;
  real viewportchr;

	printf("\n\nLog-Normal plotting ...\n\n");
// COMIENZA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	for (ifile=0; ifile<gd.nfiles; ifile++) {

	if (strnull(cmd.inputfile))
		TestData(&vnpoint[ifile]);
	else
		InputData(gd.filenames[ifile],gd.vcol1[ifile], gd.vcol2[ifile], &vnpoint[ifile]);

//	InputData(filenames[ifile],vcol1[ifile], vcol2[ifile], &vnpoint[ifile]);

	printf("\nvnpoint[%d] = %d\n",ifile,vnpoint[ifile]);

	gd.npoint=vnpoint[ifile];
/*
	xn = (realptr) allocate(vnpoint[ifile]*sizeof(real));
	yn = (realptr) allocate(vnpoint[ifile]*sizeof(real));

	DO_POINT(p, pointtab, pointtab+vnpoint[ifile]) {
		i = p-pointtab;
		xn[i]=Xval(p);
		yn[i]=Yval(p);
	}
*/

	xnc[ifile] = &gd.xval[0]; ync[ifile]=&gd.yval[0];

//	DO_POINT(p, gd.pointtab, gd.pointtab+vnpoint[ifile]) {
	for (i=0; i<vnpoint[ifile]; ++i) {
//		i = p-gd.pointtab;
		tmp = rlog10(xnc[ifile][i]);
		xnc[ifile][i] = tmp;
	}

	if (gd.x_autoscale) {
		if (ifile==0) {
			gd.xmin = xnc[ifile][0];
			gd.xmax = xnc[ifile][0]; 
			for (i = 1; i < vnpoint[ifile]; i++) {
				if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
				if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
			}
		} else
			for (i = 0; i < vnpoint[ifile]; i++) {
				if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
				if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
			}
	}

	if (gd.y_autoscale) {
		if (ifile==0) {
			gd.ymin = ync[ifile][0];
			gd.ymax = ync[ifile][0]; 
			for (i = 1; i < vnpoint[ifile]; i++) {
				if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
				if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
			}
		} else
			for (i = 0; i < vnpoint[ifile]; i++) {
				if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
				if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
			}
	}
	
	}

	printf("\nscaling done ...\n");

// TERMINA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	pllsty(cmd.axestype);								// Line type (8 styles). See 'plstyl'
	plwid(cmd.axeswidth);								// Line width
    plcol0(cmd.axescolor);								// Set color of frame - viewport

	pladv(0);										// Set to zero to advance to the next subpage

// Movidos aqui para que queden dentro del area de vision...
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, 1.1*viewportchr);		// Set the scaled size of the labels
	plvsta();										// Sets up a standard viewport, leaving apropriate margins

	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);					// Sets up the world-coordinates windows

    plsxax(cmd.nlxdigmax, 0);						// To set digmax on x-axis
    plsyax(cmd.nlydigmax, 0);						// To set digmax on y-axis
	plcol0(cmd.nlcolor);							// Set color for numeric labels. Also set color for axes...
	plschr(0.0,cmd.nlsize);

	if (strnull(cmd.framestyle)) {

	if (cmd.frame)
		if (!cmd.gridlines)
			plbox("bcnstl", 0.0, 0, "bcnstv", 0.0, 0);	// Specify whether a frame is drawn, axis, ticks, labels
		else
			plbox("bcgnstl", 0.0, 0, "bcgnstv", 0.0, 0);	// Specify whether a frame is drawn, axis, ticks, labels
	else
		plbox("bnstl", 0.0, 0, "bnstv", 0.0, 0);		// Draw only x and y axis not a frame

	if (cmd.xaxis)
		plaxes(gd.xorigin, gd.yorigin, "anst", 0.0, 0, "", 0.0, 0);// Draw x-axis (y=yorigin)
														// In log scale can not be done when y=0
	if (cmd.yaxis)
		plaxes(gd.xorigin, gd.yorigin, "", 0.0, 0, "anstv", 0.0, 0);// Draw y-axis (x=xorigin). 

	} else {
		strcat(gd.frame_xopt, "l");
		if (cmd.frame)
			plbox(gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
		if (cmd.xaxis)
			plaxes(gd.xorigin,gd.yorigin,gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,"",gd.frame_ytick,gd.frame_nysub);
		if (cmd.yaxis)
			plaxes(gd.xorigin,gd.yorigin,"",gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
	}

	plwid(defpenwidth);							// Return to default pen width 

    plcol0(cmd.labelcolor);							// Set labels color
	plschr(0.0, cmd.labelfontsize);					// Set the scaled size of the labels
	plwid(cmd.labelfontweight);									// Set font weight
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);				// text labels for bottom, left hand and top of the viewport
	plwid(defpenwidth);							// Return to default pen width 

	printf("\nsetting viewport ... done ...\n");

	for (ifile=0; ifile<gd.nfiles; ifile++) {
		if (gd.plotjoined[ifile]) {
			plcol0(cmd.linecolor+ifile);								// Line color
//			pllsty(cmd.linetype+ifile);							// Line type (8 styles). See 'plstyl'
			pllsty(gd.linetype[ifile]);							// Line type (8 styles). See 'plstyl'
			plwid(cmd.linewidth);							// Line width
			plline(vnpoint[ifile], xnc[ifile], ync[ifile]);		// Draw line through data points
			plwid(defpenwidth);							// Return to default pen width 
		}
		if (gd.withsymbols[ifile]) {
			plcol0(cmd.symbolcolor+ifile);								// symbol color
			plssym(0.0, cmd.symbolsize);
			plwid(cmd.symbolweight);							// Symbol weight
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
			plwid(defpenwidth);							// Return to default pen width 
		}
//		if (cmd.withdots) {
		if (gd.withdots[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plssym(0.0, 0.0);								// Set symbol to pixel
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
		}
	}

	if (!strnull(cmd.legends))
		setlegends(gd.nfiles);
	plottext();
	printf("\n\n... done log-normal plotting! ...\n\n");	
}

//-------------------------------- TERMINA PLOT2D LOGNORMAL ------------------------------------

//-------------------------------- COMIENZA PLOT2D LOGLOG ------------------------------------

local void plot2d_loglog(void)
{
//  real *xn, *yn;
//  pointptr p;
  int i, ifile;
  real *xnc[gd.nfiles], *ync[gd.nfiles];
  int vnpoint[gd.nfiles];
  real tmp;
  real viewportchr;

	printf("\n\nLog-Log plotting ...\n\n");
// COMIENZA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	for (ifile=0; ifile<gd.nfiles; ifile++) {

	if (strnull(cmd.inputfile))
		TestData(&vnpoint[ifile]);
	else
		InputData(gd.filenames[ifile],gd.vcol1[ifile], gd.vcol2[ifile], &vnpoint[ifile]);

//	InputData(filenames[ifile],vcol1[ifile], vcol2[ifile], &vnpoint[ifile]);

	printf("\nvnpoint[%d] = %d\n",ifile,vnpoint[ifile]);

	gd.npoint=vnpoint[ifile];
/*
	xn = (realptr) allocate(vnpoint[ifile]*sizeof(real));
	yn = (realptr) allocate(vnpoint[ifile]*sizeof(real));

	DO_POINT(p, pointtab, pointtab+vnpoint[ifile]) {
		i = p-pointtab;
		xn[i]=Xval(p);
		yn[i]=Yval(p);
	}
*/

	xnc[ifile] = &gd.xval[0]; ync[ifile]=&gd.yval[0];

//	DO_POINT(p, gd.pointtab, gd.pointtab+vnpoint[ifile]) {
	for (i=0; i<vnpoint[ifile]; ++i) {
//		i = p-gd.pointtab;
		tmp = rlog10(xnc[ifile][i]);
		xnc[ifile][i] = tmp;
		tmp = rlog10(ync[ifile][i]);
		ync[ifile][i] = tmp;
	}

	if (gd.x_autoscale) {
		if (ifile==0) {
			gd.xmin = xnc[ifile][0];
			gd.xmax = xnc[ifile][0]; 
			for (i = 1; i < vnpoint[ifile]; i++) {
				if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
				if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
			}
		} else
			for (i = 0; i < vnpoint[ifile]; i++) {
				if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
				if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
			}
	}

	if (gd.y_autoscale) {
		if (ifile==0) {
			gd.ymin = ync[ifile][0];
			gd.ymax = ync[ifile][0]; 
			for (i = 1; i < vnpoint[ifile]; i++) {
				if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
				if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
			}
		} else
			for (i = 0; i < vnpoint[ifile]; i++) {
				if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
				if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
			}
	}
	
	}

	printf("\nscaling done ...\n");

// TERMINA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	pllsty(cmd.axestype);								// Line type (8 styles). See 'plstyl'
	plwid(cmd.axeswidth);								// Line width
    plcol0(cmd.axescolor);								// Set color of frame - viewport

	pladv(0);										// Set to zero to advance to the next subpage

// Movidos aqui para que queden dentro del area de vision...
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, 1.1*viewportchr);		// Set the scaled size of the labels
	plvsta();										// Sets up a standard viewport, leaving apropriate margins

	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);							// Sets up the world-coordinates windows

    plsxax(cmd.nlxdigmax, 0);									// To set digmax on x-axis
    plsyax(cmd.nlydigmax, 0);									// To set digmax on y-axis
	plcol0(cmd.nlcolor);										// Set color for numeric labels. 
															// Also set color for axes...
	plschr(0.0,cmd.nlsize);

	if (strnull(cmd.framestyle)) {

	if (cmd.frame)
		if (!cmd.gridlines)
			plbox("bcnstl", 0.0, 0, "bcnstvl", 0.0, 0);		// Frame is drawn, axis, ticks, labels
		else
			plbox("bcgnstl", 0.0, 0, "bcgnstvl", 0.0, 0);	// Frame is drawn, axis, ticks, labels
	else
		plbox("bnstl", 0.0, 0, "bnstvl", 0.0, 0);			// Draw only x and y axis not a frame

	if (cmd.xaxis)
		plaxes(gd.xorigin, gd.yorigin, "anst", 0.0, 0, "", 0.0, 0);// Draw x-axis (y=yorigin)
															// In log scale can not be done when y=0
	if (cmd.yaxis)
		plaxes(gd.xorigin, gd.yorigin, "", 0.0, 0, "anstv", 0.0, 0);// Draw y-axis (x=xorigin). 
															// In log scale can not be done when x=0
	} else {
		strcat(gd.frame_xopt, "l");
		strcat(gd.frame_yopt, "l");
		if (cmd.frame)
			plbox(gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
		if (cmd.xaxis)
			plaxes(gd.xorigin,gd.yorigin,gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,"",gd.frame_ytick,gd.frame_nysub);
		if (cmd.yaxis)
			plaxes(gd.xorigin,gd.yorigin,"",gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
	}

	plwid(defpenwidth);										// Return to default pen width 

    plcol0(cmd.labelcolor);										// Set labels color
	plschr(0.0, cmd.labelfontsize);								// Set the scaled size of the labels
	plwid(cmd.labelfontweight);									// Set font weight
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);						// text labels for b, l hand and t of the viewport
	plschr(0.0, 1.0);								// Set the scaled size of the labels

	printf("\nsetting viewport ... done ...\n");

	for (ifile=0; ifile<gd.nfiles; ifile++) {
		if (gd.plotjoined[ifile]) {
			plcol0(cmd.linecolor+ifile);						// Line color
//			pllsty(cmd.linetype+ifile);							// Line type (8 styles). See 'plstyl'
			pllsty(gd.linetype[ifile]);							// Line type (8 styles). See 'plstyl'
			plwid(cmd.linewidth);								// Line width
			plline(vnpoint[ifile], xnc[ifile], ync[ifile]);	// Draw line through data points
			plwid(defpenwidth);								// Return to default pen width 
		}
		if (gd.withsymbols[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plwid(cmd.symbolweight);							// Symbol weight
			plssym(0.0, cmd.symbolsize);
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
			plwid(defpenwidth);								// Return to default pen width 
		}
//		if (cmd.withdots) {
		if (gd.withdots[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plssym(0.0, 0.0);								// Set symbol to pixel
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
		}
	}

	if (!strnull(cmd.legends))
		setlegends(gd.nfiles);
	plottext();
	printf("\n\n... done log-log plotting! ...\n\n");	
}

//-------------------------------- TERMINA PLOT2D LOGLOG ------------------------------------

//-------------------------------- COMIENZA PLOT2D GRAPHICSARRAY ------------------------------------

local void plot2d_graphicsarray(void)
{
//  real *xn, *yn;
//  pointptr p;
  int i, ifile;
  real *xnc[gd.nfiles], *ync[gd.nfiles];
  int vnpoint[gd.nfiles];
  real pxmin, pxmax, pymin, pymax;
  real viewportchr;

	printf("\n\nGraphics array plotting (%dx%d)...\n\n", gd.graphicsarray_nrows, gd.graphicsarray_ncols);
// COMIENZA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

for (ifile=0; ifile<gd.nfiles; ifile++) {

	InputData(gd.filenames[ifile],gd.vcol1[ifile], gd.vcol2[ifile], &vnpoint[ifile]);

	printf("\nvnpoint[%d] = %d\n",ifile,vnpoint[ifile]);

	gd.npoint=vnpoint[ifile];
/*
	xn = (realptr) allocate(vnpoint[ifile]*sizeof(real));
	yn = (realptr) allocate(vnpoint[ifile]*sizeof(real));

	DO_POINT(p, pointtab, pointtab+vnpoint[ifile]) {
		i = p-pointtab;
		xn[i]=Xval(p);
		yn[i]=Yval(p);
	}
*/

	xnc[ifile] = &gd.xval[0]; ync[ifile]=&gd.yval[0];

	if (gd.x_autoscale) {
		gd.xmin = xnc[ifile][0];
		gd.xmax = xnc[ifile][0]; 
		for (i = 1; i < vnpoint[ifile]; i++) {
			if (gd.xmin>xnc[ifile][i]) gd.xmin=xnc[ifile][i];
			if (gd.xmax<xnc[ifile][i]) gd.xmax=xnc[ifile][i];
		}
	}

	if (gd.y_autoscale) {
		gd.ymin = ync[ifile][0];
		gd.ymax = ync[ifile][0]; 
		for (i = 1; i < vnpoint[ifile]; i++) {
			if (gd.ymin>ync[ifile][i]) gd.ymin=ync[ifile][i];
			if (gd.ymax<ync[ifile][i]) gd.ymax=ync[ifile][i];
		}
	}

	printf("\nscaling done ...\n");

// TERMINA BLOQUE DE LECTURA DE ARCHIVOS Y GENERACION DEL ARREGLO DE X Y Y'S

	pllsty(cmd.axestype);								// Line type (8 styles). See 'plstyl'
	plwid(cmd.axeswidth);								// Line width
    plcol0(cmd.axescolor);								// Set color of frame - viewport

//    plenv(xmin, xmax, ymin, ymax, 0, 0);			// Set viewport. Axes scaled separately. Box drawn.

// The following commands are equivalent to the previous one

	pladv(0);										// Set to zero to advance to the next subpage

// Movidos aqui para que queden dentro del area de vision...
	viewportchr=MAX(cmd.labelfontsize,cmd.nlsize);
	plschr(0.0, 1.1*viewportchr);		// Set the scaled size of the labels
	plvsta();										// Sets up a standard viewport, leaving apropriate margins

	plwind(gd.xmin, gd.xmax, gd.ymin, gd.ymax);					// Sets up the world-coordinates windows

    plsxax(cmd.nlxdigmax, 0);						// To set digmax on x-axis
    plsyax(cmd.nlydigmax, 0);						// To set digmax on y-axis
	plcol0(cmd.nlcolor);							// Set color for numeric labels. Also set color for axes...
	plschr(0.0,cmd.nlsize);

	if (strnull(cmd.framestyle)) {

	if (cmd.frame)
		if (!cmd.gridlines)
			plbox("bcnst", 0.0, 0, "bcnstv", 0.0, 0);		// Frame is drawn, axis, ticks, labels
		else
			plbox("bcgnst", 0.0, 0, "bcgnstv", 0.0, 0);		// Frame is drawn, axis, ticks, labels
	else
		plbox("bnst", 0.0, 0, "bnstv", 0.0, 0);				// Draw only x and y axis not a frame

	if (cmd.xaxis)
		plaxes(gd.xorigin, gd.yorigin, "anst", 0.0, 0, "", 0.0, 0);// Draw x-axis (y=yorigin)
	if (cmd.yaxis)
		plaxes(gd.xorigin, gd.yorigin, "", 0.0, 0, "anstv", 0.0, 0);// Draw y-axis (x=xorigin)

	} else {
		if (cmd.frame)
			plbox(gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
		if (cmd.xaxis)
			plaxes(gd.xorigin,gd.yorigin,gd.frame_xopt,gd.frame_xtick,gd.frame_nxsub,"",gd.frame_ytick,gd.frame_nysub);
		if (cmd.yaxis)
			plaxes(gd.xorigin,gd.yorigin,"",gd.frame_xtick,gd.frame_nxsub,gd.frame_yopt,gd.frame_ytick,gd.frame_nysub);
	}

	plwid(defpenwidth);										// Return to default pen width 

    plcol0(cmd.labelcolor);										// Set labels color
	plschr(0.0, cmd.labelfontsize);								// Set the scaled size of the labels
	plwid(cmd.labelfontweight);									// Set font weight
    pllab(cmd.xlabel, cmd.ylabel, cmd.plotlabel);						// text labels for b, l hand and t of the viewport
	plwid(defpenwidth);										// Return to default pen width 

	printf("\nsetting viewport ... done ...\n");

		if (cmd.plotjoined) {
			plcol0(cmd.linecolor+ifile);						// Line color
//			pllsty(linetype+ifile);							// Line type (8 styles). See 'plstyl'
//			pllsty(cmd.linetype);								// Line type (8 styles). See 'plstyl'
			pllsty(gd.linetype[0]);								// Line type (8 styles). See 'plstyl'
			plwid(cmd.linewidth);								// Line width
			plline(vnpoint[ifile], xnc[ifile], ync[ifile]);	// Draw line through data points
			plwid(defpenwidth);								// Return to default pen width 
		}
		if (gd.withsymbols[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plssym(0.0, cmd.symbolsize);
			plwid(cmd.symbolweight);							// Symbol weight
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
			plwid(defpenwidth);								// Return to default pen width 
		}
		if (gd.withdots[ifile]) {
			plcol0(cmd.symbolcolor+ifile);						// symbol color
			plssym(0.0, 0.0);								// Set symbol to pixel
//			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], cmd.symboltype+ifile);
			plpoin(vnpoint[ifile], xnc[ifile], ync[ifile], gd.symboltype[ifile]);
		}

	if (!strnull(cmd.legends))
		setlegends(gd.nfiles);
	plottext();
	printf("\n\n... done array plotting %d ...\n\n", ifile);
}
}

//-------------------------------- TERMINA PLOT2D GRAPHICSARRAY ------------------------------------


// DIBUJO DE LEGENDAS PARA LOS SIMBOLOS

#define TOPLEFT		0
#define TOPRIGHT	1
#define BOTTONLEFT	2
#define BOTTONRIGHT	3

local void setlegends(int nlegends) {
  int i, j;
  real xu, yu, yfs, xstart, ystart, lsize, dyi, fl, fs, lwx, lwy;
  real xsym[1], ysym[1];

// Setting reference point
	xu = (gd.xmax-gd.xmin)/10.;
	yu = (gd.ymax-gd.ymin)/10.;
	yfs = 2.5*(gd.ymax-gd.ymin)/32.4;
	xstart = gd.xmin+1.0*xu;					// (x,y) por defecto. Puede ser dado en linea de comandos.
	ystart = gd.ymax-0.6*yu;
	lsize = xu/2.0;
	fl = 0.65;
	fs =1.0;
	lwx = 5.0*xu;
	lwy = 5.0*yu;

//	printf("\n\nAqui voy [1-1]\n");

	plschr(0.0, cmd.legendsfontsize);								// Set the scaled size of the labels
	plwid(cmd.legendsfontweight);									// Set font weight

    for (i = 0; i < nlegends; i++) {

	switch(cmd.legendspos) {

		case TOPLEFT:

//	printf("\n\nAqui voy [1-2]\n",gd.plotjoined[i]);

			dyi = (real)(i)*yfs;
//			if (cmd.plotjoined) {
			if (gd.plotjoined[i]) {
//	printf("\n\nAqui voy [1-2]\n");
				plcol0(cmd.linecolor+i);							// Line color
//				pllsty(cmd.linetype+i);								// Line type (8 styles). See 'plstyl'
				pllsty(gd.linetype[i]);								// Line type (8 styles). See 'plstyl'
				plwid(cmd.linewidth);								// Line width
				pljoin(xstart,		 ystart-dyi, xstart+lsize, ystart-dyi);
				plwid(defpenwidth);								// Return to default pen width 
			}
//			if (cmd.withsymbols) {
			if (gd.withsymbols[i]) {
				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, cmd.symbolsize);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-dyi;
				plwid(cmd.symbolweight);							// Symbol weight
//				plpoin(1, xsym, ysym, cmd.symboltype+i);
				plpoin(1, xsym, ysym, gd.symboltype[i]);
				plwid(defpenwidth);								// Return to default pen width 
			}
//			if (cmd.withdots) {
			if (gd.withdots[i]) {

//	printf("\n\nAqui voy [1-3]\n");

				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, 0.0);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-dyi;
				plpoin(1, xsym, ysym, 1);

//	printf("\n\nAqui voy [1-4]\n");

			}
			plptex(xstart+fs*xu, ystart-dyi, 0.1, 0.0, 0., gd.legendnames[i]);

//	printf("\n\nAqui voy [1-5]\n");

			break;

		case TOPRIGHT:
			dyi = (real)(i)*yfs;
//			if (cmd.plotjoined) {
//			if (gd.plotjoined[nlegends]) {
			if (gd.plotjoined[i]) {
				plcol0(cmd.linecolor+i);							// Line color
//				pllsty(cmd.linetype+i);								// Line type (8 styles). See 'plstyl'
				pllsty(gd.linetype[i]);								// Line type (8 styles). See 'plstyl'
				plwid(cmd.linewidth);								// Line width
				pljoin(xstart+lwx,		 ystart-dyi, xstart+lwx+lsize, ystart-dyi);
				plwid(defpenwidth);								// Return to default pen width 
			}
//			if (cmd.withsymbols) {
			if (gd.withsymbols[i]) {
				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, cmd.symbolsize);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-dyi;
				plwid(cmd.symbolweight);							// Symbol weight
//				plpoin(1, xsym, ysym, cmd.symboltype+i);
				plpoin(1, xsym, ysym, gd.symboltype[i]);
				plwid(defpenwidth);								// Return to default pen width
			}
//			if (cmd.withdots) {
			if (gd.withdots[i]) {
				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, 0.0);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-dyi;
				plpoin(1, xsym, ysym, 1);
			}
			plptex(xstart+lwx+fs*xu, ystart-dyi, 0.1, 0.0, 0., gd.legendnames[i]);
			break;


		case BOTTONLEFT:
			dyi = (real)(i)*yfs;
//			if (cmd.plotjoined) {
//			if (gd.plotjoined[nlegends]) {
			if (gd.plotjoined[i]) {
				plcol0(cmd.linecolor+i);							// Line color
//				pllsty(cmd.linetype+i);								// Line type (8 styles). See 'plstyl'
				pllsty(gd.linetype[i]);								// Line type (8 styles). See 'plstyl'
				plwid(cmd.linewidth);								// Line width
				pljoin(xstart,		 ystart-lwy-dyi, xstart+lsize, ystart-lwy-dyi);
				plwid(defpenwidth);								// Return to default pen width
			}
//			if (cmd.withsymbols) {
			if (gd.withsymbols[i]) {
				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, cmd.symbolsize);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plwid(cmd.symbolweight);							// Symbol weight
//				plpoin(1, xsym, ysym, cmd.symboltype+i);
				plpoin(1, xsym, ysym, gd.symboltype[i]);
				plwid(defpenwidth);								// Return to default pen width
			}
//			if (cmd.withdots) {
			if (gd.withdots[i]) {
				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, 0.0);
				xsym[0]=xstart+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plpoin(1, xsym, ysym, 1);
			}
			plptex(xstart+fs*xu, ystart-lwy-dyi, 0.1, 0.0, 0., gd.legendnames[i]);
			break;

		case BOTTONRIGHT:
			dyi = (real)(i)*yfs;
//			if (cmd.plotjoined) {
//			if (gd.plotjoined[nlegends]) {
			if (gd.plotjoined[i]) {
				plcol0(cmd.linecolor+i);							// Line color
//				pllsty(cmd.linetype+i);								// Line type (8 styles). See 'plstyl'
				pllsty(gd.linetype[i]);								// Line type (8 styles). See 'plstyl'
				plwid(cmd.linewidth);								// Line width
				pljoin(xstart+lwx,		 ystart-lwy-dyi, xstart+lwx+lsize, ystart-lwy-dyi);
				plwid(defpenwidth);								// Return to default pen width
			}
//			if (cmd.withsymbols) {
			if (gd.withsymbols[i]) {
				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, cmd.symbolsize);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plwid(cmd.symbolweight);							// Symbol weight
//				plpoin(1, xsym, ysym, cmd.symboltype+i);
				plpoin(1, xsym, ysym, gd.symboltype[i]);
				plwid(defpenwidth);								// Return to default pen width
			}
//			if (cmd.withdots) {
			if (gd.withdots[i]) {
				plcol0(cmd.symbolcolor+i);							// symbol color
				plssym(0.0, 0.0);
				xsym[0]=xstart+lwx+fl*xu;
				ysym[0]=ystart-lwy-dyi;
				plpoin(1, xsym, ysym, 1);
			}
			plptex(xstart+lwx+fs*xu, ystart-lwy-dyi, 0.1, 0.0, 0., gd.legendnames[i]);
			break;
		}
    }

	plwid(defpenwidth);										// Return to default pen width 
}

#undef TOPLEFT
#undef TOPRIGHT
#undef BOTTONLEFT
#undef BOTTONRIGHT

//#undef OFFSET

local void plottext(void) {
	if (!strnull(cmd.text1)) {
		plcol0(cmd.text1color);
		plschr(0.0, cmd.text1size);
		plwid(cmd.text1weight);								// Set fontweight
		plmtex(cmd.text1side,cmd.text1disp,cmd.text1pos,cmd.text1just,cmd.text1);
		plwid(defpenwidth);								// Return to default pen width
	}
	if (!strnull(cmd.text2)) {
		plcol0(cmd.text2color);
		plschr(0.0, cmd.text2size);
		plwid(cmd.text2weight);								// Set fontweight
		plmtex(cmd.text2side,cmd.text2disp,cmd.text2pos,cmd.text2just,cmd.text2);
		plwid(defpenwidth);								// Return to default pen width
	}
	if (!strnull(cmd.text3)) {
		plcol0(cmd.text3color);
		plschr(0.0, cmd.text3size);
		plwid(cmd.text3weight);								// Set fontweight
		plmtex(cmd.text3side,cmd.text3disp,cmd.text3pos,cmd.text3just,cmd.text3);
		plwid(defpenwidth);								// Return to default pen width
	}
	if (!strnull(cmd.text4)) {
		plcol0(cmd.text4color);
		plschr(0.0, cmd.text4size);
		plwid(cmd.text4weight);								// Set fontweight
		plptex(cmd.text4x,cmd.text4y,cmd.text4dx,cmd.text4dy,cmd.text4just,cmd.text4);
		plwid(defpenwidth);								// Return to default pen width
	}
}


