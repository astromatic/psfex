/*
 				cplot.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	CheckPlot functions
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for cplot.c
*
*	Last modify:	10/09/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _CPLOT_H_
#define _CPLOT_H_

/*------------------------------- constants ---------------------------------*/

#define		CPLOT_DEFRESX	800	/* Default X resol. for PNG and JPG */
#define		CPLOT_DEFRESY	600	/* Default X resol. for PNG and JPG */
#define		CPLOT_AAFAC	3	/* Anti-aliasing factor */
#define		CPLOT_NPOINTDEF	1024	/* Default number of points to plot */
#define		CPLOT_FGRIDLINES   9	/* Number of grid lines per axis */
#define		CPLOT_NDISTGRID	32	/* # of distort steps in each CCD dim*/
#define		CPLOT_ASTNSUBPLOTS 3	/* Number of subplot/dim/detector*/
#define		CPLOT_NTYPES	 128 /* Number of CPLOT types (typedef below)*/
#define		CPLOT_NSHADES	  32	/* Number of shading levels */

/*---------------------------- return messages ------------------------------*/
/*-------------------------------- macros -----------------------------------*/
/*--------------------------------- typedefs --------------------------------*/
typedef enum {CPLOT_NONE, CPLOT_FWHM, CPLOT_ELLIPTICITY, CPLOT_MOFFATRESI,
		CPLOT_ASYMRESI, CPLOT_COUNTS, CPLOT_COUNTFRAC}
		cplotenum;

typedef enum {CPLOT_NULL, CPLOT_XWIN, CPLOT_TK, CPLOT_XTERM, CPLOT_PLMETA,
	CPLOT_PS, CPLOT_PSC, CPLOT_XFIG, CPLOT_LJIIP, CPLOT_LJHPGL, CPLOT_IMP,
	CPLOT_PBM, CPLOT_PNG, CPLOT_JPEG, CPLOT_PSTEX, CPLOT_AQT, CPLOT_PDF,
	CPLOT_SVG} cplotdevenum;

typedef struct {cplotdevenum device; char *devname; char *extension;}
		devicestruct;

/*---------------------------------- svgp -----------------------------------*/
/*------------------------------- functions ---------------------------------*/

extern int		cplot_countfrac(fieldstruct *field),
			cplot_counts(fieldstruct *field),
			cplot_check(cplotenum cplottype),
			cplot_drawbounds(wcsstruct *wcsin, wcsstruct *wcsout),
			cplot_drawloccoordgrid(wcsstruct *wcs, double xmin,
					double xmax, double ymin, double ymax),
			cplot_ellipticity(fieldstruct *field),
			cplot_end(cplotenum cplottype),
			cplot_fwhm(fieldstruct *field),
			cplot_init(char *name, int nx, int ny,
				cplotenum cplottype),
			cplot_moffatresi(fieldstruct *field),
			cplot_asymresi(fieldstruct *field);
			
char			*cplot_degtosexal(char *str, double alpha,double step),
			*cplot_degtosexde(char *str, double delta,double step);
#endif

