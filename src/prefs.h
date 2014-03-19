/*
*				prefs.h
*
* Include file for prefs.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	PSFEx is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	PSFEx is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with PSFEx.  If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		26/02/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _CATOUT_H_
#include "catout.h"
#endif

#ifndef	_CHECK_H_
#include "check.h"
#endif

#ifndef _CONTEXT_H_
#include "context.h"
#endif

#ifndef _CPLOT_H_
#include "cplot.h"
#endif

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _PREFS_H_
#define _PREFS_H_

/*----------------------------- Internal constants --------------------------*/

#define	MAXCHARL	16384		/* max. nb of chars in a string list */
#define	MAXLIST		64		/* max. nb of list members */
#define	MAXLISTSIZE	2000000		/* max size of list */

/* NOTES:
One must have:	MAXLIST >= 1 (preferably >= 16!)
*/
/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		**command_line;			/* Command line */
  int		ncommand_line;			/* nb of params */
  char		prefs_name[MAXCHAR];		/* prefs filename */
  char		*(incat_name[MAXFILE]);		/* Filename(s) of input cats */
  int		ncat;				/* Number of input images */
  char		psf_dir[MAXCHAR];		/* PSF output dir */
  char		psf_suffix[MAXCHAR];		/* Suffix for PSF filenames */
  int		psf_size[2], npsf_size;		/* PSF size */
  enum {NEWBASIS_NONE, NEWBASIS_PCAINDEPENDENT, NEWBASIS_PCACOMMON}
		newbasis_type;			/* Type of new basis */
  int		newbasis_number;		/* Number of PCs */
/* Point-source sample */
  double	minsn;				/* Minimum S/N for patterns */
  double	maxellip;			/* Maximum (A-B)/(A+B) */
  double	maxvar;				/* Maximum FWHM variability */
  double	fwhmrange[2];			/* Allowed FWHM range */
  int		nfwhmrange;	       		/* nb of params */
  int		flag_mask;			/* Rejection mask on SEx FLAGS*/
  int		wflag_mask;			/* Rej. mask on FLAGS_WEIGHT */
  int		imaflag_mask;			/* Rej. mask on IMAFLAGS_ISO */
  int		nmax;				/* Max. nb of samples per set*/
  double	prof_accuracy;			/* Required PSF accuracy */
  double	psf_step;			/* Oversampling (pixels) */
  double	psf_pixsize[2];			/* Eff. pixel size (pixels) */
  int		npsf_pixsize;			/* nb of params */
  int		badpix_flag;			/* Filter bad pixels? */
  int		badpix_nmax;			/* Max number of bad pixels */
  char		photflux_key[MAXCHAR];		/* Name of phot. flux key */
  char		photflux_rkey[MAXCHAR];		/* Reduced phot. flux key */
  char		photfluxerr_key[MAXCHAR];	/* Name of phot. flux err. key*/
  char		photfluxerr_rkey[MAXCHAR];	/* Reduced phot. flux err. key*/
/* Vector basis */
  basistypenum	basis_type;			/* PSF vector basis set */
  int		basis_number;			/* nb of supersampled pixels */
  char		basis_name[MAXCHAR];		/* PSF vector basis filename */
  double	basis_scale;			/* Gauss-Laguerre beta param */
/* Re-centering */
  char		*(center_key[2]);		/* Names of centering keys */
  int		ncenter_key;			/* nb of params */
  int		autoselect_flag;		/* Auto. select FWHMs ? */
  int		recenter_flag;			/* Recenter PSF-candidates? */
/* Check-images */
  checkenum	check_type[MAXCHECK];		/* check-image types */
  int		ncheck_type;			/* nb of params */
  char		*(check_name[MAXCHECK]);	/* check-image names */
  int		ncheck_name;			/* nb of params */
  int		check_cubeflag;			/* check-images as datacubes?*/
/* PSF variability */
  enum {VAR_NONE, VAR_SEEING}	var_type;	/* PSF variability type */
  char		*(context_name[MAXCONTEXT]);	/* Names of context-keys */
  int		ncontext_name;			/* nb of params */
  int		context_group[MAXCONTEXT];	/* Context group */
  int		ncontext_group;			/* nb of params */
  int		context_nsnap;			/* nb of snapshots / context */
  int		group_deg[MAXCONTEXT];		/* Degree for each group */
  int		ngroup_deg;			/* nb of params */
  enum	{HIDDEN_MEF_INDEPENDENT, HIDDEN_MEF_COMMON}
		hidden_mef_type;		/* Mosaic handling for hiddens*/
  enum	{STABILITY_EXPOSURE, STABILITY_SEQUENCE}
		stability_type;			/* PSF stability range*/
  enum	{PSF_MEF_INDEPENDENT, PSF_MEF_COMMON}
		psf_mef_type;			/* Mosaic handling for PSF */
/* Homogenisation kernel vector basis */
  enum	{HOMOBASIS_NONE, HOMOBASIS_GAUSSLAGUERRE}
		homobasis_type;			/* Homo. kernel basis set */
  int		homobasis_number;		/* nb of supersampled pixels */
  char		homokernel_dir[MAXCHAR];	/* Homo. kernel output dir */
  char		homokernel_suffix[MAXCHAR];	/* Homo. kernel file suffix */
  double	homobasis_scale;		/* Gauss-Laguerre beta param */
  double	homopsf_params[2];		/* Idealised Moffat PSF params*/
  int		nhomopsf_params;		/* nb of params */
/* Output catalogs */
  char		outcat_name[MAXCHAR];		/* Output filename */
  cattypenum	outcat_type;			/* Output catalog type */
  int		outcatpipe_flag;		/* Pipe output catalogs? */
/* Check-plots */
  cplotenum	cplot_device[MAXCHECK];		/* check-plot format */
  int		ncplot_device;			/* nb of params */
  cplotenum	cplot_type[MAXCHECK];		/* check-plot types */
  int		ncplot_type;			/* nb of params */
  char		*(cplot_name[MAXCHECK]);	/* check-plot names */
  int		ncplot_name;			/* nb of params */
  int		cplot_flag;			/* = 0 if no check-plot */
  int		cplot_res[2];			/* X,Y check-plot resolution */
  int		ncplot_res;			/* nb of params */
  int		cplot_antialiasflag;		/* Anti-aliasing on/off */
/* Multithreading */
  int		nthreads;			/* Number of active threads */
/* Misc */
  enum {QUIET, NORM, LOG, FULL}	verbose_type;	/* How much it displays info */
  int		xml_flag;			/* Write XML file? */
  char		xml_name[MAXCHAR];		/* XML file name */
  char		xsl_name[MAXCHAR];		/* XSL file name (or URL) */
  char		sdate_start[12];		/* PSFEx start date */
  char		stime_start[12];		/* PSFEx start time */
  char		sdate_end[12];			/* PSFEx end date */
  char		stime_end[12];			/* PSFEx end time */
  double	time_diff;			/* Execution time */
  }	prefstruct;

extern prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern char	*list_to_str(char *listname);

extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);
#endif
