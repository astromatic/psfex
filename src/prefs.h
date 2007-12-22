/*
 				prefs.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global type definitions.
*
*	Last modify:	22/12/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef	_CHECK_H_
#include "check.h"
#endif

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _PREFS_H_
#define _PREFS_H_

/*----------------------------- Internal constants --------------------------*/

#define	MAXLIST		32		/* max. nb of list members */

/* NOTES:
One must have:	MAXLIST >= 1 (preferably >= 16!)
*/
/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		**command_line;			/* Command line */
  int		ncommand_line;			/* nb of params */
  char		prefs_name[MAXCHAR];		/* prefs filename */
  char		psf_name[MAXCHAR];		/* PSF filename */
  int		retisize[2], nretisize;		/* Retina size */
  enum {NEWBASIS_NONE, NEWBASIS_PCAMULTI, NEWBASIS_PCASINGLE}
		newbasis_type;			/* Type of new basis */
  int		newbasis_number;		/* Number of PCs */
/* Point-source sample */
  double	minsn;				/* Minimum S/N for patterns */
  double	maxelong;			/* Maximum A/B for patterns */
  double	maxvar;				/* Maximum FWHM variability */
  double	fwhmrange[2];			/* Allowed FWHM range */
  int		nfwhmrange;	       		/* nb of params */
  int		flag_mask;			/* Rejection mask on SEx FLAGS*/
  double	prof_accuracy;			/* Required PSF accuracy */
  double	psf_step;			/* Oversampling (pixels) */
  int		badpix_flag;			/* Filter bad pixels? */
  int		badpix_nmax;			/* Max number of bad pixels */
/* Vector basis */
  basistypenum	basis_type;			/* PSF vector basis set */
  int		basis_number;			/* nb of supersampled pixels */
  char		basis_name[MAXCHAR];		/* PSF vector basis filename */
  double	basis_scale;			/* Gauss-Laguerre beta param */
/* Re-centering */
  int		autoselect_flag;		/* Auto. select FWHMs ? */
  int		recenter_flag;			/* Recenter PSF-candidates? */
/* Check-images */
  checkenum	check_type[MAXCHECK];		/* check-image types */
  int		ncheck_type;			/* nb of params */
  char		*(check_name[MAXCHECK]);	/* check-image names */
  int		ncheck_name;			/* nb of params */
  int		check_cubeflag;			/* check-images as datacubes?*/
/* PSF variability */
  char		*(context_name[MAXCONTEXT]);	/* Names of context-keys */
  int		ncontext_name;			/* nb of params */
  int		context_group[MAXCONTEXT];	/* Context group */
  int		ncontext_group;			/* nb of params */
  int		context_nsnap;			/* nb of snapshots / context */
  int		group_deg[MAXCONTEXT];		/* Degree for each group */
  int		ngroup_deg;			/* nb of params */
/* Homogenisation kernel vector basis */
  enum	{HOMOBASIS_NONE, HOMOBASIS_GAUSSLAGUERRE}
		homobasis_type;			/* Homo. kernel basis set */
  int		homobasis_number;		/* nb of supersampled pixels */
  char		homokernel_name[MAXCHAR];	/* Homo. kernel filename */
  double	homobasis_scale;		/* Gauss-Laguerre beta param */
  double	homopsf_params[2];		/* Idealised Moffat PSF params*/
  int		nhomopsf_params;		/* nb of params */
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

  prefstruct		prefs;

/*-------------------------------- protos -----------------------------------*/
extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);
#endif
