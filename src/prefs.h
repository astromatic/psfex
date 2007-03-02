/*
 				types.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global type definitions.
*
*	Last modify:	01/03/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef	_CHECK_H_
#include "check.h"
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
  double	minsn;				/* Minimum S/N for patterns */
  double	maxelong;			/* Maximum A/B for patterns */
  double	maxvar;				/* Maximum FWHM variability */
  double	fwhmrange[2];			/* Allowed FWHM range */
  int		nfwhmrange;	       		/* nb of params */
  double	prof_accuracy;			/* Required PSF accuracy */
  double	psf_step;			/* Oversampling (pixels) */
  int		nsuper;				/* nb of supersampled pixels */
  int		autoselect_flag;		/* Auto. select FWHMs ? */
  int		recenter_flag;			/* Recenter PSF-candidates? */
  checkenum	check_type[MAXCHECK];		/* check-image types */
  int		ncheck_type;			/* nb of params */
  char		*(check_name[MAXCHECK]);	/* check-image names */
  int		ncheck_name;			/* nb of params */
  char		*(context_name[MAXCONTEXT]);	/* Names of context-keys */
  int		ncontext_name;			/* nb of params */
  int		context_group[MAXCONTEXT];	/* Context group */
  int		ncontext_group;			/* nb of params */
  int		context_nsnap;			/* nb of snapshots / context */
  int		group_deg[MAXCONTEXT];		/* Degree for each group */
  int		ngroup_deg;			/* nb of params */
  int		badpix_flag;			/* Filter bad pixels? */
  int		badpix_nmax;			/* Max number of bad pixels */
  int		pc_flag;			/* Include PCs? */
  char		pc_name[MAXCHAR];		/* PC filename */
  int		pc_npc;				/* Max. number of PCs */
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
