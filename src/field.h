/*
*				field.h
*
* Include file for field.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		25/06/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _CONTEXT_H_
#include "context.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _SAMPLE_H_
#include "sample.h"
#endif

#ifndef _PSFMEF_H_
#define _PSFMEF_H_

/*----------------------------- Internal constants --------------------------*/
#define	COUNT_LOADED	1		/* Count detections that are loaded */
#define	COUNT_ACCEPTED	2		/* Count detections that are accepted */

/*------------------------------ Type definitions ---------------------------*/
/*--------------------------- structure definitions -------------------------*/

typedef struct field
  {
  char		catname[MAXCHAR];	/* Input catalog filename */
  char		*rcatname;		/* "Reduced" catalog name */
  char		rtcatname[MAXCHAR];	/* "Reduced", no trail catalog name */
  char		ident[MAXCHAR];		/* Field identifier (read from FITS) */
  int		next;			/* Number of extensions */
  int		ndet;			/* Number of detections (info only) */
  psfstruct	**psf;			/* Array of PSFs */
  wcsstruct	**wcs;			/* Array of WCS structures */
  setstruct	*set;			/* Array of catalogues */
  catstruct	**ccat;			/* Pointers to check-image files */
  double	meanwcspos[NAXIS];	/* Mean pixel coordinate */
  double	meanwcsscale[NAXIS];	/* Mean pixel scale */
  double	maxradius;		/* Maxium radius */
  int		**lcount;		/* Count detections that are loaded */
  int		**acount;		/* Count detections that are accepted */
  int		**count;		/* Count detections in stats */
  double	**modchi2;		/* Sum of chi2's per image area */
  double	**modresi;		/* Sum of res. indices per image area */
  }	fieldstruct;

/*---------------------------------- protos --------------------------------*/
extern fieldstruct	*field_init(char *catname);

extern void		field_count(fieldstruct **fields, setstruct *set,
				int counttype),
			field_end(fieldstruct *field),
			field_locate(fieldstruct *field),
			field_psfsave(fieldstruct *field, char *filename),
			field_stats(fieldstruct **fields, setstruct *set);

void			context_apply(contextstruct *context, psfstruct *psf,
				fieldstruct **fields, int ext, int catindex,
				int ncat),
			context_end(contextstruct *context);

#endif

