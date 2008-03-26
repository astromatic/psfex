 /*
 				field.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for field.c.
*
*	Last modify:	21/03/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _CONTEXT_H_
#include "context.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _PSFMEF_H_
#define _PSFMEF_H_

/*----------------------------- Internal constants --------------------------*/

/*------------------------------ Type definitions ---------------------------*/
/*--------------------------- structure definitions -------------------------*/

typedef struct field
  {
  char		catname[MAXCHAR];	/* Input catalog filename */
  char		*rcatname;		/* "Reduced" catalog name */
  int		next;			/* Number of extensions */
  psfstruct	**psf;			/* Array of PSFs */
  wcsstruct	**wcs;			/* Array of WCS structures */
  setstruct	*set;			/* Array of catalogues */
  double	meanwcspos[NAXIS];	/* Mean pixel coordinate */
  double	meanwcsscale[NAXIS];	/* Mean pixel scale */
  double	maxradius;		/* Maxium radius */
  }	fieldstruct;

/*---------------------------------- protos --------------------------------*/
extern fieldstruct	*field_init(char *catname);

extern double		dhmedian(double *ra, int n);

extern void		field_end(fieldstruct *field),
			field_locate(fieldstruct *field),
			field_psfsave(fieldstruct *field, char *filename);

void			context_apply(contextstruct *context, psfstruct *psf,
				fieldstruct **fields, int ext, int ncat),
			context_end(contextstruct *context);

#endif

