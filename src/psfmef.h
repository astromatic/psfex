 /*
 				psfmef.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for psfmef.c.
*
*	Last modify:	13/03/2008
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

typedef struct psfmef
  {
  char		catname[MAXCHAR];	/* Input catalog filename */
  char		*rcatname;		/* "Reduced" catalog name */
  int		next;			/* Number of extensions */
  psfstruct	**psf;			/* Array of PSFs */
  wcsstruct	**wcs;			/* Array of WCS structures */
  }	psfmefstruct;

/*---------------------------------- protos --------------------------------*/
extern void	psfmef_end(psfmefstruct *psfmef),
		psfmef_save(psfmefstruct *psfmef, char *filename);

extern psfmefstruct	*psfmef_init(char *catname);


void		context_apply(contextstruct *context, psfstruct *psf,
			psfmefstruct **psfmefs, int ext, int ncat),
		context_end(contextstruct *context);

#endif

