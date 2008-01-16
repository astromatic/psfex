 /*
 				homo.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for homo.c.
*
*	Last modify:	15/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _HOMO_H_
#define _HOMO_H_

/*----------------------------- Internal constants --------------------------*/

#define		HOMO_NSNAP	5	/* Number of points per PSFVar dim. */

/*----------------------------- Global variables ---------------------------*/
/*---------------------------------- protos --------------------------------*/
extern void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
			int homobasis_number, double homobasis_scale,
			int ext, int next),
		psf_savehomo(psfstruct *psf, char *filename, int ext, int next);
#endif

