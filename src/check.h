 /*
 				check.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for producing check-images.
*
*	Last modify:	12/11/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef	_SAMPLE_H_
#include "sample.h"
#endif

#ifndef	_PSF_H_
#include "psf.h"
#endif

#ifndef _CHECK_H_
#define _CHECK_H_

/*----------------------------- Internal constants --------------------------*/

#define		MAXCHECK	16		/* max. # of CHECKimages */

/*----------------------------- Type definitions --------------------------*/
typedef enum {PSF_NONE, PSF_BASIS, PSF_CHI, PSF_PROTO, PSF_RESIDUALS,
		PSF_RAWDATA, PSF_SAMPLES, PSF_SNAPSHOTS, PSF_SNAPSHOTS_IMRES,
		PSF_WEIGHTS, PSF_MOFFAT,PSF_SUBMOFFAT,PSF_SUBSYM}
	checkenum;

/*---------------------------------- protos --------------------------------*/
extern void		psf_writecheck(psfstruct *psf, setstruct *set,
					char *filename, checkenum checktype,
				int ext, int next, int cubeflag);

#endif

