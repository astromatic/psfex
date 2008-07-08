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
*	Last modify:	08/07/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef	_FIELD_H_
#include "field.h"
#endif

#ifndef _CHECK_H_
#define _CHECK_H_

/*----------------------------- Internal constants --------------------------*/

#define		MAXCHECK	16		/* max. # of CHECKimages */

/*----------------------------- Type definitions --------------------------*/
typedef enum {PSF_NONE, PSF_BASIS, PSF_CHI, PSF_PROTO, PSF_RESIDUALS,
		PSF_SAMPLES, PSF_SNAPSHOTS, PSF_SNAPSHOTS_IMRES,
		PSF_WEIGHTS, PSF_MOFFAT,PSF_SUBMOFFAT,PSF_SUBSYM}
	checkenum;

/*---------------------------------- protos --------------------------------*/
extern void		check_write(fieldstruct *field,	char *checkname,
				checkenum checktype, int ext, int next,
				int cubeflag);

#endif

