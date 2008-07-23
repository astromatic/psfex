 /*
 				pca.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for pca.c.
*
*	Last modify:	23/07/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PCA_H_
#define	_PCA_H_

#ifndef _PSF_H_
#include "psf.h"
#endif

/*----------------------------- Internal constants --------------------------*/

#define		PCA_NSNAP	5	/* Number of points per PSFVar dim. */
#define		PCA_NITER	100	/* Max nb of iter. in pc_find() */
#define		PCA_CONVEPS	1e-6	/* pc_find() converg. criterion */

/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
extern double	*pca_oncomps(psfstruct **psfs, int ncat, int npc),
		pca_findpc(double *covmat, float *vec, int nmat);

extern float	*pca_onsnaps(psfstruct **psfs, int ncat, int npc);

#endif
