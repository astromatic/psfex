 /*
 				diagnostic.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for diagnostic.c.
*
*	Last modify:	26/02/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _DIAGNOSTIC_H_
#define _DIAGNOSTIC_H_

/*----------------------------- Internal constants --------------------------*/
#define		PSF_DIAGMAXITER	100	/* Max. nb of iterations in fitting */
#define		PSF_DIAGNPARAM	7	/* Number of fitted parameters */
/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
extern void	psf_diagnostic(psfstruct *psf),
		psf_diagprintout(int n_par, double *par, int m_dat,
			double *fvec, void *data, int iflag,int iter,int nfev),
		psf_diagresi(double *par, int m_dat, double *fvec, void *data,
			int *info),
		psf_moffat(psfstruct *psf, moffatstruct *moffat);

extern double	psf_normresi(double *par, psfstruct *psf);

#endif

