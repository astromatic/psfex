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
*	Last modify:	25/02/2002
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
#define		PSF_NSNAP	7	/* Number of PSF snapshots/dimension */
/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
extern void	psf_diagnostic(psfstruct *psf, out_data_struct *out),
		psf_diagprintout(int n_par, double *par, int m_dat,
			double *fvec, void *data, int iflag,int iter,int nfev),
		psf_diagresi(double *par, int m_dat, double *fvec, void *data,
			int *info);

#endif

