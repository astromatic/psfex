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
*	Last modify:	17/07/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _DIAGNOSTIC_H_
#define _DIAGNOSTIC_H_

/*----------------------------- Internal constants --------------------------*/
#define		PSF_DIAGMAXITER	200	/* Max. nb of iterations in fitting */
#define		PSF_DIAGNPARAM	7	/* Number of fitted parameters */
#define		PSF_FWHMMIN	0.1	/* Minimum FWHM for fit (model pixels)*/
#define		PSF_FWHMMAX	100.0	/* Maximum FWHM for fit (model pixels)*/
#define		PSF_BETAMIN	0.1	/* Minimum Moffat beta for fit */
/*----------------------------- Global variables ---------------------------*/

double	moffat_parammin[PSF_DIAGNPARAM], moffat_parammax[PSF_DIAGNPARAM];

/*---------------------------------- protos --------------------------------*/
extern void	psf_boundtounbound(double *param),
		psf_diagnostic(psfstruct *psf),
		psf_diagprintout(int n_par, double *par, int m_dat,
			double *fvec, void *data, int iflag,int iter,int nfev),
		psf_diagresi(double *par, double *fvec, int m, int n,
			void *adata),
		psf_moffat(psfstruct *psf, moffatstruct *moffat),
		psf_unboundtobound(double *param);

extern double	psf_normresi(double *par, psfstruct *psf),
		psf_symresi(psfstruct *psf);

#endif

