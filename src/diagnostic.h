/*
*				diagnostic.h
*
* Include file for diagnostic.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2006-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif
#ifndef _DIAGNOSTIC_H_
#define _DIAGNOSTIC_H_

/*----------------------------- Internal constants --------------------------*/
#define		PSF_DIAGMAXITER	1000	/* Max. nb of iterations in fitting */
#define		PSF_DIAGNPARAM	7	/* Number of fitted parameters */
#define		PSF_FWHMMIN	0.1	/* Minimum FWHM for fit (model pixels)*/
#define		PSF_BETAMIN	0.5	/* Minimum Moffat beta for fit */
#define		PSF_NSUBPIX	5	/* Oversamp. factor to mimick top-hat */

/*-------------------------------- macros -----------------------------------*/

#define         PSFEX_POW(x,a)	(x>0.01? exp(a*log(x)) : pow(x,a))

/*----------------------------- Global variables ---------------------------*/

float		moffat_parammin[PSF_DIAGNPARAM],moffat_parammax[PSF_DIAGNPARAM];

/*---------------------------------- protos --------------------------------*/
extern void	psf_boundtounbound(float *param, double *dparam),
		psf_compdiag(psfstruct *psf, moffatstruct *moffat,
			double *dpos, int oversamp),
		psf_diagnostic(psfstruct *psf),
		psf_diagprintout(int n_par, float *par, int m_dat,
			float *fvec, void *data, int iflag,int iter,int nfev),
		psf_diagresi(double *par, double *fvec, int m, int n,
			void *adata),
		psf_moffat(psfstruct *psf, moffatstruct *moffat),
		psf_unboundtobound(double *dparam, float *param),
		psf_wcsdiagnostic(psfstruct *psf, wcsstruct *wcs);

extern double	psf_noiseqarea(psfstruct *psf),
		psf_normresi(float *par, psfstruct *psf),
		psf_symresi(psfstruct *psf);

#endif

