/*
*				psf.h
*
* Include file for psf.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		20/02/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _POLY_H_
#include "wcs/poly.h"
#endif

#ifndef _SAMPLE_H_
#include "sample.h"
#endif

#ifndef _PSF_H_
#define _PSF_H_

/*----------------------------- Internal constants --------------------------*/

#define	PSF_NODIAG	0	/* Don't do diagnostics */
#define	PSF_DIAG	1	/* Do diagnostics */
#define	PSF_FREEDFACTOR	1.1	/* Margin against overfitting (10%) */
#define	PSF_NMASKDIM	3	/* Number of dimensions for PSF data */
#define	PSF_MAXSHIFT	3.0	/* Max shift from initial guess (pixels)*/
#define	PSF_MINSHIFT	1e-4	/* Min shift from previous guess (pixels)*/
#define PSF_NITER	40	/* Maximum number of iterations in fit */
#define	PSF_NSNAPMAX	16	/* Maximum number of PSF snapshots/dimension */
#define	GAUSS_LAG_OSAMP	3	/* Gauss-Laguerre oversampling factor */
#define	PSF_AUTO_FWHM	3.0	/* FWHM theshold for PIXEL-AUTO mode */
#define	PSF_NORTHOSTEP	16	/* Number of PSF orthonor. snapshots/dimension*/

/*----------------------------- Type definitions --------------------------*/
typedef enum {BASIS_NONE, BASIS_PIXEL, BASIS_GAUSS_LAGUERRE, BASIS_FILE,
		BASIS_PIXEL_AUTO}
        basistypenum;
/*--------------------------- structure definitions -------------------------*/

typedef struct moffat
  {
  double	context[POLY_MAXDIM];	/* Context coordinates */
  float		amplitude;	/* Central amplitude */
  float		xc[2];		/* Center coordinates */
  float		fwhm_min;	/* FWHM along the minor axis */
  float		fwhm_max;	/* FWHM along the major axis */
  float		theta;		/* Position angle of the major axis / NAXIS1 */
  float		beta;		/* Moffat beta parameter */
  float		residuals;	/* Normalized residuals */
  float		symresiduals;	/* Normalized symmetry residuals */
  float		noiseqarea;	/* Noise equivalent area (pixels^2) */
  int		nsubpix;	/* Number of supersampled pixels */
  }	moffatstruct;

typedef struct psf
  {
  int		dim;		/* Dimensionality of the tabulated data */
  int		*size;		/* PSF dimensions */
  int		npix;		/* Total number of involved PSF pixels */
  float		*comp; 		/* Complete pix. data (PSF components) */
  float		*loc;		/* Local PSF */
  float		*resi;		/* Map of residuals */
  char		**contextname;	/* Array of context key-names */
  double	*contextoffset;	/* Offset to apply to context data */
  double	*contextscale;	/* Scaling to apply to context data */
  int		cx,cy;		/* Indices of X and Y mapping contexts */
  struct poly	*poly;		/* Polynom describing the PSF variations */
  float		pixstep;	/* Mask oversampling (pixel). */
  float		pixsize[2];	/* Effective pixel size on each axis (pixel) */
  int		samples_loaded;	/* Number of detections loaded */
  int		samples_accepted;/* Number of detections accepted */
  double	chi2;		/* chi2/d.o.f. */
  float		fwhm;		/* Initial guess of the FWHM */
  int		*pixmask;	/* Pixel mask for local bases */
  float		*basis;		/* Basis vectors */
  float		*basiscoeff;	/* Basis vector coefficients */
  int		nbasis;		/* Number of basis vectors */
  int		ndata;		/* Size of the design matrix along data axis */
  int		nsnap;		/* Total number of snapshots */
  int		nmed;		/* Median position amongst snapshots */
  int		nsubpix;	/* Number of intrapixel samples per axis */
  moffatstruct	*moffat;	/* Array of Moffat fits to PSF */
  moffatstruct	*pfmoffat;	/* Array of pixel-free Moffat fits to PSF */
  float		moffat_fwhm_min;
  float		moffat_fwhm;	/* Central Moffat FWHM */
  float		moffat_fwhm_max;
  float		moffat_fwhm_wcs_min;
  float		moffat_fwhm_wcs;	/* Average Moffat FWHM in arcsec*/
  float		moffat_fwhm_wcs_max;
  float		moffat_ellipticity_min;
  float		moffat_ellipticity;	/* Central Moffat ellipticity */
  float		moffat_ellipticity_max;
  float		moffat_ellipticity1_min;
  float		moffat_ellipticity1;	/* Central Moffat e1 */
  float		moffat_ellipticity1_max;
  float		moffat_ellipticity2_min;
  float		moffat_ellipticity2;	/* Central Moffat e2 */
  float		moffat_ellipticity2_max;
  float		moffat_beta_min;
  float		moffat_beta;	/* Central Moffat beta */
  float		moffat_beta_max;
  float		moffat_residuals_min;
  float		moffat_residuals;/* Central Moffat residuals */
  float		moffat_residuals_max;
  float		moffat_score_min;
  float		moffat_score;	/* Central pixel-free Moffat score */
  float		moffat_score_max;
  float		pfmoffat_fwhm_min;
  float		pfmoffat_fwhm;	/* Central pixel-free Moffat FWHM */
  float		pfmoffat_fwhm_max;
  float		pfmoffat_fwhm_wcs_min;
  float		pfmoffat_fwhm_wcs; /* Average pixel-free Moffat FWHM in arcsec*/
  float		pfmoffat_fwhm_wcs_max;
  float		pfmoffat_ellipticity_min;
  float		pfmoffat_ellipticity;	/* Central pix-free Moffat ellipticity*/
  float		pfmoffat_ellipticity_max;
  float		pfmoffat_ellipticity1_min;
  float		pfmoffat_ellipticity1;	/* Central pix-free Moffat e1 */
  float		pfmoffat_ellipticity1_max;
  float		pfmoffat_ellipticity2_min;
  float		pfmoffat_ellipticity2;	/* Central pix-free Moffat e2 */
  float		pfmoffat_ellipticity2_max;
  float		pfmoffat_beta_min;
  float		pfmoffat_beta;	/* Central pixel-free Moffat beta */
  float		pfmoffat_beta_max;
  float		pfmoffat_residuals_min;
  float		pfmoffat_residuals;/* Central pixel-free Moffat residuals */
  float		pfmoffat_residuals_max;
  float		sym_residuals_min;
  float		sym_residuals;/* Symmetry residuals */
  float		sym_residuals_max;
  float		noiseqarea_min;
  float		noiseqarea;	/* Noise equivalent area */
  float		noiseqarea_max;
  float		pixscale_wcs_min;
  float		pixscale_wcs;	/* Average pixel scale in arcsec */
  float		pixscale_wcs_max;
  float		*homo_kernel;		/* PSF homogenization kernel */
  double	homopsf_params[2];	/* Idealised Moffat PSF params*/
  int		homobasis_number;	/* nb of supersampled pixels */
  }	psfstruct;


/*---------------------------------- protos --------------------------------*/
extern void	psf_build(psfstruct *psf, double *pos),
		psf_clip(psfstruct *psf),
		psf_end(psfstruct *psf),
		psf_make(psfstruct *psf, setstruct *set, double prof_accuracy),
		psf_makebasis(psfstruct *psf, setstruct *set,
			basistypenum basis_type,  int nvec),
		psf_makeresi(psfstruct *psf, setstruct *set, int centflag,
			double prof_accuracy),
		psf_makemask(psfstruct *psf, setstruct *set, double chithresh),
		psf_orthopoly(psfstruct *psf, setstruct *set),
		psf_save(psfstruct *psf,  char *filename, int ext, int next);

extern int	psf_pshapelet(float **shape, int w, int h, int nmax,
			double beta),
		psf_readbasis(psfstruct *psf, char *filename, int ext),
		psf_refine(psfstruct *psf, setstruct *set);

extern double	psf_chi2(psfstruct *psf, setstruct *set),
		psf_clean(psfstruct *psf, setstruct *set, double prof_accuracy);

extern psfstruct	*psf_copy(psfstruct *psf),
			*psf_inherit(contextstruct *context, psfstruct *psf),
			*psf_init(contextstruct *context, int *size,
				float psfstep, float *pixsize, int nsample),
			*psf_load(char *filename);

#endif

