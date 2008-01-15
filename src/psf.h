 /*
 				psf.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for psf.c.
*
*	Last modify:	15/01/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _POLY_H_
#include "poly.h"
#endif

#ifndef _SAMPLE_H_
#include "sample.h"
#endif

#ifndef _PSF_H_
#define _PSF_H_

/*----------------------------- Internal constants --------------------------*/

#define	PSF_FREEDFACTOR	1.0	/* Margin against overfitting */
#define	PSF_NMASKDIM	3	/* Number of dimensions for PSF data */
#define	PSF_MAXSHIFT	3.0	/* Max shift from initial guess (pixels)*/
#define	PSF_MINSHIFT	1e-4	/* Min shift from previous guess (pixels)*/
#define PSF_NITER	40	/* Maximum number of iterations in fit */
#define	PSF_NSNAPMAX	16	/* Maximum number of PSF snapshots/dimension */
#define	GAUSS_LAG_OSAMP	3	/* Gauss-Laguerre oversampling factor */
#define	PSF_AUTO_FWHM	3.0	/* FWHM theshold for PIXEL-AUTO mode */

/*----------------------------- Type definitions --------------------------*/
typedef enum {BASIS_NONE, BASIS_PIXEL, BASIS_GAUSS_LAGUERRE, BASIS_FILE,
		BASIS_PIXEL_AUTO}
        basistypenum;
/*--------------------------- structure definitions -------------------------*/

typedef struct moffat
  {
  double	context[POLY_MAXDIM];	/* Context coordinates */
  double	amplitude;	/* Central amplitude */
  double	xc[2];		/* Center coordinates */
  double	fwhm_min;	/* FWHM along the minor axis */
  double	fwhm_max;	/* FWHM along the major axis */
  double	theta;		/* Position angle of the major axis / NAXIS1 */
  double	beta;		/* Moffat beta parameter */
  double	residuals;	/* Normalized residuals */
  double	symresiduals;	/* Normalized symmetry residuals */
  }	moffatstruct;

typedef struct psf
  {
  int		dim;		/* Dimensionality of the tabulated data */
  int		*size;		/* PSF  dimensions */
  int		npix;		/* Total number of involved PSF pixels */
  float		*comp; 		/* Complete pix. data (principal components) */
  float		*loc;		/* Local PSF */
  float		*resi;		/* Map of residuals */
  char		**contextname;	/* Array of context key-names */
  double	*contextoffset;	/* Offset to apply to context data */
  double	*contextscale;	/* Scaling to apply to context data */
  struct poly	*poly;		/* Polynom describing the PSF variations */
  float		pixstep;	/* Mask oversampling (pixel). */
  int		samples_loaded;	/* Number of detections loaded */
  int		samples_accepted;/* Number of detections accepted */
  double	chi2;		/* chi2/d.o.f. */
  float		fwhm;		/* Initial guess of the FWHM */
  int		*pixmask;	/* Pixel mask for local bases */
  float		*basis;		/* Basis vectors */
  int		nbasis;		/* Number of basis vectors */
  int		ndata;		/* Size of the design matrix along data axis */
  int		nsnap;		/* Total number of snapshots */
  int		nmed;		/* Median position amongst snapshots */
  moffatstruct	*moffat;	/* Array of Moffat fits to PSF */
  double	moffat_fwhm_min;
  double	moffat_fwhm;	/* Central Moffat FWHM */
  double	moffat_fwhm_max;
  double	moffat_elongation_min;
  double	moffat_elongation;/* Central Moffat elongation */
  double	moffat_elongation_max;
  double	moffat_beta_min;
  double	moffat_beta;	/* Central Moffat beta */
  double	moffat_beta_max;
  double	moffat_residuals_min;
  double	moffat_residuals;/* Central Moffat residuals */
  double	moffat_residuals_max;
  double	moffat_symresiduals_min;
  double	moffat_symresiduals;/* Symmetry residuals */
  double	moffat_symresiduals_max;
  float		*homo_kernel;		/* PSF homogenization kernel */
  double	homopsf_params[2];	/* Idealised Moffat PSF params*/
  int		homobasis_number;	/* nb of supersampled pixels */
  }	psfstruct;


/*---------------------------------- protos --------------------------------*/
extern void	psf_build(psfstruct *psf, double *pos),
		psf_clip(psfstruct *psf),
		psf_end(psfstruct *psf),
		psf_make(psfstruct *psf, setstruct *set),
		psf_makebasis(psfstruct *psf, setstruct *set,
			basistypenum basis_type,  int nvec),
		psf_makeresi(psfstruct *psf, setstruct *set, int centflag,
			float psf_extraccu),
		psf_makemask(psfstruct *psf, setstruct *set, double chithresh),
		psf_refine(psfstruct *psf, setstruct *set),
		psf_save(psfstruct *psf,  char *filename, int ext, int next);

extern int	psf_pshapelet(float **shape, int w, int h, int nmax,
			double beta),
		psf_readbasis(psfstruct *psf, char *filename, int ext);

extern double	psf_chi2(psfstruct *psf, setstruct *set),
		psf_clean(psfstruct *psf, setstruct *set);

extern psfstruct	*psf_init(char **names, int *group, int ndim,
				int *dim, int ngroup,
				int wpsf, int hpsf, float psfstep,int nsample),
			*psf_load(char *filename);


#endif

