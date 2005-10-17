/*
 				sample.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Type definitions related to samples
*
*	Last modify:	10/07/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _SAMPLE_H_
#define _SAMPLE_H_

/*--------------------------------- constants -------------------------------*/

#define	LSAMPLE_DEFSIZE	1000		/* Sample stacksize at the beginning */

/*--------------------------- structure definitions -------------------------*/

typedef struct sample
  {
  float		*vig;			/* Vignette array */
  float		*vigresi;		/* Chi-map of the PSF-residuals */
  float		*vigweight;		/* Vignette-weight array */
  float		*retina;		/* Retina array */
  float		*retiweight;		/* Retina-weight array */
  float		*pos;			/* Extra parameters array */
  float		norm;			/* Normalisation */
  double	x,y;			/* x,y position estimate in frame */
  float		dx,dy;			/* x,y shift / vignet center */
  float		backnoise2;		/* Variance of the background noise */
  float		gain;			/* conversion factor (e-/ADU) */
  float		chi2;			/* Chi2 of the fit */
  double	*context;		/* Context vector */
  }	samplestruct;

typedef struct set
  {
  char		*head;			/* Table structure */
  struct sample	*sample;		/* Array of samples */
  int		nsample;		/* Number of samples in stack */
  int		nsamplemax;		/* Max number of samples in stack */
  int		*vigsize;		/* Dimensions of vignette frames */
  int		vigdim;			/* Dimensionality of the vignette */
  int		nvig;			/* Number of pixels of the vignette */
  int		retidim;		/* Dimensionality of the retina */
  int		nreti;			/* Number of pixels of the retina */
  int		*retisize;		/* Dimensions of retina frames */
  int		ncontext;		/* Number of contexts */
  char		**contextname;		/* List of context keywords used */
  double	*contextoffset;		/* Offset to apply to context data */
  double	*contextscale;		/* Scaling to apply to context data */
  float		fwhm;			/* FWHM of the PSF core */
  }	setstruct;

/*-------------------------------- protos -----------------------------------*/

samplestruct	*remove_sample(setstruct *set, int isample);

setstruct	*init_set(void),
		*load_samples(char **filename, int ncat, int ext, int next),
		*read_samples(setstruct *set, char *filename,
			float frmin, float frmax, int ext, int next);

void		end_set(setstruct *set),
		free_samples(setstruct *set),
 		malloc_samples(setstruct *set, int nsample),
		make_weights(setstruct *set, samplestruct *sample),
		realloc_samples(setstruct *set, int nsample),
		update_retina(setstruct *set, samplestruct *sample,
			float pixstep);

#endif
