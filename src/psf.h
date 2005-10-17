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
*	Last modify:	18/12/2002
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PSF_H_
#define _PSF_H_

/*----------------------------- Internal constants --------------------------*/

#define	PSF_FREEDFACTOR	1.0	/* Margin against overfitting */
#define	PSF_NMASKDIM	3	/* Number of dimensions for PSF data */
#define	PSF_MAXSHIFT	3.0	/* Max shift from initial guess (pixels)*/
#define	PSF_MINSHIFT	1e-4	/* Min shift from previous guess (pixels)*/
#define PSF_NITER	40	/* Maximum number of iterations in fit */

/*--------------------------- structure definitions -------------------------*/

typedef struct code
  {
  float		*pc;
  float		**param;
  int		*parammod;
  int		ncode;
  int		nparam;
  }     codestruct;

typedef struct psf
  {
  int		dim;	/* Dimensionality of the tabulated data */
  int		*size;	/* PSF  dimensions */
  int		npix;	/* Total number of involved PSF pixels */
  float		*comp; 	/* Complete pix. data (principal components) */
  float		*loc;	/* Local PSF */
  float		*resi;	/* Map of residuals */
  char		**contextname;	/* Array of context key-names */
  double	*contextoffset;	/* Offset to apply to context data */
  double	*contextscale;	/* Scaling to apply to context data */
  struct poly	*poly;		/* Polynom describing the PSF variations */
  float		pixstep;	/* Mask oversampling (pixel). */
  }	psfstruct;


typedef struct pc
  {
  char		name[MAXCHAR];	/* PC filename */
  int		npc;		/* Number of Principal Components */
  int		dim;	/* Dimensionality of the tabulated data */
  int		*size;	/* PC  dimensions */
  int		npix;	/* Total number of involved PC pixels */
  float		*comp; 	/* Complete pix. data (principal components) */
  double	*mx2,*my2,*mxy;	/* 2nd order moments for each component */
  double	*flux;		/* Flux of each component */
  double	*bt;		/* B/T for each component */
  codestruct	*code;
  }	pcstruct;

typedef struct out_data 
{ int samples_loaded, samples_accepted;
  double chi2;
  float fwhm;
} out_data_struct;


/*---------------------------------- protos --------------------------------*/
extern void	psf_build(psfstruct *psf, double *pos),
		psf_end(psfstruct *psf),
		psf_make(psfstruct *psf, setstruct *set),
		psf_makeresi(psfstruct *psf, setstruct *set, int centflag),
		psf_makemask(psfstruct *psf, setstruct *set, double chithresh),
		psf_refine(psfstruct *psf, setstruct *set, int npsf),
		psf_save(psfstruct *psf, pcstruct *pcc, pcstruct *pc,
			char *filename, int ext, int next, out_data_struct *out_data);

extern double	psf_clean(psfstruct *psf, setstruct *set, int clean_flag);

extern psfstruct	*psf_init(char **names, int *group, int ndim,
				int *dim, int ngroup,
				int wpsf, int hpsf, float psfstep,int nsample),
			*psf_load(char *filename);


extern void	matinv(double *mat, int nmat),
		pc_end(pcstruct *pc);

extern pcstruct	*pc_convolve(pcstruct *pc, psfstruct *psf),
		*pc_load(char *filename),
		*pc_orthogon(pcstruct *pc2, pcstruct *pc, float pixstep);

#endif

