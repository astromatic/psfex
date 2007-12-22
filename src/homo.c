  /*
 				homo.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	PSF homogenisation stuff.
*
*	Last modify:	22/12/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"levmar/lm.h"
#include	"diagnostic.h"
#include	"fft.h"
#include	"homo.h"
#include	"prefs.h"
#include	"poly.h"
#include	"psf.h"
#include	"vignet.h"


/****** psf_homo *******************************************************
PROTO	void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
		int homobasis_number, double homobasis_scale,
		int ext, int next)
PURPOSE	Compute an homogenization kernel based on an idealised PSF.
INPUT	Pointer to the PSF structure.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 22/12/2007
 ***/
void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
		int homobasis_number, double homobasis_scale,
		int ext, int next)
  {
   moffatstruct		*moffat;
   double		dpos[POLY_MAXDIM];
   float		*basis, *bigpsf, *fbigpsf, *bigconv, *cube;
   int			bigsize[2],
			b,i, nbasis, npix,nbigpix;

  npix = psf->size[0]*psf->size[1];

/* Generate real PSF image at center with extra padding for convolution */
  for (i=0; i<psf->poly->ndim; i++)
     dpos[i] = 0.5;
  psf_build(psf, dpos);
  bigsize[0] = psf->size[0]*2;
  bigsize[1] = psf->size[1]*2;
  nbigpix = bigsize[0]*bigsize[1];
  QCALLOC(bigpsf, float, nbigpix);
  vignet_copy(psf->loc, psf->size[0],psf->size[1],
	bigpsf, bigsize[0],bigsize[1], 0,0, VIGNET_CPY);
  fft_shift(bigpsf, bigsize[0], bigsize[1]);
  fbigpsf = fft_rtf(bigpsf, bigsize[0], bigsize[1]);
  free(bigpsf);

/* Create kernel basis */
  nbasis = psf_pshapelet(&basis, psf->size[0],psf->size[1],
	homobasis_number, sqrt(homobasis_number+1.0)*homobasis_scale);

/* Convolve kernel basis vectors with real PSF */
  QMALLOC(bigconv, float, nbigpix);
  fft_init(prefs.nthreads);
  for (b=0; b<nbasis; b++)
    {
    vignet_copy(&basis[b*npix], psf->size[0],psf->size[1],
	bigconv, bigsize[0],bigsize[1], 0,0, VIGNET_CPY);
    fft_conv(bigconv, fbigpsf, bigsize[0], bigsize[1]);
    vignet_copy(bigconv, bigsize[0],bigsize[1],
	&basis[b*npix], psf->size[0],psf->size[1], 0,0, VIGNET_CPY);
    }
  fft_end(prefs.nthreads);
  free(fbigpsf);

/* Create idealized PSF */
  QCALLOC(moffat, moffatstruct, 1);
  moffat->amplitude = 1.0;
  moffat->fwhm_min = moffat->fwhm_max = homopsf_params[0];
  moffat->theta = 0.0;
  moffat->beta = homopsf_params[1];
  psf_moffat(psf, moffat);

  free(basis);
  free(moffat);

  return;
  }


