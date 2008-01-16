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
*	Last modify:	15/01/2008
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
#include	"diagnostic.h"
#include	"fft.h"
#include	"homo.h"
#include	"prefs.h"
#include	"poly.h"
#include	"psf.h"
#include	"vignet.h"
#include	ATLAS_LAPACK_H


/****** psf_homo *******************************************************
PROTO	void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
		int homobasis_number, double homobasis_scale,
		int ext, int next)
PURPOSE	Compute an homogenization kernel based on an idealised PSF.
INPUT	Pointer to the PSF structure.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 15/01/2008
 ***/
void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
		int homobasis_number, double homobasis_scale,
		int ext, int next)
  {
   moffatstruct		*moffat;
   polystruct		*poly;
   double		dpos[POLY_MAXDIM],
			*amat, *bmat, *contcross, *coeff,
			dstep,dstart, dval;
   float		*basis,*basisc,*basis1,*basis2, *bigpsf, *fbigpsf,
			*bigconv, *mofpix,*mofpixt,
			*kernorm,*kernel,*kernelt,
			a,b;
   int			bigsize[2],
			d,i,j,n,p, nt, nbasis, npix,nbigpix, ndim, ncoeff;

  NFPRINTF(OUTPUT,"Computing the PSF homogenization kernel...");

  npix = psf->size[0]*psf->size[1];
  poly = psf->poly;
  ndim = poly->ndim;
  coeff = poly->coeff;
  ncoeff = poly->ncoeff;

/* Computing context cross-products */
  QCALLOC(contcross, double, ncoeff*ncoeff);
  nt = 1;
  for (d=0; d<ndim; d++)
    nt *= HOMO_NSNAP;
  dstep = 1.0/HOMO_NSNAP;
  dstart = (1.0-dstep)/2.0;
  for (d=0; d<ndim; d++)
    dpos[d] = -dstart;
  for (n=0; n<nt; n++)
    {
    poly_func(poly, dpos);
    for (j=0; j<ncoeff; j++)
      for (i=j; i<ncoeff; i++)
        contcross[j*ncoeff+i] += coeff[j]*coeff[i];
    for (d=0; d<ndim; d++)
      if (dpos[d]<dstart-0.01)
        {
        dpos[d] += dstep;
        break;
        }
      else
        dpos[d] = -dstart;
    }

  free(contcross);

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

/* Create kernel basis */
  nbasis = psf_pshapelet(&basis, psf->size[0],psf->size[1],
	homobasis_number, sqrt(homobasis_number+1.0)*homobasis_scale);
  basisc = NULL;		/* to avoid gcc -Wall warnings */
  QMEMCPY(basis, basisc, float, nbasis*npix);

/* Convolve kernel basis vectors with real PSF */
  QMALLOC(bigconv, float, nbigpix);
  fft_init(prefs.nthreads);
  for (i=0; i<nbasis; i++)
    {
    vignet_copy(&basis[i*npix], psf->size[0],psf->size[1],
	bigconv, bigsize[0],bigsize[1], 0,0, VIGNET_CPY);
    fft_conv(bigconv, fbigpsf, bigsize[0], bigsize[1]);
    vignet_copy(bigconv, bigsize[0],bigsize[1],
	&basisc[i*npix], psf->size[0],psf->size[1], 0,0, VIGNET_CPY);
    }
  fft_end(prefs.nthreads);
  free(bigpsf);

/* Create idealized PSF */
  QCALLOC(moffat, moffatstruct, 1);
  moffat->xc[0] = (double)(psf->size[0]/2);
  moffat->xc[1] = (double)(psf->size[1]/2);
  moffat->amplitude = 1.0;
  moffat->fwhm_min = moffat->fwhm_max = homopsf_params[0];
  moffat->theta = 0.0;
  moffat->beta = homopsf_params[1];
  psf_moffat(psf, moffat);
  mofpix = psf->loc;
  free(moffat);

/* Compute the normal equation matrix */
  QCALLOC(amat, double, nbasis*nbasis);
  QCALLOC(bmat, double, nbasis);
  for (j=0;j<nbasis; j++)
    for (i=0; i<nbasis; i++)
      {
      basis1 = basisc + j*npix;
      basis2 = basisc + i*npix;
      a = 0.0;
      for (p=npix; p--;)
        a += *(basis1++)**(basis2++);
      amat[j*nbasis+i] = a;
      }
  for (j=0; j<nbasis; j++)
    {
    basis1 = basisc + j*npix;
    mofpixt = mofpix;
    b = 0.0;
    for (p=npix; p--;)
      b += *(basis1++)**(mofpixt++);
    bmat[j] = b;
    }
  free(basisc);

  clapack_dpotrf(CblasRowMajor, CblasUpper, nbasis, amat, nbasis);
  clapack_dpotrs(CblasRowMajor, CblasUpper, nbasis, 1, amat, nbasis,
	bmat, nbasis);

  QCALLOC(kernel, float, npix);
  for (j=0; j<nbasis; j++)
    {
    basis1 = basis + j*npix;
    b = bmat[j];
    kernelt = kernel;
    for (p=npix; p--;)
      *(kernelt++) += b**(basis1++);
    }

  free(amat);
  free(bmat);
  free(basis);

/* Normalize kernel (with respect to continuum) */
/*
  vignet_copy(kernel, psf->size[0],psf->size[1],
	bigconv, bigsize[0],bigsize[1], 0,0, VIGNET_CPY);
  fft_conv(bigconv, fbigpsf, bigsize[0], bigsize[1]);
  QMALLOC(kernorm, float, npix);
  vignet_copy(bigconv, bigsize[0],bigsize[1],
	kernorm, psf->size[0],psf->size[1], 0,0, VIGNET_CPY);
*/
  a = 0.0;
  kernelt = kernel;
  for (p=npix; p--;)
    a += *(kernelt++);
  if (a>0.0)
    {
    a = 1.0/a;
    kernelt = kernel;
    for (p=npix; p--;)
      *(kernelt++) *= a;
    }

  free(fbigpsf);
  free(bigconv);
//  free(kernorm);

/* Save homogenization kernel */
  psf->homo_kernel = kernel;
  psf->homopsf_params[0] = homopsf_params[0];
  psf->homopsf_params[1] = homopsf_params[1];
  psf->homobasis_number = homobasis_number;
  psf_savehomo(psf, filename, ext, next);

  return;
  }


/****** psf_savehomo **********************************************************
PROTO   void	psf_savehomo(psfstruct *psf, char *filename, int ext, int next)
PURPOSE Save the PSF homogenization kernel data as a FITS file.
INPUT   Pointer to the PSF structure,
	Filename,
	Extension number,
	Number of extensions.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 26/12/2007
 ***/
void	psf_savehomo(psfstruct *psf, char *filename, int ext, int next)
  {
   static catstruct	*cat;
   tabstruct		*tab;
   char			str[88];
   int			i, temp;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  if (!ext)
    {
    cat = new_cat(1);
    init_cat(cat);
    strcpy(cat->filename, filename);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);
    if (next>1)
      save_tab(cat, cat->tab);
    }
  tab = new_tab("HOMO_DATA");
  addkeywordto_head(tab, "POLNAXIS", "Number of context parameters");
  fitswrite(tab->headbuf, "POLNAXIS", &psf->poly->ndim, H_INT, T_LONG);
  for (i=0; i<psf->poly->ndim; i++)
    {
    sprintf(str, "POLGRP%1d", i+1);
    addkeywordto_head(tab, str, "Polynom group for this context parameter");
    temp = psf->poly->group[i]+1;
    fitswrite(tab->headbuf, str, &temp, H_INT, T_LONG);
    sprintf(str, "POLNAME%1d", i+1);
    addkeywordto_head(tab, str, "Name of this context parameter");
    fitswrite(tab->headbuf, str, psf->contextname[i], H_STRING, T_STRING);
    sprintf(str, "POLZERO%1d", i+1);
    addkeywordto_head(tab, str, "Offset value for this context parameter");
    fitswrite(tab->headbuf, str, &psf->contextoffset[i], H_EXPO, T_DOUBLE);
    sprintf(str, "POLSCAL%1d", i+1);
    addkeywordto_head(tab, str, "Scale value for this context parameter");
    fitswrite(tab->headbuf, str, &psf->contextscale[i], H_EXPO, T_DOUBLE);
    }

  addkeywordto_head(tab, "POLNGRP", "Number of context groups");
  fitswrite(tab->headbuf, "POLNGRP", &psf->poly->ngroup, H_INT, T_LONG);
  for (i=0; i<psf->poly->ngroup; i++)
    {
    sprintf(str, "POLDEG%1d", i+1);
    addkeywordto_head(tab, str, "Polynom degree for this context group");
    fitswrite(tab->headbuf, str, &psf->poly->degree[i], H_INT, T_LONG);
    }

/* Add and write important scalars as FITS keywords */
  /* -- FM -- : write fwhm too */
  addkeywordto_head(tab, "PSF_FWHM", "PSF FWHM");
  fitswrite(tab->headbuf, "PSF_FWHM", &psf->homopsf_params[0],
	H_FLOAT,T_DOUBLE);
  addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data");
  fitswrite(tab->headbuf, "PSF_SAMP", &psf->pixstep, H_FLOAT, T_FLOAT);
  tab->bitpix = BP_FLOAT;
  tab->bytepix = t_size[T_FLOAT];
  tab->naxis = 2;
  tab->naxisn[0] = psf->size[0];
  tab->naxisn[1] = psf->size[1];
  tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
  tab->bodybuf = (char *)psf->homo_kernel; 
  if (next == 1)
    prim_head(tab);
  fitswrite(tab->headbuf, "XTENSION", "IMAGE   ", H_STRING, T_STRING);

  save_tab(cat, tab);
/* But don't touch my arrays!! */
  tab->bodybuf = NULL;
  free_tab(tab);

  if (ext==next-1)
    free_cat(&cat , 1);

  return;
  }

