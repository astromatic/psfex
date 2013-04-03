/*
*				homo.c
*
* PSF homogenisation.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2008-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		20/11/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
#include	"wcs/poly.h"
#include	"psf.h"
#include	"vignet.h"

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
//#define MATSTORAGE_PACKED 1
#endif

/****** psf_homo *******************************************************
PROTO	void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
		int homobasis_number, double homobasis_scale,
		int ext, int next)
PURPOSE	Compute an homogenization kernel based on an idealised PSF.
INPUT	Pointer to the PSF structure.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 10/07/2012
 ***/
void	psf_homo(psfstruct *psf, char *filename, double *homopsf_params,
		int homobasis_number, double homobasis_scale,
		int ext, int next)
  {
   moffatstruct		*moffat;
   polystruct		*poly;
   double		dpos[POLY_MAXDIM],
			*amat,*amatt, *bmat,*bmatt, *cross,*tcross, *coeff,
			dstep,dstart, dval;
   float		*basis,*basisc,*basis1,*basis2, *bigbasis, *fbigbasis,
			*bigconv, *target,
			*kernel,*kernelt,
			a;
   int			bigsize[2],
			c,c1,c2,c3,c4,d,f1,f2,i,j,j1,j2,n,p,
			nt, npix,nbigpix, ndim, nbasis,ncoeff,nfree;

//  NFPRINTF(OUTPUT,"Computing the PSF homogenization kernel...");

  npix = psf->size[0]*psf->size[1];
  poly = psf->poly;
  ndim = poly->ndim;
  coeff = poly->basis;
  ncoeff = poly->ncoeff;

/* Create kernel basis */
  nbasis = psf_pshapelet(&basis, psf->size[0],psf->size[1],
	homobasis_number, sqrt(homobasis_number+1.0)*homobasis_scale);
  nfree = nbasis*ncoeff;

/* Create idealized PSF */
  QCALLOC(moffat, moffatstruct, 1);
  moffat->xc[0] = (double)(psf->size[0]/2);
  moffat->xc[1] = (double)(psf->size[1]/2);
  moffat->amplitude = 1.0;
  moffat->fwhm_min = moffat->fwhm_max = homopsf_params[0];
  moffat->theta = 0.0;
  moffat->beta = homopsf_params[1];
  moffat->nsubpix = 1;
  psf_moffat(psf, moffat);
  target = NULL;		/* to avoid gcc -Wall warnings */
  QMEMCPY(psf->loc, target, float, npix);
  free(moffat);

/* Prepare a padded space for PSF */
  bigsize[0] = psf->size[0]*2;
  bigsize[1] = psf->size[1]*2;
  nbigpix = bigsize[0]*bigsize[1];
  QMALLOC(bigbasis, float, nbigpix);
/* bigconv will receive the result of the convolution */
  QMALLOC(bigconv, float, nbigpix);

/* Convolve kernel basis vectors with PSF components and compute X-products*/
  QMALLOC(basisc, float, nfree*npix);
  QMALLOC(cross, double, nfree*nfree);
  QMALLOC(tcross, double, nfree);
  fft_init(prefs.nthreads);
  f1 = 0;
  for (i=0; i<nbasis; i++)
    {
    vignet_copy(&basis[i*npix], psf->size[0],psf->size[1],
	bigbasis, bigsize[0],bigsize[1], 0,0, VIGNET_CPY);
    fft_shift(bigbasis, bigsize[0], bigsize[1]);
    fbigbasis = fft_rtf(bigbasis, bigsize[0], bigsize[1]);
    for (c=0; c<ncoeff; c++, f1++)
      {
      vignet_copy(&psf->comp[c*npix], psf->size[0],psf->size[1],
	bigconv, bigsize[0],bigsize[1], 0,0, VIGNET_CPY);
      fft_conv(bigconv, fbigbasis, bigsize[0], bigsize[1]);
      vignet_copy(bigconv, bigsize[0],bigsize[1],
	&basisc[f1*npix], psf->size[0],psf->size[1], 0,0, VIGNET_CPY);
      basis2 = basisc;
      for (f2=0; f2<=f1; f2++)
        {
        basis1 = basisc + f1*npix;
        dval = 0.0;
        for (p=npix; p--;)
          dval += *(basis1++)**(basis2++);
        cross[f1+f2*nfree] = cross[f2+f1*nfree] = dval;
        }
/*---- Cross-product with target PSF */
      basis1 = basisc + f1*npix;
      basis2 = target;
      dval = 0.0;
      for (p=npix; p--;)
        dval += *(basis1++)**(basis2++);
      tcross[f1] = dval;
      }
    free(fbigbasis);
    }

  fft_end(prefs.nthreads);
  free(basisc);
  free(bigbasis);
  free(bigconv);
  free(target);

  QCALLOC(amat, double, nfree*nfree);
  QCALLOC(bmat, double, nfree);
/* Proceed through realizations of the PSF throughout the grid of parameters */
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
    amatt = amat;
    bmatt = bmat;
    for (j1=0; j1<nbasis; j1++)
      {
      for (c1=0; c1<ncoeff; c1++)
        {
        for (j2=0; j2<nbasis; j2++)
          {
          for (c2=0; c2<ncoeff; c2++)
            {
            dval = 0.0;
            for (c3=0; c3<ncoeff; c3++)
              {
              f1 = j1*ncoeff + c3;
              for (c4=0; c4<ncoeff; c4++)
                dval += coeff[c3]*coeff[c4]*cross[f1*nfree+j2*ncoeff + c4];
              }
            *(amatt++) += coeff[c1]*coeff[c2]*dval;
            }
          }
        dval = 0.0;
        for (c3=0; c3<ncoeff; c3++)
          dval += coeff[c3]*tcross[j1*ncoeff+c3];
        *(bmatt++) += coeff[c1]*dval;
        }
      }

    for (d=0; d<ndim; d++)
      if (dpos[d]<dstart-0.01)
        {
        dpos[d] += dstep;
        break;
        }
      else
        dpos[d] = -dstart;
    }

  free(cross);
  free(tcross);

#if defined(HAVE_LAPACKE)
 #ifdef MATSTORAGE_PACKED
  if (LAPACKE_dppsv(LAPACK_COL_MAJOR,'L',nfree,1,amat,bmat,nfree) != 0)
 #else
  if (LAPACKE_dposv(LAPACK_COL_MAJOR,'L',nfree,1,amat,nfree,bmat,nfree) != 0)
 #endif
#else
  if (clapack_dposv(CblasRowMajor,CblasUpper,nfree,1,amat,nfree,bmat,nfree)!=0)
#endif
    warning("Not a positive definite matrix", " in homogenization solver");

  QCALLOC(kernel, float, npix*ncoeff);
  bmatt = bmat;
  for (j=0; j<nbasis; j++)
    {
    kernelt = kernel;
    for (c=ncoeff; c--;)
      {
      basis1 = basis + j*npix;
      dval = *(bmatt++);
      for (p=npix; p--;)
        *(kernelt++) += dval**(basis1++);
      }
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
  dval = 0.0;
  kernelt = kernel;
  for (p=npix; p--;)
    dval += *(kernelt++);
  if (dval>0.0)
    {
    a = 1.0/dval;
    kernelt = kernel;
    for (p=npix*ncoeff; p--;)
      *(kernelt++) *= a;
    }

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
VERSION 23/03/2011
 ***/
void	psf_savehomo(psfstruct *psf, char *filename, int ext, int next)
  {
   static catstruct	*cat;
   tabstruct		*tab;
   polystruct		*poly;
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

  poly = psf->poly;
  tab = new_tab("HOMO_DATA");
  addkeywordto_head(tab, "POLNAXIS", "Number of context parameters");
  fitswrite(tab->headbuf, "POLNAXIS", &poly->ndim, H_INT, T_LONG);
  for (i=0; i<poly->ndim; i++)
    {
    sprintf(str, "POLGRP%1d", i+1);
    addkeywordto_head(tab, str, "Polynom group for this context parameter");
    temp = poly->group[i]+1;
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
  fitswrite(tab->headbuf, "POLNGRP", &poly->ngroup, H_INT, T_LONG);
  for (i=0; i<poly->ngroup; i++)
    {
    sprintf(str, "POLDEG%1d", i+1);
    addkeywordto_head(tab, str, "Polynom degree for this context group");
    fitswrite(tab->headbuf, str, &poly->degree[i], H_INT, T_LONG);
    }

/* Add and write important scalars as FITS keywords */
  /* -- FM -- : write fwhm too */
  addkeywordto_head(tab, "PSF_FWHM", "FWHM of target PSF");
  fitswrite(tab->headbuf, "PSF_FWHM", &psf->homopsf_params[0],
	H_FLOAT,T_DOUBLE);
  addkeywordto_head(tab, "PSF_BETA", "Moffat Beta of target PSF");
  fitswrite(tab->headbuf, "PSF_BETA", &psf->homopsf_params[1],
	H_FLOAT,T_DOUBLE);
  addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data");
  fitswrite(tab->headbuf, "PSF_SAMP", &psf->pixstep, H_FLOAT, T_FLOAT);
  tab->bitpix = BP_FLOAT;
  tab->bytepix = t_size[T_FLOAT];
  if (poly->ncoeff>1)
    {
    tab->naxis = 3;
    QREALLOC(tab->naxisn, int, tab->naxis);
    tab->naxisn[0] = psf->size[0];
    tab->naxisn[1] = psf->size[1];
    tab->naxisn[2] = poly->ncoeff;
    tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1]*tab->naxisn[2];
    }
  else
    {
    tab->naxis = 2;
    tab->naxisn[0] = psf->size[0];
    tab->naxisn[1] = psf->size[1];
    tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
    }
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

