/*
*				psfbasis.c
*
* Manage PSF vector bases.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2019 IAP/CNRS/SorbonneU
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
*	Last modified:		08/05/2019
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"prefs.h"
#include	"context.h"
#include	"misc.h"
#include	"wcs/poly.h"
#include	"psf.h"
#include	"sample.h"
#include	"vignet.h"

static double	psfbasis_laguerre(double x, int p, int q);

/****** psfbasis_make *********************************************************
PROTO	void	psfbasis_make(psfstruct *psf, setstruct *set,
			basistypenum basis_type, int nvec)
PURPOSE	Generate basis vectors for representing the PSF.
INPUT	Pointer to the PSF,
	Pointer to the sample set,
	Basis type,
	Basis number.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 05/10/2010
 ***/
void	psfbasis_make(psfstruct *psf, setstruct *set,
			basistypenum basis_type, int nvec)
  {
  double	xc,yc, x,y, rmax2;
  float		*psforder,*psfordert,*ppix,*basis,
		psfthresh;
  int		*psfmask,
		i, ix,iy, ixmin,iymin,ixmax,iymax,irad, npsf,npix;

  npix = psf->size[0]*psf->size[1];

  switch(basis_type)
    {
    case BASIS_NONE:
      break;
    case BASIS_PIXEL:
/*---- The number of elements is set to the square of the input number */
      npsf = nvec*nvec;
      if (npsf>npix)
        npsf=npix;
/*---- First Select the brightest pixels */
//      NFPRINTF(OUTPUT,"Selecting pixels...");
      psforder = (float *)NULL;		/* To avoid gcc -Wall warnings */
      QMEMCPY(psf->comp, psforder, float, npix);
      for (psfordert=psforder, i=npix; i--; psfordert++)
        *psfordert = fabs(*psfordert);
      fqmedian(psforder, npix);
      psfthresh = psforder[npix-npsf];
      free(psforder);

/*---- Mark pixels which have to be reexamined */
      QCALLOC(psf->pixmask, int, npix);
      psfmask = psf->pixmask;
      npsf = 0;
      irad = (int)((set->vigsize[1]-1)/(2*psf->pixstep));
      iymin = psf->size[1]/2 - irad;
      iymax = psf->size[1]/2 + irad;
      irad = (int)((set->vigsize[0]-1)/(2*psf->pixstep));
      ixmin = psf->size[0]/2 - irad;
      ixmax = psf->size[0]/2 + irad;
      ppix=psf->comp;
      xc = (double)(psf->size[0]/2);
      yc = (double)(psf->size[1]/2);
      y = -yc;
      rmax2 = (psf->size[0]<psf->size[1]? (double)(psf->size[0]/2)
				: (double)(psf->size[1]/2))+0.5;
      rmax2 *= rmax2;
      for (iy=psf->size[1]; iy--; y+=1.0)
        {
        i = iy*psf->size[0];
        x = -xc;
        if (iy>=iymin && iy<=iymax)
          for (ix=psf->size[0]; ix--; i++, x+=1.0)
            if (fabs(ppix[i])>=psfthresh && ix>=ixmin && ix<=ixmax
		&& x*x+y*y<rmax2)
              {
              npsf++;
              psfmask[i] = 1;
              }
        }
      psf->nbasis = npsf;
/*---- Prepare a PSF mask that will contain Dirac peaks only... */
      QCALLOC(psf->basis, float, npsf*npix);
      basis = psf->basis;
      psfmask = psf->pixmask;
      for (i=npix; i--; basis++)
        if (*(psfmask++))
          {
          *basis = 1.0;
          basis += npix;
          }

/*-- Size of the compressed design matrix along the "data" axis */
      psf->ndata = (1+(int)(INTERPW*psf->pixstep))
		*(1+(int)(INTERPW*psf->pixstep))+1;
      break;

    case BASIS_GAUSS_LAGUERRE:
      psf->nbasis = psfbasis_pshapelet(&psf->basis, psf->size[0],psf->size[1],
		nvec, sqrt(nvec+1.0)*prefs.basis_scale);
      break;
    case BASIS_FILE:
      psf->nbasis = psfbasis_read(psf, prefs.basis_name, 0);
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: unknown PSF vector basis in ",
			"psfbasis_make()");
    }

  return;
  }


/****** psfbasis_laguerre ****************************************************
PROTO	double	psfbasis_laguerre(double x, int p, int q)
PURPOSE	Return Laguerre polynomial value.
INPUT	x,
	p,
	q.
OUTPUT  Value of the Laguerre polynomial.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2007
 ***/
static double	psfbasis_laguerre(double x, int p, int q)
  {
   double	dn,dq, lpm1,lpm2, l;
   int		n;

  dq = q - 1.0;
  if (p==0)
    return 1.0;
  else if (p==1)
    return (2.0 - x + dq);
  else
    {
    l = 0.0;
    lpm2 = 1.0;
    lpm1 = 2.0 - x + dq;
    dn = 2.0;
    for (n=p-1; n--; dn+=1.0)
      {
      l = (2.0+(dq-x)/dn)*lpm1 - (1.0+dq/dn)*lpm2;
      lpm2 = lpm1;
      lpm1 = l;
      }
    }

  return l;
  }


/****** psfbasis_pshapelet ***************************************************
PROTO	int psfbasis_pshapelet(float **shape, int w, int h, int nmax,
		double beta)
PURPOSE	Compute Polar shapelet basis set.
INPUT	Pointer to the array of image vectors (which will be allocated),
	Image vector width,
	Image vector height,
	Shapelet n_max,
	beta parameter.
OUTPUT  Total number of image vectors generated.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 08/05/2019
 ***/
int psfbasis_pshapelet(float **basis, int w, int h, int nmax, double beta)
  {
   char		str[128];
   double	*fr2,*fr2t,*fexpr2,*fexpr2t,*ftheta,*fthetat,
		dm,fac, xc,yc, x,y, x1,y1, r2,rmax2, invbeta2, val,
		ostep,ostep2,odx;
   float	*basist;
   int		i,j,k, m,n,p, kmax,hnmm, ix,iy, idx,idy;

  kmax = (nmax+1)*(nmax+2)/2;

  invbeta2 = 1.0/(beta*beta);
  ostep = 1.0/(GAUSS_LAG_OSAMP);
  ostep2 = ostep*ostep;
  odx = 0.5*(ostep - 1.0);
  xc =(double)(w/2);
  yc = (double)(h/2);
  rmax2 = (xc<yc? xc: yc);
  rmax2 *= rmax2*invbeta2;

/* Precompute some slow functions */
  QMALLOC(fr2, double, w*h*GAUSS_LAG_OSAMP*GAUSS_LAG_OSAMP);
  QMALLOC(fexpr2, double, w*h*GAUSS_LAG_OSAMP*GAUSS_LAG_OSAMP);
  QMALLOC(ftheta, double, w*h*GAUSS_LAG_OSAMP*GAUSS_LAG_OSAMP);
  fr2t = fr2;
  fexpr2t = fexpr2;
  fthetat = ftheta;
  y = odx - yc;
  for (iy=h; iy--; y+=1.0)
    {
    x = odx - xc;
    for (ix=w; ix--; x+=1.0)
      {
      y1 = y;
      for (idy=GAUSS_LAG_OSAMP; idy--; y1+=ostep)
        {
        x1 = x;
        for (idx=GAUSS_LAG_OSAMP; idx--; x1+=ostep)
          {
          *(fr2t++) = r2 = (x1*x1+y1*y1)*invbeta2;
          *(fexpr2t++) = exp(-r2/2.0);
          *(fthetat++) = atan2(y1,x1);
          }
        }
      }
    }

  QCALLOC(*basis, float, w*h*kmax);
  basist = *basis;
  k=1;
  for (n=0; n<=nmax; n++)
    {
    for (m=n%2; m<=n; m+=2)
      {
      sprintf(str, "Generating basis vector #%d/%d", k++, kmax);
//      NFPRINTF(OUTPUT, str);
      dm = (double)m;
/*---- Compute ((n+m)/2)!/((n-m)/2)! */
      hnmm = (n-m)/2;
      fac = 1.0;
      for (p=(n+m)/2; p>=hnmm; p--)
        if (p)
          fac *= (double)p;
      fac = sqrt(1.0/(PI*fac))/beta;
      if ((hnmm%2))
        fac = -fac;
      fr2t = fr2;
      fexpr2t = fexpr2;
      fthetat = ftheta;
      for (i=w*h; i--;)
        {
        val = 0.0;
        for (j=GAUSS_LAG_OSAMP*GAUSS_LAG_OSAMP; j--; fr2t++)
          val += fac*pow(*fr2t, dm/2.0)*psfbasis_laguerre(*fr2t, hnmm, m)
		**(fexpr2t++)*cos(dm**(fthetat++));
        *(basist++) = val*ostep2;
        }
      if (m!=0)
        {
        fr2t = fr2;
        fexpr2t = fexpr2;
        fthetat = ftheta;
        for (i=w*h; i--;)
          {
          val = 0.0;
          for (j=GAUSS_LAG_OSAMP*GAUSS_LAG_OSAMP; j--; fr2t++)
            val += fac*pow(*fr2t, dm/2.0)*psfbasis_laguerre(*fr2t, hnmm, m)
		**(fexpr2t++)*sin(dm**(fthetat++));
          *(basist++) = val*ostep2;
          }
        k++;
        }
      }
    }

  free(fr2);
  free(fexpr2);
  free(ftheta);

  return kmax;
  }


/****** psfbasis_read ********************************************************
PROTO   int psfbasis_read(psfstruct *psf, char *filename, int ext)
PURPOSE Read a set of basis functions for the PSF from a 3D FITS-file.
INPUT   Pointer to the PSF structure,
	FITS filename,
	Extension number.
OUTPUT  Number of basis vectors read.
NOTES   The maximum degrees and number of dimensions allowed are set in poly.h.
AUTHOR  E. Bertin (IAP)
VERSION 13/11/2007
 ***/
int	psfbasis_read(psfstruct *psf, char *filename, int ext)
  {
   catstruct	*cat;
   tabstruct	*tab, *firstab;
   PIXTYPE	*pixin;
   int		n, next, extp1, ntabp1, npixin,npixout,ncomp;

/*-- Read input FITS file */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: No such catalog: ", filename);
/* Go to the right extension */
  tab = cat->tab;
  ntabp1 = cat->ntab+1;
  firstab = NULL;
  extp1 = ext+1;
  for (next=0; ntabp1-- && next<extp1; tab = tab->nexttab)
    if (tab->naxis>=2)
      {
      if (!next)
        firstab = tab;
      next++;
      }
  if (!ntabp1)
    {
    if (!next)
      error(EXIT_FAILURE, "No image data in ", filename);
    if (next>extp1)
      warning("Not enough extensions, using only 1st datacube of ",
		filename);
    }

  tab = tab->prevtab;
  npixin = tab->naxisn[0]*tab->naxisn[1];
  npixout = psf->size[0]*psf->size[1];
  QMALLOC(pixin, PIXTYPE, npixin);
  ncomp = tab->tabsize/tab->bytepix/npixin;
  QMALLOC(psf->basis, float, ncomp*npixout);
  QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, tab->cat->filename);
  for (n=0; n<ncomp; n++)
    {
    read_body(tab, pixin, npixin);
    vignet_copy(pixin, tab->naxisn[0], tab->naxisn[1],
		&psf->basis[n*npixout], psf->size[0], psf->size[1], 0, 0,
		VIGNET_CPY);
    }
  free(pixin);
  free_cat(&cat, 1);

  return ncomp;
  }


