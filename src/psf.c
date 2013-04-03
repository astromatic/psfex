/*
*				psf.c
*
* PSF management and modelling.
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
*	Last modified:		02/04/2013
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
#include	"prefs.h"
#include	"context.h"
#include	"misc.h"
#include	"wcs/poly.h"
#include	"psf.h"
#include	"sample.h"
#include	"vignet.h"

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
//#define MATSTORAGE_PACKED 1
#endif

static double	psf_laguerre(double x, int p, int q);

/****** psf_clean *************************************************************
PROTO	double	psf_clean(psfstruct *psf, setstruct *set)
PURPOSE	Filter out PSF candidates
INPUT	Pointer to the PSF,
	Pointer to the sample set,
	PSF accuracy.
OUTPUT	Reduced chi2.
NOTES	-
AUTHOR	E. Bertin (IAP)
VERSION	19/02/2009
 ***/
double	psf_clean(psfstruct *psf, setstruct *set, double prof_accuracy)
  {
#define	EPS	(1e-4)  /* a small number */
   samplestruct	*sample;
   double	chi2,chimean,chivar,chisig,chisig1,chival, locut,hicut;
   float	*chi, *chit,*chit2,
		chimed, chi2max;
   int		i, n, nsample;

/* First compute residuals for each sample (chi^2) */
//  NFPRINTF(OUTPUT,"Computing residuals...");
  psf_makeresi(psf, set, prefs.recenter_flag, prof_accuracy);

/* Store the chi's (sqrt(chi2) pdf close to Gaussian) */
//  NFPRINTF(OUTPUT,"Computing Chi2 statistics...");
  nsample = set->nsample;
  QMALLOC(chi, float, nsample);
  chit = chi;
  for (sample=set->sample, n=nsample; n--; sample++)
    *(chit++) = (float)sqrt(sample->chi2);
/* Produce k-sigma-clipped statistiscs */
  locut = -BIG;
  hicut = BIG;
  chisig = BIG;
  chisig1 = 1.0;
  chivar = 0.0;
  for (i=100; i-- && chisig>=0.1 && fabs(chisig/chisig1-1.0)>EPS;)
    {
    chisig1 = chisig;
    chimed = fast_median(chi, nsample);
    chimean = chivar = 0.0;
    chit2 = chit = chi;
    for (n=nsample; n--;)
      {
      chival = *(chit++);
      if (chival>locut && chival<hicut)
        {
        chimean += (*(chit2++) = chival);
        chivar += chival*chival;
        }
      else
        nsample--;
      }

    chimean /= (double)nsample;
    chisig = sqrt((chivar-chimean*chimean*nsample)/(nsample-(nsample>1?1:0)));
    locut = chimed - 3.0*chisig;
    hicut = chimed + 3.0*chisig;
    }

  free(chi);
/*
  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "<Chi2/dof> = %.3f\n",chivar/(nsample-(nsample>1?1:0)));
*/
  chi2 = chivar/(nsample-(nsample>1?1:0));

/* Clip outliers */
//  NFPRINTF(OUTPUT,"Filtering PSF-candidates...");
  chi2max = (float)hicut;
  chi2max *= chi2max;
  nsample=set->nsample;
  for (sample=set->sample, n=0; n<nsample;)
  if ((sample++)->chi2>chi2max)
    {
    sample=remove_sample(set, n);
    nsample--;
    }
  else
    n++;

  return chi2;
#undef EPS
  }


/****** psf_chi2 **************************************************************
PROTO	double	psf_chi2(psfstruct *psf, setstruct *set)
PURPOSE	Return the reduced chi2 of PSF-fitting.
INPUT	Pointer to the PSF,
	Pointer to the sample set.
OUTPUT	Reduced chi2.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/02/2009
 ***/
double	psf_chi2(psfstruct *psf, setstruct *set)
  {
   samplestruct	*sample;
   double	chi2;
   int		n, nsample;

/* First compute residuals for each sample (chi^2) */
//  NFPRINTF(OUTPUT,"Computing residuals...");
  psf_makeresi(psf, set, prefs.recenter_flag, prefs.prof_accuracy);

/* Store the chi's (sqrt(chi2) pdf close to gaussian) */
//  NFPRINTF(OUTPUT,"Computing Chi2 statistics...");
  nsample = set->nsample;
  chi2 = 0.0;
  for (sample=set->sample, n=nsample; n--; sample++)
    chi2 += sample->chi2;
  chi2 /= (nsample-(nsample>1?1:0));

  return chi2;
  }


/****** psf_clip **************************************************************
PROTO	void	psf_clip(psfstruct *psf)
PURPOSE	Apply soft-clipping to the PSF model using circular boundaries.
INPUT	Pointer to the PSF.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	12/11/2007
 ***/
void	psf_clip(psfstruct *psf)
  {
   double	xc,yc, x,y, r2,rmin2,rmax2,dr2;
   float	*pix;
   int		p, ix,iy, npsf;

  xc = (double)(psf->size[0]/2);
  yc = (double)(psf->size[1]/2);
  rmax2 = (psf->size[0]<psf->size[1]? (double)(psf->size[0]/2)
				: (double)(psf->size[1]/2))+0.5;
  dr2 = (psf->fwhm / psf->pixstep);
  if (dr2<1.0)
    dr2 = 1.0;
  if (dr2 >= rmax2)
    dr2 = rmax2/2.0;
  rmin2 = rmax2 - dr2;
  rmin2 *= rmin2;
  rmax2 *= rmax2;
  dr2 = rmax2 - rmin2;
  npsf = psf->poly->ncoeff;
  pix = psf->comp;
  for (p=npsf; p--;)
    {
    y = -yc;
    for (iy=psf->size[1]; iy--; y+=1.0)
      {
      x = -xc;
      for (ix=psf->size[0]; ix--; x+=1.0, pix++)
        if ((r2=x*x+y*y)>=rmin2)
          {
          if (r2<rmax2)
            *pix *= (rmax2-r2) / dr2;
          else
            *pix = 0.0;
          }
      }
    }

  return;
  }


/****** psf_init **************************************************************
PROTO	psfstruct *psf_init(contextstruct *context, int *size,
			float psfstep, float *pixsize, int nsample)
PURPOSE	Allocate and initialize a PSF structure.
INPUT	Pointer to context structure,
	PSF model image size (pixels),
	PSF pixel step,
	array of 2 effective pixel sizes (along x and along y),
	number of samples.
OUTPUT  psfstruct pointer.
NOTES   The maximum degrees and number of dimensions allowed are set in poly.h.
AUTHOR  E. Bertin (IAP)
VERSION 19/02/2013
 ***/
psfstruct	*psf_init(contextstruct *context, int *size,
			float psfstep, float *pixsize, int nsample)
  {
   psfstruct	*psf;
   static char	str[MAXCHAR];
   char		**names2, **names2t;
   double	psfelemdens;
   int		*group2, *dim2,
		d, ndim,ndim2,ngroup2, npix, nsnap;

/* Allocate memory for the PSF structure itself */
  QCALLOC(psf, psfstruct, 1);
  psf->dim = PSF_NMASKDIM;	/* This is constant */
  QMALLOC(psf->size, int, psf->dim);

/* The polynom */
  names2 = NULL;
  group2 = dim2 = NULL;
  ndim2 = ndim = context->ncontext;
  if (ndim)
    {
    QMEMCPY(context->group, group2, int, ndim);
    QMEMCPY(context->name, names2, char *, ndim);
    }
  if ((ngroup2=context->ngroup))
    QMEMCPY(context->degree, dim2, int, context->ngroup);

  psf->poly = poly_init(group2, ndim2, dim2, ngroup2);

/* Add additional constraint for supersampled PSFs */
  psfelemdens = (psfstep>0.0 && psfstep<1.0)? 1.0/psfstep*psfstep : 1.0;
/*-- Compute the maximum advised number of degrees of freedom */
  if (ngroup2)
    while (psf->poly->ncoeff*psfelemdens*PSF_FREEDFACTOR+0.499 > nsample)
      {
      poly_end(psf->poly);
      if (ngroup2)
        {
/*------ If still too many degrees of freedom, try to lower degrees */
        d=ngroup2%10;
        sprintf(str, "%d%s", ngroup2, d==1?"st":(d==2?"nd":(d==3?"rd":"th")));
        if (!(--dim2[ngroup2-1]))
          {
/*------ If degree is 0, just remove all the group components */
          for (d=0; d<ndim2; d++)
            if (group2[d]==ngroup2 && d!=(--ndim2))
              {
              names2[d]=names2[ndim2];
              group2[d]=group2[ndim2];
              }
          warning(str, " context group removed (not enough samples)");
          if (!(--ngroup2))
            ndim2 = 0;
          }
        else
          warning(str, " context group-degree lowered (not enough samples)");
        psf->poly = poly_init(group2, ndim2, dim2, ngroup2);
        }
      if (!ngroup2)
        break;	/* No sample at all!*/
      }

  psf->pixstep = psfstep;
  psf->pixsize[0] = pixsize[0];
  psf->pixsize[1] = pixsize[1];
  psf->npix = psf->size[0] = size[0];
  psf->npix *= (psf->size[1] = size[1]);
  psf->npix *= (psf->size[2] = psf->poly->ncoeff);
  QCALLOC(psf->comp, float, psf->npix);
  npix = psf->size[0]*psf->size[1];
  QCALLOC(psf->loc, float, npix);
  QCALLOC(psf->resi, float, npix);

/* Context arrays */
  nsnap = 1;
  psf->nsnap = prefs.context_nsnap;
  psf->cx = psf->cy = -1;
  if (ndim2)
    {
    QMALLOC(psf->contextoffset, double, ndim2);
    QMALLOC(psf->contextscale, double, ndim2);
    QMALLOC(psf->contextname, char *, ndim2);
    for (names2t=names2, d=0; d<ndim2; d++)
      {
      nsnap *= psf->nsnap;
      QMALLOC(psf->contextname[d], char, 80);
      strcpy(psf->contextname[d], *(names2t++));
/*---- Identify first spatial coordinates among contexts */
      if (!strcmp(psf->contextname[d], "X_IMAGE")
		|| !strcmp(psf->contextname[d], "XWIN_IMAGE")
		|| !strcmp(psf->contextname[d], "XPSF_IMAGE")
		|| !strcmp(psf->contextname[d], "XMODEL_IMAGE")
		|| !strcmp(psf->contextname[d], "XPEAK_IMAGE"))
        psf->cx = d;
      else if (!strcmp(psf->contextname[d], "Y_IMAGE")
		|| !strcmp(psf->contextname[d], "YWIN_IMAGE")
		|| !strcmp(psf->contextname[d], "YPSF_IMAGE")
		|| !strcmp(psf->contextname[d], "YMODEL_IMAGE")
		|| !strcmp(psf->contextname[d], "YPEAK_IMAGE"))
        psf->cy = d;
      }
    }

/* Allocate an array of Moffat function fits */
/*
  QMALLOC(psf->moffat, moffatstruct, nsnap);
  QMALLOC(psf->pfmoffat, moffatstruct, nsnap);
*/

/* Free temporary arrays */
  if (ndim)
    {
    free(names2);
    free(group2);
    free(dim2);
    }

 return psf;
  }


/****** psf_inherit ***********************************************************
PROTO	psfstruct *psf_inherit(contextstruct *context, psfstruct *psf)
PURPOSE	Initialize a PSF structure based on a preexisting PSF and a new context.
INPUT	Pointer to context structure,
	pointer to existing PSF.
OUTPUT  psfstruct pointer.
NOTES   The maximum degrees and number of dimensions allowed are set in poly.h.
AUTHOR  E. Bertin (IAP)
VERSION 03/09/2009
 ***/
psfstruct	*psf_inherit(contextstruct *context, psfstruct *psf)
  {
   psfstruct	*newpsf;
   int		c,co, ncnew,ncold, npix;

/* 10000 is just a dummy number */
  newpsf = psf_init(context, psf->size, psf->pixstep, psf->pixsize, 10000);
  newpsf->fwhm = psf->fwhm;
  npix = psf->size[0]*psf->size[1];
  if (psf->pixmask)
    QMEMCPY(psf->pixmask, newpsf->pixmask, int, npix);
  if (psf->basis)
    QMEMCPY(psf->basis, newpsf->basis, float, psf->nbasis*npix);
  ncold = psf->poly->ndim;
  ncnew = newpsf->poly->ndim;
  for (c=0; c<ncnew; c++)
    for (co=0; co<ncold; co++)
      if (!strcmp(newpsf->contextname[c],psf->contextname[co]))
        {
        newpsf->contextoffset[c] = psf->contextoffset[co];
        newpsf->contextscale[c] = psf->contextscale[co];
        }

  return newpsf;
  }


/****** psf_end ***************************************************************
PROTO   void psf_end(psfstruct *psf)
PURPOSE Free a PSF structure and everything it contains.
INPUT   psfstruct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 14/10/2009
 ***/
void	psf_end(psfstruct *psf)
  {
   int	d, ndim;

  ndim = psf->poly->ndim;
  for (d=0; d<ndim; d++)
    free(psf->contextname[d]);
  free(psf->contextname);
  free(psf->contextoffset);
  free(psf->contextscale);
  poly_end(psf->poly);
  free(psf->pixmask);
  free(psf->basis);
  free(psf->basiscoeff);
  free(psf->comp);
  free(psf->loc);
  free(psf->resi);
  free(psf->size);
  free(psf->moffat);
  free(psf->pfmoffat);
  free(psf->homo_kernel);
  free(psf);

  return;
  }


/****** psf_copy **************************************************************
PROTO   psfstruct *psf_copy(psfstruct *psf)
PURPOSE Copy a PSF structure and everything it contains.
INPUT   psfstruct pointer.
OUTPUT  psfstruct pointer.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 25/09/2011
 ***/
psfstruct *psf_copy(psfstruct *psf)
  {
   psfstruct	*newpsf;
   int		d, ndim,npix,nsnap;

  QMALLOC(newpsf, psfstruct, 1);
  *newpsf = *psf;
  ndim = psf->poly->ndim;
  QMEMCPY(psf->size, newpsf->size, int, psf->dim);
  QMALLOC(newpsf->contextname, char *, ndim);
  for (d=0; d<ndim; d++)
    {
    QMALLOC(newpsf->contextname[d], char, 80);
    strncpy(newpsf->contextname[d], psf->contextname[d], 80);
    }
  if (ndim)
    {
    QMEMCPY(psf->contextoffset, newpsf->contextoffset, double, ndim);
    QMEMCPY(psf->contextscale, newpsf->contextscale, double, ndim);
    }
  newpsf->poly = poly_copy(psf->poly);
  npix = psf->size[0]*psf->size[1];
  if (psf->pixmask)
    QMEMCPY(psf->pixmask, newpsf->pixmask, int, npix);
  if (psf->basis)
    QMEMCPY(psf->basis, newpsf->basis, float, psf->nbasis*npix);
  if (psf->basiscoeff)
    QMEMCPY(psf->basiscoeff, newpsf->basiscoeff, float,
	psf->nbasis*psf->poly->ncoeff);
  QMEMCPY(psf->comp, newpsf->comp, float, psf->npix);
  QMEMCPY(psf->loc, newpsf->loc, float, npix);
  QMEMCPY(psf->resi, newpsf->resi, float, npix);
  nsnap = 1;
  for (d=0; d<ndim; d++)
    {
    nsnap *= psf->nsnap;
    }

  if (psf->moffat)
    QMEMCPY(psf->moffat, newpsf->moffat, moffatstruct, nsnap);
  if (psf->pfmoffat)
    QMEMCPY(psf->pfmoffat, newpsf->pfmoffat, moffatstruct, nsnap);
  if (psf->homo_kernel)
    QMEMCPY(psf->homo_kernel, newpsf->homo_kernel, float, psf->npix);

  return newpsf;
  }


/****** psf_make **************************************************************
PROTO	void	psf_make(psfstruct *psf, setstruct *set, double prof_accuracy)
PURPOSE	Make the PSF.
INPUT	Pointer to the PSF,
	Pointer to the sample set,
	PSF accuracy.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 20/11/2012
 ***/
void	psf_make(psfstruct *psf, setstruct *set, double prof_accuracy)
  {
   polystruct	*poly;
   samplestruct	*sample;
   double	*pstack,*wstack, *basis, *pix,*wpix, *coeff, *pos, *post;
   float	*comp,*image,*imaget, *weight,*weightt,
		backnoise2, gain, norm, norm2, noise2, profaccu2, pixstep, val;
   int		i,c,n, ncoeff,npix,nsample;

  poly = psf->poly;

/* First copy the offset and scaling information from the set structure */
  for (i=0; i<poly->ndim; i++)
    {
    psf->contextoffset[i] = set->contextoffset[i];
    psf->contextscale[i] = set->contextscale[i];
    }

  nsample = set->nsample;
  if (!nsample)
    return;

  ncoeff = poly->ncoeff;
  npix = psf->size[0]*psf->size[1];
  QCALLOC(image, float, nsample*npix);
  QMALLOC(weight, float, nsample*npix);
  QMALLOC(pstack, double, nsample);
  QMALLOC(wstack, double, nsample);
  QMALLOC(pos, double, poly->ndim?(nsample*poly->ndim):1);
  QMALLOC(basis, double, poly->ncoeff*nsample);
  pixstep = psf->pixstep>1.0? psf->pixstep : 1.0;
  post = pos;
  for (n=0; n<nsample; n++)
    {
    sample = &set->sample[n];
/*-- Normalize approximately the image and produce a weight-map */
    norm = sample->norm;
    norm2 = norm*norm;
    profaccu2 = (float)(prof_accuracy*prof_accuracy)*norm2;
    gain = sample->gain;
    backnoise2 = sample->backnoise2;
    imaget = image+n*npix;
    weightt = weight+n*npix;
    vignet_resample(sample->vig, set->vigsize[0], set->vigsize[1],
	imaget, psf->size[0], psf->size[1],
	sample->dx, sample->dy, psf->pixstep, pixstep);
    for (i=npix; i--;)
      {
      val = (*(imaget++) /= norm);
      noise2 = backnoise2 + profaccu2*val*val;
      if (val>0.0 && gain>0.0)
        noise2 += val/gain;
      *(weightt++) = norm2/noise2;      
      }

    for (i=0; i<poly->ndim; i++)
      *(post++) = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    }

/* Make a polynomial fit to each pixel */
  for (i=0; i<npix; i++)
    {
    imaget = image+i;
    weightt = weight+i;
    pix=pstack;
    wpix=wstack;
/*-- Stack ith pixel from each PSF candidate */
    for (n=nsample; n--;)
      {
      *(pix++) = (double)*imaget;
      *(wpix++) = (double)*weightt;
      imaget += npix; 
      weightt += npix;     
      }

/*-- Polynomial fitting */
    poly_fit(poly, i?NULL:pos, pstack, wstack, nsample, basis, 1000.0);

/*-- Store as a PSF component */
    for (coeff=poly->coeff, comp=psf->comp+i,  c=ncoeff; c--; comp+=npix)
      *comp = *(coeff++);
    }

  free(image);
  free(weight);
  free(pstack);
  free(wstack);
  free(basis);
  free(pos);

  return;
  }


/****** psf_build *************************************************************
PROTO	void	psf_build(psfstruct *psf, double *pos)
PURPOSE	Build the local PSF (function of "coordinates").
INPUT	Pointer to the PSF,
	Pointer to the (context) coordinates.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 02/04/2013
 ***/
void	psf_build(psfstruct *psf, double *pos)
  {
   double	*basis;
   float	*ppc, *pl, fac;
   int		n,p, npix;

  npix = psf->size[0]*psf->size[1];
/* Reset the Local PSF mask */
  memset(psf->loc, 0, npix*sizeof(float));

  poly_func(psf->poly, pos);
  basis = psf->poly->basis;

  ppc = psf->comp;
/* Sum each component */
  for (n = (psf->dim>2?psf->size[2]:1); n--;)
    {
    pl = psf->loc;
    fac = (float)*(basis++);
#pragma ivdep
    for (p=npix; p--;)
      *(pl++) +=  fac**(ppc++);
    }

  return;
  }


/****** psf_makeresi **********************************************************
PROTO	void	psf_makeresi(psfstruct *psf, setstruct *set, int centflag,
		double prof_accuracy)
PURPOSE	Compute PSF residuals.
INPUT	Pointer to the PSF,
	Pointer to the sample set,
	Re-centering flag (0=no),
	PSF accuracy parameter.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 10/07/2012
 ***/
void	psf_makeresi(psfstruct *psf, setstruct *set, int centflag,
		double prof_accuracy)
  {
   samplestruct		*sample;
   static double	pos[MAXCONTEXT], amat[9], bmat[3];
   double		*dresi, *dresit, *amatt,
			*cvigx,*cvigxt, *cvigy,*cvigyt,
			nm1, chi2, dx,dy, ddx,ddy, dval,dvalx,dvaly,dwval,
			radmin2,radmax2, hcw,hch, yb, mx2,my2,mxy,
			xc,yc,rmax2,x,y, mse, xi2, xyi, resival, resinorm;
   float		*vigresi, *vig, *vigw, *fresi,*fresit, *vigchi,
			*cbasis,*cbasist, *cdata,*cdatat, *cvigw,*cvigwt,
			norm, fval, vigstep, psf_extraccu2, wval, sval;
   int			i,j,n,ix,iy, ndim,npix,nsample, cw,ch,ncpix, okflag,
			accuflag, nchi2;

  accuflag = (prof_accuracy > 1.0/BIG);
  vigstep = 1/psf->pixstep;
  nsample = set->nsample;
  npix = set->vigsize[0]*set->vigsize[1];
  ndim = psf->poly->ndim;
  QCALLOC(dresi, double, npix);

  if (centflag)
    {
/*-- Compute Centering sub-vignet size (containing most of the signal) */
    cw=ch=(int)(2*set->fwhm+1.0);
    if (cw>set->vigsize[0])
      cw=set->vigsize[0];
    if (ch>set->vigsize[1])
      ch=set->vigsize[1];
/*-- Allocate memory for the sub-vignet */
    ncpix = cw*ch;
    QMALLOC(cdata, float, ncpix);
    QMALLOC(cbasis, float, ncpix);
    QMALLOC(cvigw, float, ncpix);
    QMALLOC(cvigx, double, ncpix);
    QMALLOC(cvigy, double, ncpix);
/*-- Initialize gradient image */
    hcw = (double)(cw/2);
    hch = (double)(ch/2);
    cvigxt = cvigx;
    cvigyt = cvigy;
    for (iy=0; iy<ch; iy++)
      {
      yb = iy-hch;
      for (ix=0; ix<cw; ix++)
        {
        *(cvigxt++) = ix-hcw;
        *(cvigyt++) = yb;
        }
      }
    }
  else
    {
    cvigx = cvigy = (double *)NULL;	/* To avoid gcc -Wall warnings */
    cbasis = cdata = cvigw = (float *)NULL;	/* ditto */
    cw = ch = ncpix = 0;			/* ibid */
    }

/* Set convergence boundaries */
  radmin2 = PSF_MINSHIFT*PSF_MINSHIFT;
  radmax2 = PSF_MAXSHIFT*PSF_MAXSHIFT;
  okflag = nchi2 = 0;
  mse = 0.0; 				/* To avoid gcc -Wall warnings */

/* Compute the chi2 */
  for (sample=set->sample, n=nsample; n--; sample++)
    {
/*-- Build the local PSF */
    for (i=0; i<ndim; i++)
      pos[i] = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    psf_build(psf, pos);

/*-- Delta-x and Delta-y in vignet-pixel units */
    dx = sample->dx;
    dy = sample->dy;

    if (centflag)
      {
/*---- Copy the data into the sub-vignet */
      vignet_copy(sample->vig, set->vigsize[0], set->vigsize[1],
		cdata, cw,ch, 0,0, VIGNET_CPY);
/*---- Weight the data */
      vignet_copy(sample->vigweight, set->vigsize[0], set->vigsize[1],
		cvigw, cw,ch, 0,0, VIGNET_CPY);

      for (cdatat=cdata, cvigwt=cvigw, i=ncpix; i--;)
        *(cdatat++) *= *(cvigwt++);

      for (j=0; j<PSF_NITER; j++)
        {
/*------ Map the PSF model at the current position */
        vignet_resample(psf->loc, psf->size[0], psf->size[1],
		cbasis, cw,ch, -dx*vigstep, -dy*vigstep, vigstep, 1.0);

/*------ Build the a and b matrices */
        memset(amat, 0, 9*sizeof(double));
        bmat[0] = bmat[1] = bmat[2] = mx2=my2=mxy = 0.0;
        for (cvigxt=cvigx,cvigyt=cvigy,cvigwt=cvigw,
		cbasist=cbasis,cdatat=cdata, i=ncpix; i--;)
          {
          dval = (double)*(cbasist++);
          bmat[0] += (dwval = dval*(double)*(cdatat++));
          bmat[1] += dwval*(dvalx = *(cvigxt++) - dx);
          bmat[2] += dwval*(dvaly = *(cvigyt++) - dy);
          mx2 += dval*dvalx*dvalx;
          my2 += dval*dvaly*dvaly;
          mxy += dval*dvalx*dvaly;
          amatt=amat;
          *(amatt++) += (dval *= dval*(double)*(cvigwt++));
          *(amatt++) += dval*dvalx;
          *(amatt++) += dval*dvaly;
          *(++amatt) += dval*dvalx*dvalx;
          *(++amatt) += dval*dvalx*dvaly;
          *(amatt+3) += dval*dvaly*dvaly;
          }

/*------ Solve the system */
#if defined(HAVE_LAPACKE)
 #ifdef MATSTORAGE_PACKED
        if (LAPACKE_dppsv(LAPACK_COL_MAJOR,'L',3,1,amat,bmat,3) != 0)
 #else
        if (LAPACKE_dposv(LAPACK_COL_MAJOR,'L',3,1,amat,3,bmat,3) != 0)
 #endif
#else
        if (clapack_dposv(CblasRowMajor,CblasUpper,3,1,amat,3,bmat,3) != 0)
#endif
          warning("Not a positive definite matrix", " in PSF model solver");

/*------ Convert to a shift */
        dx += 0.5*(ddx = (bmat[1]*mx2 + bmat[2]*mxy) / bmat[0]); 
        dy += 0.5*(ddy = (bmat[2]*my2 + bmat[1]*mxy) / bmat[0]); 
/*------ Exit if it converges or diverges */
        if (ddx*ddx+ddy*ddy < radmin2)
          {
          okflag = 1;
          break;
	  }
        else if (dx*dx+dy*dy > radmax2)
          break;
        }
      if (okflag)
        {
        sample->dx = dx;
        sample->dy = dy;
        }
      }


/*-- Map the PSF model at the current position */
    vignet_resample(psf->loc, psf->size[0], psf->size[1],
	sample->vigresi, set->vigsize[0], set->vigsize[1],
	-dx*vigstep, -dy*vigstep, vigstep, 1.0);
/*-- Fit the flux */
    xi2 = xyi = 0.0;
    for (cvigwt=sample->vigweight,cbasist=sample->vigresi,cdatat=sample->vig,
	i=npix; i--;)
      {
      dwval = *(cvigwt++);
      dval = (double)*(cbasist++);
      xi2 += dwval*dval*dval;
      xyi += dwval*dval*(double)*(cdatat++);
      }

    norm = (xi2>0.0)? xyi/xi2 : sample->norm;

/*-- Subtract the PSF model and compute Chi2 */
    chi2 = mse = resival = resinorm = 0.0;
    dresit = dresi;
    psf_extraccu2 = prof_accuracy*prof_accuracy*norm*norm;
    xc = (double)(set->vigsize[0]/2)+sample->dx;
    yc = (double)(set->vigsize[1]/2)+sample->dy;
    y = -yc;
    rmax2 = psf->pixstep*(psf->size[0]<psf->size[1]?
		(double)(psf->size[0]/2) : (double)(psf->size[1]/2));
    rmax2 *= rmax2;
    nchi2 = 0;
    vig = sample->vig;
    vigw = sample->vigweight;
    vigresi=sample->vigresi;
    vigchi = sample->vigchi;
    for (iy=set->vigsize[1]; iy--; y+=1.0)
      {
      x = -xc;
#pragma ivdep
      for (ix=set->vigsize[0]; ix--; x+=1.0, vig++, vigresi++, dresit++,
						vigchi++)
        {
        *vigchi = 0;
        if ((wval=*(vigw++))>0.0)
          {
          if (accuflag)
            wval = 1.0/(1.0 / wval + psf_extraccu2**vigresi**vigresi);
          *vigresi = fval = (*vig-*vigresi*norm);
          if (x*x+y*y<rmax2)
            {
            mse += fval*fval;
            nchi2++;
            chi2 += (double)(*vigchi=wval*fval*fval);
            *dresit += fval;
            sval = *vig+*vigresi*norm;
            resival += sval*fabsf(fval);
            resinorm += sval*sval;
            }
          }
        }
      }

    sample->chi2 = (nchi2> 1)? chi2/(nchi2-1) : chi2;
    sample->modresi = (resinorm > 0.0)? 2.0*resival/resinorm : resival;
    }

/* Normalize and convert to floats the Residual array */
  mse = sqrt(mse/nsample/nchi2);
/*printf("%g\n", mse);*/
  QMALLOC(fresi, float, npix); 
  nm1 = nsample > 1?  (double)(nsample - 1): 1.0;
  for (dresit=dresi,fresit=fresi, i=npix; i--;)
      *(fresit++) = sqrt(*(dresit++)/nm1);

/*-- Map the residuals to PSF coordinates */
  vignet_resample(fresi, set->vigsize[0], set->vigsize[1],
	psf->resi, psf->size[0], psf->size[1], 0.0,0.0, psf->pixstep, 1.0);

/* Free memory */
  free(dresi);
  free(fresi);
  if (centflag)
    {
    free(cvigx);
    free(cvigy);
    free(cvigw);
    free(cbasis);
    free(cdata);
    }

  return;
  }


/****** psf_refine ************************************************************
PROTO	int	psf_refine(psfstruct *psf, setstruct *set)
PURPOSE	Refine PSF by solving a system to recover "aliased" components.
INPUT	Pointer to the PSF,
	Pointer to the sample set.
OUTPUT  RETURN_OK if a PSF is succesfully computed, RETURN_ERROR otherwise.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 02/04/2013
 ***/
int	psf_refine(psfstruct *psf, setstruct *set)
  {
   polystruct		*poly;
   samplestruct		*sample;
   double		pos[MAXCONTEXT];
   char			str[MAXCHAR];
   double		*desmat,*desmatt,*desmatt2, *desmat0,*desmat02,
			*bmat,*bmatt, *basis,*basist, *basist2,
			*sigvig,*sigvigt, *alphamat,*alphamatt,
			*betamat,*betamatt,*betamat2, *coeffmat,*coeffmatt,
			dx,dy, dval, norm, tikfac;
   float		*vig,*vigt,*vigt2, *wvig,
			*vecvig,*vecvigt, *ppix, *vec, *bcoeff,
			vigstep;
   int			*desindex,*desindext,*desindext2,
			*desindex0,*desindex02;
   int			i,j,jo,k,l,c,n, npix,nvpix, ndata,ncoeff,nsample,npsf,
			ncontext, nunknown, matoffset, dindex;

/* Exit if no pixel is to be "refined" or if no sample is available */
  if (!set->nsample || !psf->basis)
    return RETURN_ERROR;

  npix = psf->size[0]*psf->size[1];
  nvpix = set->vigsize[0]*set->vigsize[1];
  vigstep = 1/psf->pixstep;

  npsf = psf->nbasis;
  ndata = psf->ndata? psf->ndata : set->vigsize[0]*set->vigsize[1]+1;
  poly = psf->poly;
  ncontext = set->ncontext;
  ncoeff = poly->ncoeff;
  nsample = set->nsample;
  nunknown = ncoeff*npsf;

/* Prepare a vignet that will contain each projected basis vector */
  QCALLOC(vecvig, float, nvpix);

//  NFPRINTF(OUTPUT,"Processing samples...");
  matoffset =nunknown-ncoeff;		/* Offset between matrix coeffs */
/* Set-up the (compressed) design matrix and data vector */
  QCALLOC(desmat, double, npsf*ndata);
  QCALLOC(desindex, int, npsf*ndata);
  QMALLOC(bmat, double, nvpix);
/* ... a matrix containing the context coefficient submatrix... */
  QMALLOC(coeffmat, double, ncoeff*ncoeff);
/* ... a vignet that will contain the current vignet residuals... */
  QMALLOC(vig, float, nvpix);
/* ... a vignet that will contain the current 1/sigma map... */
  QMALLOC(sigvig, double, nvpix);
/* ... and allocate some more for storing the normal equations */
  QCALLOC(alphamat, double, nunknown*nunknown);
  QCALLOC(betamat, double, nunknown);
/*
  psf_orthopoly(psf, set);
*/
/* Go through each sample */
  for (sample=set->sample, n=0; n<nsample ; n++, sample++)
    {
    sprintf(str, "Processing sample #%d", n+1);
//    NFPRINTF(OUTPUT, str);
/*-- Delta-x and Delta-y in PSF-pixel units */
    dx = -sample->dx*vigstep;
    dy = -sample->dy*vigstep;
    norm = (double)sample->norm;

/*-- Build the local PSF */
    for (i=0; i<ncontext; i++)
      pos[i] = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    psf_build(psf, pos);

/*-- Build the current context coefficient sub-matrix */
    basis = poly_ortho(poly, poly->basis, poly->orthobasis);
    for (basist=basis, coeffmatt=coeffmat, l=ncoeff; l--;)
      for (dval=*(basist++), basist2=basis, i=ncoeff; i--;)
        *(coeffmatt++) = dval**(basist2++);

/*-- Precompute the 1/sigma-map for the current sample */
    for (sigvigt=sigvig, wvig=sample->vigweight, i=nvpix; i--;)
      *(sigvigt++) = sqrt(*(wvig++));

/*-- Go through each relevant PSF pixel */
    desmatt = desmat;
    desindext = desindex;
    if (psf->pixmask)
      {
/*---- Map the PSF model at the current position */
      vignet_resample(psf->loc, psf->size[0], psf->size[1],
		vig, set->vigsize[0], set->vigsize[1], dx, dy, vigstep, 1.0);
/*---- Subtract the PSF model */
      for (vigt=vig, vigt2=sample->vig, i=nvpix; i--; vigt++)
          *vigt = (float)(*(vigt2++) - *vigt*norm);
      }
    else
/*---- Simply copy the image data */
      for (vigt=vig, vigt2=sample->vig, i=nvpix; i--;)
        *(vigt++) = (float)*(vigt2++);
    for (i=0; i<npsf; i++)
      {
/*---- Shift the current basis vector to the current PSF position */
      vignet_resample(&psf->basis[i*npix], psf->size[0], psf->size[1],
		vecvig, set->vigsize[0],set->vigsize[1], dx,dy, vigstep, 1.0);
/*---- Retrieve coefficient for each relevant data pixel */
      for (vecvigt=vecvig, sigvigt=sigvig,
		desmatt2=desmatt, desindext2=desindext, j=jo=0; j++<nvpix;)
        if (fabs(dval = *(vecvigt++) * *(sigvigt++)) > (1/BIG))
          {
          *(desmatt2++) = norm*dval;
          *(desindext2++) = (j-jo);
          jo = j;
          }

      *desindext2 = 0;

      desindext += ndata;
      desmatt += ndata;
      }

/*-- Fill the b matrix with data points */
    for (vigt=vig, sigvigt=sigvig, bmatt=bmat, j=nvpix; j--;)
      *(bmatt++) = *(vigt++) * *(sigvigt++);

/*-- Compute the matrix of normal equations */
    betamatt = betamat;
    for (desmat0=desmat, desindex0=desindex, k=0; k<npsf;
		desmat0+=ndata, desindex0+=ndata, k++)
      {
      for (desmat02=desmat0, desindex02=desindex0, j=k; j<npsf;
		desmat02+=ndata, desindex02+=ndata, j++)
        {
        dval = 0.0;
        desmatt=desmat0;
        desmatt2=desmat02;
        desindext=desindex0;
        desindext2=desindex02;
        dindex=*desindext-*desindext2;
        while (*desindext && *desindext2)
          {
          while (*desindext && dindex<0)
            {
            dindex+=*(++desindext);
            desmatt++;
            }
          while (*desindext2 && dindex>0)
            {
            dindex-=*(++desindext2);
            desmatt2++;
            }
          while (*desindext && !dindex)
            {
            dval += *(desmatt++)**(desmatt2++);
            dindex = *(++desindext)-*(++desindext2);
            }
          }
        if (fabs(dval) > (1/BIG))
          {
          alphamatt = alphamat+(j+k*npsf*ncoeff)*ncoeff;
          for (coeffmatt=coeffmat, l=ncoeff; l--; alphamatt+=matoffset)
            for (i=ncoeff; i--;)
              *(alphamatt++) += dval**(coeffmatt++);
          }
        }
      dval = 0.0;
      desmatt=desmat0;
      desindext=desindex0;
      bmatt=bmat-1;
      while (*desindext)
        dval += *(desmatt++)**(bmatt+=*(desindext++));
      for (basist=basis,i=ncoeff; i--;)
        *(betamatt++) += dval**(basist++);
      }
    }

/* Free some memory... */
  free(coeffmat);
  free(desmat);
  free(desindex);
  free(bmat);
  free(vecvig);
  free(vig);
  free(sigvig);

/* Basic Tikhonov regularisation */
  if (psf->pixmask)
    {
    tikfac= 0.01;
    tikfac = 1.0/(tikfac*tikfac);
    for (i=0; i<nunknown; i++)
      alphamat[i+nunknown*i] += tikfac;
    }

//  NFPRINTF(OUTPUT,"Solving the system...");

#if defined(HAVE_LAPACKE)
 #ifdef MATSTORAGE_PACKED
  if (LAPACKE_dppsv(LAPACK_COL_MAJOR,'L',nunknown,1,alphamat,betamat,nunknown)
	!= 0)
 #else
  if (LAPACKE_dposv(LAPACK_COL_MAJOR,'L',nunknown,1,alphamat,nunknown,
	betamat,nunknown) != 0)
 #endif
#else
  if (clapack_dposv(CblasRowMajor,CblasUpper,nunknown,1,alphamat,nunknown,
	betamat,nunknown) != 0)
#endif
    warning("Not a positive definite matrix"," in PSF model refinement solver");

/* Check whether the result is coherent or not */
#if defined(HAVE_ISNAN2) && defined(HAVE_ISINF)
  if (isnan(*betamat) || isinf(*betamat))
#else
  if ((0x7ff00000 & *(unsigned int *)((char *)betamat+4)) == 0x7ff00000)
#endif
    {
/*-- If not, exit without doing anything */
    warning("Insufficient constraints for deriving/refining PSF", "");
/*-- Free all */
    free(alphamat);
    free(betamat);
    return RETURN_ERROR;
    }

//  NFPRINTF(OUTPUT,"Updating the PSF...");
  bcoeff = NULL;		/* To avoid gcc -Wall warnings */
  if (psf->basiscoeff)
    {
    free(psf->basiscoeff);
    psf->basiscoeff = NULL;
    }
  if (!psf->pixmask)
    {
    memset(psf->comp, 0, npix*ncoeff*sizeof(float));
    QMALLOC(psf->basiscoeff, float, nunknown);
    bcoeff = psf->basiscoeff;
    }
  QMALLOC(betamat2, double, ncoeff);
  for (j=0; j<npsf; j++)
    {
    ppix = psf->comp;
    betamatt = poly_deortho(poly, betamat + j*ncoeff, betamat2);
    for (c=ncoeff; c--;)
      {
      vec = &psf->basis[j*npix];
      dval = *(betamatt++);
#pragma ivdep
      for (i=npix; i--;)
        *(ppix++) += dval**(vec++);
      if (psf->basiscoeff)
/*---- Copy the basis coefficients (to be written later in PSF file) */
        *(bcoeff++) = (float)dval;
      }
    }

/* Free all */
  free(alphamat);
  free(betamat);
  free(betamat2);

  return RETURN_OK;
  }


/****** psf_orthopoly *********************************************************
PROTO	void	psf_orthopoly(psfstruct *psf)
PURPOSE	Orthonormalize the polynomial basis over the range of possible contexts.
INPUT	PSF structure.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 02/04/2013
 ***/
void psf_orthopoly(psfstruct *psf, setstruct *set)
  {
   samplestruct	*sample;
   polystruct	*poly;
   double	pos[POLY_MAXDIM],
		*basis, *data,*datat,
		norm;
   int		c,i,n, ndim, ncoeff, ndata;

  poly = psf->poly;
  ncoeff = poly->ncoeff;
  ndim = poly->ndim;
  ndata = set->nsample;
  norm = -1.0/sqrt(ndata);

  QMALLOC(data, double, ndata*ncoeff);
/* Go through each sample */
  for (sample=set->sample, n=0; n<ndata ; n++, sample++)
    {
/*-- Get the local context coordinates */
    for (i=0; i<ndim; i++)
      pos[i] = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    poly_func(poly, pos);
    basis = poly->basis;
    datat = data + n;
/*-- Fill basis matrix as a series of row vectors */
#pragma ivdep
    for (c=ncoeff; c--; datat+=ndata)
      *datat = *(basis++)*norm;
    }

  poly_initortho(poly, data, ndata);
  free(data);

  return;
  }


/****** psf_makebasis *********************************************************
PROTO	void	psf_makebasis(psfstruct *psf, setstruct *set,
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
void	psf_makebasis(psfstruct *psf, setstruct *set,
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
      psf->nbasis = psf_pshapelet(&psf->basis, psf->size[0],psf->size[1],
		nvec, sqrt(nvec+1.0)*prefs.basis_scale);
      break;
    case BASIS_FILE:
      psf->nbasis = psf_readbasis(psf, prefs.basis_name, 0);
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: unknown PSF vector basis in ",
			"psf_makebasis()");
    }

  return;
  }


/****** psf_laguerre **********************************************************
PROTO	double	psf_laguerre(double x, int p, int q)
PURPOSE	Return Laguerre polynomial value.
INPUT	x,
	p,
	q.
OUTPUT  Value of the Laguerre polynomial.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2007
 ***/
static double	psf_laguerre(double x, int p, int q)
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


/****** psf_pshapelet *********************************************************
PROTO	int psf_pshapelet(float **shape, int w, int h, int nmax, double beta)
PURPOSE	Compute Polar shapelet basis set.
INPUT	Pointer to the array of image vectors (which will be allocated),
	Image vector width,
	Image vector height,
	Shapelet n_max,
	beta parameter.
OUTPUT  Total number of image vectors generated.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2007
 ***/
int psf_pshapelet(float **basis, int w, int h, int nmax, double beta)
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
          val += fac*pow(*fr2t, dm/2.0)*psf_laguerre(*fr2t, hnmm, m)
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
            val += fac*pow(*fr2t, dm/2.0)*psf_laguerre(*fr2t, hnmm, m)
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


/****** psf_readbasis *********************************************************
PROTO   int psf_readbasis(psfstruct *psf, char *filename, int ext)
PURPOSE Read a set of basis functions for the PSF from a 3D FITS-file.
INPUT   Pointer to the PSF structure,
	FITS filename,
	Extension number.
OUTPUT  Number of basis vectors read.
NOTES   The maximum degrees and number of dimensions allowed are set in poly.h.
AUTHOR  E. Bertin (IAP)
VERSION 13/11/2007
 ***/
int	psf_readbasis(psfstruct *psf, char *filename, int ext)
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


