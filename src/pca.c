/*
*				pca.c
*
* Principal Component Analysis.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		06/01/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"pca.h"
#include	"prefs.h"
#include	"psf.h"


/****** pca_onsnaps ***********************************************************
PROTO	float *pca_onsnaps(psfstruct **psf, int ncat, int npc)
PURPOSE	Make a Principal Component Analysis in pixel space of PSF models.
INPUT	Pointer to an array of PSF structures,
	Number of catalogues (PSFs),
	Number of principal components.
OUTPUT  Pointer to an array of principal component vectors.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/02/2009
 ***/
float *pca_onsnaps(psfstruct **psfs, int ncat, int npc)
  {
   psfstruct	*psf;
   char		str[MAXCHAR];
   double	dpos[POLY_MAXDIM],
		*covmat,*covmatt,
		dstep,dstart, dval;
   float	*basis, *pix,*pix1,*pix2;
   int		c,d,i,j,p,n,w,h, ndim,npix,nt;

/* Build models of the PSF over a range of dependency parameters */
  ndim = psfs[0]->poly->ndim;
  w = psfs[0]->size[0];
  h = psfs[0]->size[1];
  npix = w*h;
  nt = 1;
  for (d=0; d<ndim; d++)
    nt *= PCA_NSNAP;
  dstep = 1.0/PCA_NSNAP;
  dstart = (1.0-dstep)/2.0;

//  NFPRINTF(OUTPUT, "Setting-up the PCA covariance matrix");
  QCALLOC(covmat, double, npix*npix);
  for (c=0; c<ncat; c++)
    {
    psf = psfs[c];
    for (d=0; d<ndim; d++)
      dpos[d] = -dstart;
    for (n=0; n<nt; n++)
      {
      psf_build(psf, dpos);
      pix1 = pix = psf->loc;
/*---- Set-up the covariance/correlation matrix */
      covmatt = covmat;
      sprintf(str, "Setting-up the PCA covariance matrix (%.0f%%)...",
		100.0*((float)n/nt+c)/ncat);
//      NFPRINTF(OUTPUT, str);
      for (j=npix; j--;)
        {
        pix2 = pix;
        dval = (double)*(pix1++);
        for (i=npix; i--;)
          *(covmatt++) += dval**(pix2++);
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
    }
/* Do recursive PCA */
  QMALLOC(basis, float, npc*npix);
  for (p=0; p<npc; p++)
    {
    sprintf(str, "Computing Principal Component vector #%d...", p);
//    NFPRINTF(OUTPUT, str);
    pca_findpc(covmat, &basis[p*npix], npix);
    }

  free(covmat);

  return basis;
  }


/****** pca_oncomps ***********************************************************
PROTO	double *pca_oncomps(psfstruct **psfs, int next, int ncat, int npc)
PURPOSE	Make a Principal Component Analysis in image space on PSF model
	components.
INPUT	Pointer to an array of PSF structures,
	Number of extensions,
	Number of catalogues (PSFs),
	Number of principal components.
OUTPUT  Pointer to an array of principal component vectors.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/02/2009
 ***/
double *pca_oncomps(psfstruct **psfs, int next, int ncat, int npc)
  {
   psfstruct	*psf;
   char		str[MAXCHAR];
   double	dpos[POLY_MAXDIM],
		*comp,*comp1,*comp2,*compt,
		*covmat, *covmatt, *pc,
		dval, dstep,dstart;
   float	*pix, *vector;
   int		e, d,i,c,c1,c2,n,p, ndim, nt, npix, npixt;

/* Build models of the PSF over a range of dependency parameters */
  ndim = psfs[0]->poly->ndim;
  npix = psfs[0]->size[0]*psfs[0]->size[1];
  nt = 1;
  for (d=0; d<ndim; d++)
    nt *= PCA_NSNAP;
  npixt = npix*nt*next;
  dstep = 1.0/PCA_NSNAP;
  dstart = (1.0-dstep)/2.0;

//  NFPRINTF(OUTPUT, "Setting-up the PCA covariance matrix");
/* Compute PSF snapshots */
  QMALLOC(comp, double, ncat*npixt);
  compt = comp;
  for (c=0; c<ncat; c++)
    {
    sprintf(str, "Setting-up the PCA covariance matrix (%.0f%%)...",
		50.0*(float)c/ncat);
//      NFPRINTF(OUTPUT, str);
    for (e=0; e<next; e++)
      {
      psf = psfs[c*next+e];
      for (d=0; d<ndim; d++)
        dpos[d] = -dstart;
      for (n=nt; n--;)
        {
        psf_build(psf, dpos);
        pix = psf->loc;
        for (i=npix; i--;)
          *(compt++) = (double)*(pix++);
        for (d=0; d<ndim; d++)
          if (dpos[d]<dstart-0.01)
            {
            dpos[d] += dstep;
            break;
            }
          else
            dpos[d] = -dstart;
        }
      }
    }

/* Center data cloud on origin */
  for (i=0; i<npixt; i++)
    {
    dval = 0;
    compt = comp + i;
    for (c=ncat; c--; compt+=npixt)
      dval += *compt;
    dval /= (double)ncat;
    compt = comp + i;
    for (c=ncat; c--; compt+=npixt)
      *compt -= dval;
    }

/* Set-up the covariance/correlation matrix */
  QCALLOC(covmat, double, ncat*ncat);
  covmatt = covmat;
  for (c1=0; c1<ncat; c1++)
    {
    sprintf(str, "Setting-up the PCA covariance matrix (%.0f%%)...",
		50.0+50.0*((float)c1)/ncat);
//    NFPRINTF(OUTPUT, str);
    for (c2=0; c2<ncat; c2++)
      {
      comp1 = comp + c1*npixt;
      comp2 = comp + c2*npixt;
      dval = 0.0;
      for (i=npixt; i--;)
        dval += *(comp1++)**(comp2++);
      *(covmatt++) = dval;
      }
    }

  free(comp);

/* Do recursive PCA */
  QMALLOC(vector, float, ncat);
  QMALLOC(pc, double, ncat*npc);
  for (p=0; p<npc; p++)
    {
    sprintf(str, "Computing Principal Component vector #%d...", p);
//    NFPRINTF(OUTPUT, str);
    pca_findpc(covmat, vector, ncat);
    for (c=0; c<ncat; c++)
      pc[c*npc+p] = vector[c];
    }

  free(covmat);
  free(vector);

  return pc;
  }


/****** pca_findpc ************************************************************
PROTO	double pca_findpc(double *covmat, float *vec, int nmat)
PURPOSE	Find the principal component (the one with the highest eigenvalue) and
	subtract its contribution from the covariance matrix, using the
	iterative "power" method.
INPUT	Covariance matrix,
	output vector,
	Number of principal components.
OUTPUT  Eigenvalue (variance) of the PC.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 23/07/2008
 ***/
double pca_findpc(double *covmat, float *vec, int nmat)
  {
   double	*tmat,*xmat, *c,*t,*x,
		dtval,dtnorm, xval,xnorm, lambda;
   float	*v;
   int		i,j,n;

  QMALLOC(xmat, double, nmat);

/* First initialize the eigenvector to an arbitrary direction */
  QMALLOC(tmat, double, nmat);
  for (t=tmat,i=nmat; i--;)
    *(t++) = 1.0;

  dtnorm = 1.0;
  for (n=PCA_NITER; n-- && dtnorm>PCA_CONVEPS;)    
    {
/*-- Compute |x> = C|t> */
    xnorm = 0.0;
    for (c=covmat,x=xmat,j=nmat; j--;)
      {
      for (xval=0.0,t=tmat,i=nmat; i--;)
        xval += *(c++)**(t++);
      xnorm += xval*xval;
      *(x++) = xval;
      }
/*-- Compute |t> = |x>/||x|| and ||Delta t|| (for testing convergence) */
    xnorm = 1.0/sqrt(xnorm);
    dtnorm = 0.0;
    for (t=tmat,x=xmat,i=nmat; i--;)
      {
      dtval = *t;
      dtval -= (*(t++) = *(x++)*xnorm);
      dtnorm += dtval*dtval;
      }
    dtnorm = sqrt(dtnorm);
    }

  free(xmat);

/* Compute the eigenvalue lambda = <t|C|t> */
  lambda = 0.0;
  for (c=covmat,x=tmat,j=nmat; j--;)
    {
    for (xval=0.0,t=tmat,i=nmat; i--;)
      xval += *(c++)**(t++);
    lambda += xval**(x++);
    }

/* Finally subtract the contribution from the found PC: C -= lambda.|t><t| */
  for (c=covmat,x=tmat,j=nmat; j--;)
    for (xval=*(x++)*lambda,t=tmat,i=nmat; i--;)
      *(c++) -= xval**(t++);

/* Convert output vector to simple precision */
  for (v=vec,t=tmat,i=nmat; i--;)
    *(v++) = (float)*(t++);

  free(tmat);

  return lambda;
  }

