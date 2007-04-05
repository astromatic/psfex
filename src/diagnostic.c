  /*
 				diagnostic.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	PSF diagnostics.
*
*	Last modify:	04/04/2007
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
#include	"prefs.h"
#include	"poly.h"
#include	"psf.h"

/****** psf_diagnostic *******************************************************
PROTO	void	psf_diagnostic(psfstruct *psf)
PURPOSE	Free a PSF structure and everything it contains.
INPUT	Pointer to the PSF structure.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 05/04/2007
 ***/
void	psf_diagnostic(psfstruct *psf)
  {
   moffatstruct		*moffat;
   double		dpos[POLY_MAXDIM],
			param[PSF_DIAGNPARAM],
			parammin[PSF_DIAGNPARAM],parammax[PSF_DIAGNPARAM],
			*dresi,
			dstep,dstart, fwhm;
   int			i,m,n, w,h, npc,nt;

  npc = psf->poly->ndim;
  for (nt=1, i=npc; (i--)>0;)
    nt *= prefs.context_nsnap;
  w = psf->size[0];
  h = psf->size[1];
  m = w*h;
  QCALLOC(dresi, double, m);
  moffat = psf->moffat;
  dstep = 1.0/prefs.context_nsnap;
  dstart = (1.0-dstep)/2.0;
  memset(dpos, 0, POLY_MAXDIM*sizeof(double));
  for (i=0; i<npc; i++)
     dpos[i] = -dstart;
/* For each snapshot of the PSF */ 
  for (n=0; n<nt; n++)
    {
    psf_build(psf, dpos);
/*-- Initialize PSF parameters */
    fwhm = psf->fwhm / psf->pixstep;
/*-- Amplitude */
    param[0] = 	1.0/(fwhm*fwhm);
    parammin[0] = 0.0;
    parammax[0] = 1000.0;
/*-- Xcenter */
    param[1] = (w-1)/2.0;
    parammin[1] = 0.0;
    parammax[1] = w - 1.0;
/*-- Ycenter */
    param[2] = (h-1)/2.0;
    parammin[2] = 0.0;
    parammax[2] = h - 1.0;
/*-- Major axis FWHM (pixels) */
    param[3] = fwhm;
    parammin[3] = PSF_FWHMMIN;
    parammax[3] = PSF_FWHMMAX;
/*-- Major axis FWHM (pixels) */
    param[4] = psf->fwhm;
    parammin[4] = PSF_FWHMMIN;
    parammax[4] = PSF_FWHMMAX;
/*-- Position angle (deg)  */
    param[5] = 0.0;
    parammin[5] = -360.0;
    parammax[5] = 360.0;
/*-- Moffat beta */
    param[6] = 1.0;
    parammin[6] = 0.01;
    parammax[6] = 10.0;
    dlevmar_bc_dif(psf_diagresi, param, dresi,
	PSF_DIAGNPARAM, m, 
	parammin, parammax,
	PSF_DIAGMAXITER, 
	NULL, NULL, NULL, NULL, psf);
    moffat[n].amplitude = param[0]/(psf->pixstep*psf->pixstep);
    moffat[n].xc[0] = param[1];
    moffat[n].xc[1] = param[2];
    if (param[3] > param[4])
      {
      moffat[n].fwhm_max = param[3]*psf->pixstep;
      moffat[n].fwhm_min = param[4]*psf->pixstep;
      moffat[n].theta = (fmod(param[5]+360.0, 180.0));
      }
    else
      {
      moffat[n].fwhm_max = param[4]*psf->pixstep;
      moffat[n].fwhm_min = param[3]*psf->pixstep;
      moffat[n].theta = (fmod(param[5]+450.0, 180.0));
      }
    if (moffat[n].theta > 90.0)
      moffat[n].theta -= 180.0;
    moffat[n].beta = param[6];
    for (i=0; i<npc; i++)
      moffat[n].context[i] = dpos[i]*psf->contextscale[i]+psf->contextoffset[i];
    moffat[n].residuals = psf_normresi(param, psf);
    for (i=0; i<npc; i++)
      if (dpos[i]<dstart-0.01)
        {
        dpos[i] += dstep;
        break;
        }
      else
        dpos[i] = -dstart;
    }

  free(dresi);

  return;
  }


/****** psf_diagresi *********************************************************
PROTO	void psf_diagresi(double *par, double *fvec, int m, int n, void *adata)
PURPOSE	Provide a function returning residuals to lmfit.
INPUT	Pointer to the vector of parameters,
	pointer to the vector of residuals (output),
	number of parameters,
	number of data points,
	pointer to the PSF structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/04/2007
 ***/
void	psf_diagresi(double *par, double *fvec, int m, int n, void *adata)
  {
   psfstruct	*psf;
   double	dx,dy,dy2, ct,st, r2,fac,inva2,invb2, cxx,cyy,cxy;
   float	*loc;
   int		x,y,w,h;

//printf("--%g %g %g %g %g %g %g\n", par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
  psf = (psfstruct *)adata;
  loc = psf->loc;
  w = psf->size[0];
  h = psf->size[1];
  ct = cos(par[5]*PI/180.0);
  st = sin(par[5]*PI/180.0);
  fac = 4.0*(pow(2.0, par[6]>PSF_BETAMIN? 1.0/par[6] : 1.0/PSF_BETAMIN) - 1.0);
  inva2 = fac/(par[3]>PSF_FWHMMIN? par[3]*par[3] : PSF_FWHMMIN*PSF_FWHMMIN);
  invb2 = fac/(par[4]>PSF_FWHMMIN? par[4]*par[4] : PSF_FWHMMIN*PSF_FWHMMIN);
  cxx = inva2*ct*ct + invb2*st*st;
  cyy = inva2*st*st + invb2*ct*ct;
  cxy = 2.0*ct*st*(inva2 - invb2);
  for (y=0; y<h; y++)
    {
    dy = y - par[2];
    dy2 = dy*dy;
    for (x=0; x<w; x++)
      {
      dx = x-par[1];
      r2 = cxx*dx*dx + cyy*dy2 + cxy*dx*dy;
if (isnan(r2))
printf("%g ", r2);
      *(fvec++) = *(loc++) - par[0]*pow(1.0 + r2, -par[6]);
      }
    }
//printf("\n");
  return;
  }


/****** psf_normresi *********************************************************
PROTO	double psf_normresi(double *par, psfstruct *psf)
PURPOSE	Compute a normalized estimate of residuals w/to a Moffat function.
INPUT	Pointer to the vector of fitted parameters,
	pointer to the PSF structure.
OUTPUT	Normalized residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/04/2007
 ***/
double	psf_normresi(double *par, psfstruct *psf)
  {
   double	dx,dy,dy2, ct,st, r2,fac,inva2,invb2, cxx,cyy,cxy,
		resi,val,norm;
   float	*loc;
   int		x,y,w,h;

  resi = norm = 0.0;
  loc = psf->loc;
  w = psf->size[0];
  h = psf->size[1];
  ct = cos(par[5]*PI/180.0);
  st = sin(par[5]*PI/180.0);
  fac = 4.0*(pow(2.0, par[6]>PSF_BETAMIN? 1.0/par[6] : 1.0/PSF_BETAMIN) - 1.0);
  inva2 = fac/(par[3]>PSF_FWHMMIN? par[3]*par[3] : PSF_FWHMMIN*PSF_FWHMMIN);
  invb2 = fac/(par[4]>PSF_FWHMMIN? par[4]*par[4] : PSF_FWHMMIN*PSF_FWHMMIN);
  cxx = inva2*ct*ct + invb2*st*st;
  cyy = inva2*st*st + invb2*ct*ct;
  cxy = 2.0*ct*st*(inva2 - invb2);
  for (y=0; y<h; y++)
    {
    dy = y - par[2];
    dy2 = dy*dy;
    for (x=0; x<w; x++)
      {
      dx = x-par[1];
      r2 = cxx*dx*dx + cyy*dy2 + cxy*dx*dy;
      val = (double)*(loc++);
      norm += val*val;
      resi += fabs(val*(val - par[0]*pow(1.0 + r2, -par[6])));
      }
    }

  return norm > 0.0? resi / norm : 1.0;
  }


/****** psf_moffat *********************************************************
PROTO	void	psf_moffat(psfstruct *psf, moffatstruct *moffat)
PURPOSE	Generate an image of a Moffat PSF.
INPUT	Pointer to the PSF structure,
	pointer to the Moffat structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/03/2007
 ***/
void	psf_moffat(psfstruct *psf, moffatstruct *moffat)
  {
   double	dx,dy,dy2,r2, xc,yc, ct,st, inva2,invb2, cxx,cyy,cxy,
		amp, beta, fac;
   float	*loc;
   int		x,y,w,h;

  xc = moffat->xc[0];
  yc = moffat->xc[1];
  amp = moffat->amplitude*psf->pixstep*psf->pixstep;
  beta = moffat->beta;
  ct = cos(moffat->theta*PI/180.0);
  st = sin(moffat->theta*PI/180.0);
  fac = 4*(pow(2, 1.0/beta) - 1.0);
  inva2 = fac/(moffat->fwhm_max*moffat->fwhm_max)*psf->pixstep*psf->pixstep;
  invb2 = fac/(moffat->fwhm_min*moffat->fwhm_min)*psf->pixstep*psf->pixstep;
  cxx = inva2*ct*ct + invb2*st*st;
  cyy = inva2*st*st + invb2*ct*ct;
  cxy = 2.0*ct*st*(inva2 - invb2);
  loc = psf->loc;
  w = psf->size[0];
  h = psf->size[1];
  for (y=0; y<h; y++)
    {
    dy = y - yc;
    dy2 = dy*dy;
    for (x=0; x<w; x++)
      {
      dx = x-xc;
      r2 = cxx*dx*dx + cyy*dy2 + cxy*dx*dy;
      *(loc++) = amp*pow(1.0 + r2, -beta);
      }
    }

  return;
  }


