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
*	Last modify:	02/03/2007
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
#include	"lmfit/lmmin.h"
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
VERSION 02/03/2007
 ***/
void	psf_diagnostic(psfstruct *psf)
  {
   lm_control_type	control;
   moffatstruct		*moffat;
   double		dpos[POLY_MAXDIM],
			diag[PSF_DIAGNPARAM], qtf[PSF_DIAGNPARAM],
			wa1[PSF_DIAGNPARAM],wa2[PSF_DIAGNPARAM],
			wa3[PSF_DIAGNPARAM],
			*param, *dresi, *fjac, *wa4,
			dstep,dstart,
			sigma, fwhmmin,fwhmmax,twot,c2t,s2t,moffac;
   int			ipvt[PSF_DIAGNPARAM],
			i,m,n, w,h, npc,nt;

/* Initialize fitting */
  control.ftol =      1.0e-10;
  control.xtol =      1.0e-10;
  control.gtol =      1.0e-10;
  control.maxcall =   PSF_DIAGMAXITER;
  control.epsilon =   1.0e-8;
  control.stepbound = 100.0;

  npc = psf->poly->ndim;
  for (nt=prefs.context_nsnap*prefs.context_nsnap, i=npc-2; (i--)>0;)
    nt *= prefs.context_nsnap;
  w = psf->size[0];
  h = psf->size[1];
  m = w*h;
  QMALLOC(param, double, PSF_DIAGNPARAM);
  QMALLOC(fjac, double, PSF_DIAGNPARAM*m);
  QMALLOC(wa4, double, m);
  QMALLOC(dresi, double, m);
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
    control.info = 0;
    control.nfev = 0;
/*-- Initialize PSF parameters */
    param[1] = (w-1)/2.0;			/* x_center */
    param[2] = (h-1)/2.0;			/* y_center */
    sigma = psf->fwhm / (2.35*psf->pixstep);
    param[0] = 	1.0/(2*PI*sigma*sigma);		/* amplitude */
    param[3] = param[4] = 1.0/sigma;		/* Cxx and Cyy */
    param[5] = 0.0;				/* Cxy */
    param[6] = 1.0;				/* Moffat beta */
/*-- This goes through the modified legacy interface */
    lm_lmdif(m, PSF_DIAGNPARAM, param, dresi,
	control.ftol, control.xtol, control.gtol,
	control.maxcall*(PSF_DIAGNPARAM+1), control.epsilon, diag, 1,
	control.stepbound, &(control.info),
	&(control.nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
	psf_diagresi, psf_diagprintout, psf);
    param[3] *= param[3];
    param[4] *= param[4];
    twot = atan2(param[5], param[3]-param[4]);
    moffac = 2.0*sqrt(pow(2.0, 1.0/param[6]) - 1.0)*psf->pixstep;
    if (fabs(c2t=cos(twot)) > 0.0)
      {
      fwhmmax = moffac*sqrt(2.0/(param[3]+param[4] + (param[3]-param[4])/c2t));
      fwhmmin = moffac*sqrt(2.0/(param[3]+param[4] + (param[4]-param[3])/c2t));
      }
    else if (fabs(s2t=sin(twot)) > 0.0)
      {
      fwhmmax = moffac*sqrt(2.0/(param[3]+param[4] + param[5]/s2t));
      fwhmmin = moffac*sqrt(2.0/(param[3]+param[4] - param[5]/c2t));
      }
    else
      fwhmmin = fwhmmax = moffac*sqrt(2.0/(param[3]+param[4]));
    for (i=0; i<npc; i++)
      moffat[n].context[i] = dpos[i]*psf->contextscale[i]+psf->contextoffset[i];
    moffat[n].amplitude = param[0];
    moffat[n].xc[0] = param[1];
    moffat[n].xc[1] = param[2];
    moffat[n].fwhm_min = fwhmmin < fwhmmax? fwhmmin : fwhmmax;
    moffat[n].fwhm_max = fwhmmax > fwhmmin? fwhmmax : fwhmmin;
    moffat[n].theta = 90.0/PI*twot;
    moffat[n].beta = param[6];
    moffat[n].residuals = psf_normresi(param, psf);
/*
printf("%g %g  %g %g\n", moffat[n].context[0], moffat[n].context[1], moffat[n].beta, moffat[n].residuals);
*/
    for (i=0; i<npc; i++)
      if (dpos[i]<dstart-0.01)
        {
        dpos[i] += dstep;
        break;
        }
      else
        dpos[i] = -dstart;
    }

  free(param);
  free(fjac);
  free(dresi);
  free(wa4);

  return;
  }


/****** psf_diagprintout *****************************************************
PROTO	void psf_diagprintout(int n_par, double* par, int m_dat, double* fvec,
		void *data, int iflag, int iter, int nfev)
PURPOSE	Provide a function to print out results to lmfit.
INPUT	Number of fitted parameters,
	pointer to the vector of parameters,
	number of data points,
	pointer to the vector of residuals (output),
	pointer to the data structure (unused),
	0 (init) 1 (outer loop) 2(inner loop) -1(terminated),
	outer loop counter,
	number of calls to evaluate().
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/02/2007
 ***/
void	psf_diagprintout(int n_par, double* par, int m_dat, double* fvec,
		void *data, int iflag, int iter, int nfev)
  {
  if (0)
    lm_print_default(n_par, par, m_dat, fvec, data, iflag, iter, nfev);

  return;
  }


/****** psf_diagresi *********************************************************
PROTO	void psf_diagresi(double *par, int m_dat, double *fvec, void *data,
		int *info)
PURPOSE	Provide a function returning residuals to lmfit.
INPUT	Pointer to the vector of parameters,
	number of data points,
	pointer to the vector of residuals (output),
	pointer to the PSF structure,
	pointer to the info structure (unused).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/02/2007
 ***/
void	psf_diagresi(double *par, int m_dat, double *fvec, void *data,
		int *info)
  {
   psfstruct	*psf;
   double	dx,dy,dy2,r2;
   float	*loc;
   int		x,y,w,h;

  psf = (psfstruct *)data;
  loc = psf->loc;
  w = psf->size[0];
  h = psf->size[1];
  for (y=0; y<h; y++)
    {
    dy = y - par[2];
    dy2 = dy*dy;
    for (x=0; x<w; x++)
      {
      dx = x-par[1];
      r2 = par[3]*par[3]*dx*dx+par[4]*par[4]*dy2+par[5]*dx*dy;
      *(fvec++) = *(loc++) - par[0]*pow(1.0 + r2, -par[6]);
      }
    }

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
VERSION	26/02/2007
 ***/
double	psf_normresi(double *par, psfstruct *psf)
  {
   double	dx,dy,dy2,r2, resi, val, norm;
   float	*loc;
   int		x,y,w,h;

  resi = norm = 0.0;
  loc = psf->loc;
  w = psf->size[0];
  h = psf->size[1];
  for (y=0; y<h; y++)
    {
    dy = y - par[2];
    dy2 = dy*dy;
    for (x=0; x<w; x++)
      {
      dx = x-par[1];
      r2 = par[3]*par[3]*dx*dx+par[4]*par[4]*dy2+par[5]*dx*dy;
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
VERSION	26/02/2007
 ***/
void	psf_moffat(psfstruct *psf, moffatstruct *moffat)
  {
   double	dx,dy,dy2,r2, xc,yc, ct,st, inva2,invb2, step2, cxx,cyy,cxy,
		amp, beta;
   float	*loc;
   int		x,y,w,h;

  xc = moffat->xc[0];
  yc = moffat->xc[1];
  step2 = psf->pixstep*psf->pixstep;
  amp = moffat->amplitude*psf->pixstep*psf->pixstep;
  beta = moffat->beta;
  ct = cos(moffat->theta*PI/180.0);
  st = sin(moffat->theta*PI/180.0);
  inva2 = step2/(moffat->fwhm_max*moffat->fwhm_max);
  invb2 = step2/(moffat->fwhm_min*moffat->fwhm_min);
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


