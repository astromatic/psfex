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
*	Last modify:	16/12/2008
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
#include	"fastromath.h"
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
VERSION 16/12/2008
 ***/
void	psf_diagnostic(psfstruct *psf)
  {
   moffatstruct		*moffat;
   double		dpos[POLY_MAXDIM];
   float		param[PSF_DIAGNPARAM],
			*dresi,
			dstep,dstart, fwhm, temp;
   int			i,m,n, w,h, npc,nt, nmed;

  nmed = 0;
  npc = psf->poly->ndim;
  for (nt=1, i=npc; (i--)>0;)
    {
    nt *= psf->nsnap;
    nmed += nmed*psf->nsnap + (psf->nsnap-1)/2;
    }
  psf->nmed = nmed;

  w = psf->size[0];
  h = psf->size[1];
  m = w*h;
  QCALLOC(dresi, float, m);
  moffat = psf->moffat;
  dstep = 1.0/psf->nsnap;
  dstart = (1.0-dstep)/2.0;
  memset(dpos, 0, POLY_MAXDIM*sizeof(float));
  for (i=0; i<npc; i++)
     dpos[i] = -dstart;

  psf->moffat_fwhm_min = psf->moffat_ellipticity_min = psf->moffat_beta_min
		= psf->moffat_residuals_min = psf->moffat_symresiduals_min
		= BIG;
  psf->moffat_fwhm_max = psf->moffat_ellipticity_max = psf->moffat_beta_max
		= psf->moffat_residuals_max = psf->moffat_residuals_max
		= -BIG;

/* For each snapshot of the PSF */ 
  for (n=0; n<nt; n++)
    {
    psf_build(psf, dpos);
/*-- Initialize PSF parameters */
    fwhm = psf->fwhm / psf->pixstep;
/*-- Amplitude */
    param[0] = 	1.0/(psf->fwhm*psf->fwhm);
    moffat_parammin[0] = param[0]/10.0;
    moffat_parammax[0] = param[0]*10.0;
/*-- Xcenter */
    param[1] = (w-1)/2.0;
    moffat_parammin[1] = 0.0;
    moffat_parammax[1] = w - 1.0;
/*-- Ycenter */
    param[2] = (h-1)/2.0;
    moffat_parammin[2] = 0.0;
    moffat_parammax[2] = h - 1.0;
/*-- Major axis FWHM (pixels) */
    param[3] = fwhm;
    moffat_parammin[3] = fwhm/3.0;
    moffat_parammax[3] = fwhm*3.0;
/*-- Major axis FWHM (pixels) */
    param[4] = fwhm;
    moffat_parammin[4] = fwhm/3.0;
    moffat_parammax[4] = fwhm*3.0;
/*-- Position angle (deg)  */
    param[5] = 0.0;
    moffat_parammin[5] = moffat_parammax[5] = 0.0;
/*-- Moffat beta */
    param[6] = 3.0;
    moffat_parammin[6] = PSF_BETAMIN;
    moffat_parammax[6] = 10.0;
    psf_boundtounbound(param);
    slevmar_dif(psf_diagresi, param, dresi,
	PSF_DIAGNPARAM, m, 
	PSF_DIAGMAXITER, 
	NULL, NULL, NULL, NULL, psf);
    psf_unboundtobound(param);
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
    moffat[n].symresiduals = psf_symresi(psf);
    if ((temp=sqrt(psf->moffat[n].fwhm_min*psf->moffat[n].fwhm_max))
		< psf->moffat_fwhm_min)
      psf->moffat_fwhm_min = temp;
    if (temp > psf->moffat_fwhm_max)
      psf->moffat_fwhm_max = temp;
    if ((temp=(psf->moffat[n].fwhm_max-psf->moffat[n].fwhm_min)
		/ (psf->moffat[n].fwhm_max+psf->moffat[n].fwhm_min))
		< psf->moffat_ellipticity_min)
      psf->moffat_ellipticity_min = temp;
    if (temp > psf->moffat_ellipticity_max)
      psf->moffat_ellipticity_max = temp;
    if (psf->moffat[n].beta < psf->moffat_beta_min)
      psf->moffat_beta_min = psf->moffat[n].beta;
    if (psf->moffat[n].beta > psf->moffat_beta_max)
      psf->moffat_beta_max = psf->moffat[n].beta;
    if (psf->moffat[n].residuals < psf->moffat_residuals_min)
      psf->moffat_residuals_min = psf->moffat[n].residuals;
    if (psf->moffat[n].residuals > psf->moffat_residuals_max)
      psf->moffat_residuals_max = psf->moffat[n].residuals;
    if (psf->moffat[n].symresiduals < psf->moffat_symresiduals_min)
      psf->moffat_symresiduals_min = psf->moffat[n].symresiduals;
    if (psf->moffat[n].symresiduals > psf->moffat_symresiduals_max)
      psf->moffat_symresiduals_max = psf->moffat[n].symresiduals;

    for (i=0; i<npc; i++)
      if (dpos[i]<dstart-0.01)
        {
        dpos[i] += dstep;
        break;
        }
      else
        dpos[i] = -dstart;
    }

  psf->moffat_fwhm = sqrt(psf->moffat[nmed].fwhm_min
		* psf->moffat[nmed].fwhm_max);
  psf->moffat_ellipticity =
	(psf->moffat[nmed].fwhm_max-psf->moffat[nmed].fwhm_min)
	/ (psf->moffat[nmed].fwhm_max+psf->moffat[nmed].fwhm_min);
  psf->moffat_beta = psf->moffat[nmed].beta;
  psf->moffat_residuals = psf->moffat[nmed].residuals;
  psf->moffat_symresiduals = psf->moffat[nmed].symresiduals;

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
VERSION	03/09/2009
 ***/
void	psf_diagresi(float *par, float *fvec, int m, int n, void *adata)
  {
   psfstruct	*psf;
   float	*loc, *fvect,
		dx,dy,dy2, ct,st, fac,inva2,invb2, cxx,cyy,cxy, a, beta,
		dx0,dy0, dxstep,dystep;
   int		i, x,y, xd,yd, w,h;

//printf("--%g %g %g %g %g %g %g\n", par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
  psf = (psfstruct *)adata;
  psf_unboundtobound(par);
  w = psf->size[0];
  h = psf->size[1];
  ct = cosf(par[5]*PI/180.0);
  st = sinf(par[5]*PI/180.0);
  fac = 4.0*(powf(2.0, par[6]>PSF_BETAMIN? 1.0/par[6] : 1.0/PSF_BETAMIN) - 1.0);
  inva2 = fac/(par[3]>PSF_FWHMMIN? par[3]*par[3] : PSF_FWHMMIN*PSF_FWHMMIN);
  invb2 = fac/(par[4]>PSF_FWHMMIN? par[4]*par[4] : PSF_FWHMMIN*PSF_FWHMMIN);
  cxx = inva2*ct*ct + invb2*st*st;
  cyy = inva2*st*st + invb2*ct*ct;
  cxy = 2.0*ct*st*(inva2 - invb2);
  a = par[0] / (PSF_NSUBPIX*PSF_NSUBPIX);
  beta = -par[6];
  dxstep = psf->pixsize[0]/(PSF_NSUBPIX*psf->pixstep);
  dystep = psf->pixsize[1]/(PSF_NSUBPIX*psf->pixstep);
  dy0 = -par[2] - 0.5*(PSF_NSUBPIX - 1.0)*dystep;
  loc = psf->loc;
  fvect = fvec;
  for (i=w*h; i--;)
    *(fvect++) = -*(loc++);
  for (yd=PSF_NSUBPIX; yd--; dy0+=dystep)
    {
    dx0 = -par[1] - 0.5*(PSF_NSUBPIX - 1.0)*dxstep;
    for (xd=PSF_NSUBPIX; xd--; dx0+=dxstep)
      {
      fvect = fvec;
      dy = dy0;
      for (y=h; y--; dy+=1.0)
        {
        dy2 = dy*dy;
        dx = dx0;
        for (x=w; x--; dx+=1.0)
          *(fvect++) += a*powf(1.0+cxx*dx*dx + cyy*dy2 + cxy*dx*dy, beta);
        }
      }
    }


  psf_boundtounbound(par);

  return;
  }


/****** psf_normresi *********************************************************
PROTO	double psf_normresi(double *par, psfstruct *psf)
PURPOSE	Compute a normalized estimate of residuals w.r.t. a Moffat function.
INPUT	Pointer to the vector of fitted parameters,
	pointer to the PSF structure.
OUTPUT	Normalized residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/09/2009
 ***/
double	psf_normresi(float *par, psfstruct *psf)
  {
   double	norm, resi;
   float	*loc, *fvec,*fvect,
		dx,dy,dy2, ct,st, fac,inva2,invb2, cxx,cyy,cxy, a, beta,
		dx0,dy0, dxstep,dystep, val;
   int		i, x,y, xd,yd, w,h;

  w = psf->size[0];
  h = psf->size[1];
  ct = cosf(par[5]*PI/180.0);
  st = sinf(par[5]*PI/180.0);
  QMALLOC(fvec, float, w*h*sizeof(float));
  fac = 4.0*(pow(2.0, par[6]>PSF_BETAMIN? 1.0/par[6] : 1.0/PSF_BETAMIN) - 1.0);
  inva2 = fac/(par[3]>PSF_FWHMMIN? par[3]*par[3] : PSF_FWHMMIN*PSF_FWHMMIN);
  invb2 = fac/(par[4]>PSF_FWHMMIN? par[4]*par[4] : PSF_FWHMMIN*PSF_FWHMMIN);
  cxx = inva2*ct*ct + invb2*st*st;
  cyy = inva2*st*st + invb2*ct*ct;
  cxy = 2.0*ct*st*(inva2 - invb2);
  a = par[0] / (PSF_NSUBPIX*PSF_NSUBPIX);
  beta = -par[6];
  dxstep = psf->pixsize[0]/(PSF_NSUBPIX*psf->pixstep);
  dystep = psf->pixsize[1]/(PSF_NSUBPIX*psf->pixstep);
  dy0 = -par[2] - 0.5*(PSF_NSUBPIX - 1.0)*dystep;
  fvect = fvec;
  loc = psf->loc;
  for (i=w*h; i--;)
    *(fvect++) = -*(loc++);
  for (yd=PSF_NSUBPIX; yd--; dy0+=dystep)
    {
    dx0 = -par[1] - 0.5*(PSF_NSUBPIX - 1.0)*dxstep;
    for (xd=PSF_NSUBPIX; xd--; dx0+=dxstep)
      {
      fvect = fvec;
      dy = dy0;
      for (y=h; y--; dy+=1.0)
        {
        dy2 = dy*dy;
        dx = dx0;
        for (x=w; x--; dx+=1.0)
          *(fvect++) += a*powf(1.0+cxx*dx*dx + cyy*dy2 + cxy*dx*dy, beta);
        }
      }
    }

  resi = norm = 0.0;
  fvect = fvec;
  loc = psf->loc;
  for (i=w*h; i--;)
    {
    val = *(loc++);
    resi += fabsf(*(fvect++)*val);
    norm += val*val;
    }

  free(fvec);

  return norm > 0.0? resi / norm : 1.0;
  }


/****** psf_symresi *********************************************************
PROTO	double psf_symresi(psfstruct *psf)
PURPOSE	Compute a normalized estimate of PSF assymetry.
INPUT	Pointer to the PSF structure.
OUTPUT	Normalized residuals.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/04/2007
 ***/
double	psf_symresi(psfstruct *psf)
  {
   double	resi,val,valsym,valmean,norm;
   float	*loc,*locsym;
   int		i;

  loc = psf->loc;
  locsym = loc + psf->size[0]*psf->size[1];
  resi = norm = 0.0;
  for (i=psf->size[0]*psf->size[1]; i--;)
    {
    val = (double)*(loc++);
    valsym = (double)*(--locsym);
    valmean = val + valsym;
    norm += valmean*valmean;
    resi += fabs(valmean * (val - valsym));
    }

  return norm > 0.0? 2.0 * resi / norm : 1.0;
  }


/****** psf_moffat *********************************************************
PROTO	void	psf_moffat(psfstruct *psf, moffatstruct *moffat)
PURPOSE	Generate an image of a Moffat PSF.
INPUT	Pointer to the PSF structure,
	pointer to the Moffat structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/09/2009
 ***/
void	psf_moffat(psfstruct *psf, moffatstruct *moffat)
  {
   float	*loc,
		dx,dy,dy2, ct,st, fac,inva2,invb2, cxx,cyy,cxy, a, beta,
		dx0,dy0, dxstep,dystep, xc,yc;
   int		i, x,y, xd,yd, w,h;

  w = psf->size[0];
  h = psf->size[1];
  xc = moffat->xc[0];
  yc = moffat->xc[1];
  a = moffat->amplitude*psf->pixstep*psf->pixstep / (PSF_NSUBPIX*PSF_NSUBPIX);
  beta = -moffat->beta;
  ct = cos(moffat->theta*PI/180.0);
  st = sin(moffat->theta*PI/180.0);
  fac = 4*(pow(2, -1.0/beta) - 1.0);
  inva2 = fac/(moffat->fwhm_max*moffat->fwhm_max)*psf->pixstep*psf->pixstep;
  invb2 = fac/(moffat->fwhm_min*moffat->fwhm_min)*psf->pixstep*psf->pixstep;
  cxx = inva2*ct*ct + invb2*st*st;
  cyy = inva2*st*st + invb2*ct*ct;
  cxy = 2.0*ct*st*(inva2 - invb2);
  dxstep = psf->pixsize[0]/(PSF_NSUBPIX*psf->pixstep);
  dystep = psf->pixsize[1]/(PSF_NSUBPIX*psf->pixstep);
  dy0 = -yc - 0.5*(PSF_NSUBPIX - 1.0)*dystep;
  loc = psf->loc;
  for (i=w*h; i--;)
    *(loc++) = 0.0;
  for (yd=PSF_NSUBPIX; yd--; dy0+=dystep)
    {
    dx0 = -xc - 0.5*(PSF_NSUBPIX - 1.0)*dxstep;
    for (xd=PSF_NSUBPIX; xd--; dx0+=dxstep)
      {
      loc = psf->loc;
      dy = dy0;
      for (y=h; y--; dy+=1.0)
        {
        dy2 = dy*dy;
        dx = dx0;
        for (x=w; x--; dx+=1.0)
          *(loc++) += a*pow(1.0+cxx*dx*dx + cyy*dy2 + cxy*dx*dy, beta);
        }
      }
    }

  return;
  }


/****** psf_boundtounbound **************************************************
PROTO	void psf_boundtounbound(profitstruct *profit)
PURPOSE	Convert parameters from bounded to unbounded space.
INPUT	Pointer to the vector of parameters.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	17/07/2007
 ***/
void    psf_boundtounbound(float *param)
  {
   double       num,den;
   int          p;

  for (p=0; p<PSF_DIAGNPARAM; p++)
    if (moffat_parammin[p]!=moffat_parammax[p])
      {
      num = param[p] - moffat_parammin[p];
      den = moffat_parammax[p] - param[p];
      param[p] = num>1e-100? (den>1e-100? log(num/den): 200.0) : -200.0;
      }

  return;

  }


/****** psf_unboundtobound **************************************************
PROTO	void profit_unboundtobound(profitstruct *profit)
PURPOSE	Convert parameters from unbounded to bounded space.
INPUT	Pointer to the vector of parameters.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	17/07/2007
 ***/
void    psf_unboundtobound(float *param)
  {
   int          p;

  for (p=0; p<PSF_DIAGNPARAM; p++)
    if (moffat_parammin[p]!=moffat_parammax[p])
      param[p] = (moffat_parammax[p] - moffat_parammin[p])
                / (1.0 + exp(-(param[p]>200.0? 200.0 : param[p])))
                + moffat_parammin[p];

  return;
  }

