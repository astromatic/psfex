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
*	Last modify:	25/02/2007
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

int	diag_w, diag_h;

/****** psf_diagnostic *******************************************************
PROTO	void	psf_diagnostic(psfstruct *psf, out_data_struct *out)
PURPOSE	Free a PSF structure and everything it contains.
INPUT	Pointer to the PSF structure,
	pointer to the output data structure.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 25/02/2007
 ***/
void	psf_diagnostic(psfstruct *psf, out_data_struct *out)
  {
   lm_control_type	control;
   double		dpos[POLY_MAXDIM],
			*param, *dloc, *dresi, *fjac,
			diag[PSF_DIAGNPARAM], qtf[PSF_DIAGNPARAM],
			wa1[PSF_DIAGNPARAM],wa2[PSF_DIAGNPARAM],
			wa3[PSF_DIAGNPARAM],wa4[PSF_DIAGNPARAM],
			dstep,dstart;
   int			ipvt[PSF_DIAGNPARAM],
			i,m,n, w,h, npc,nt;

/* Initialize fitting */
  control.ftol =      1.0e-7;
  control.xtol =      1.0e-7;
  control.gtol =      1.0e-7;
  control.maxcall =   PSF_DIAGMAXITER;
  control.epsilon =   1.0e-2;
  control.stepbound = 100.0;

  npc = psf->poly->ndim;
  for (nt=PSF_NSNAP*PSF_NSNAP, i=npc-2; (i--)>0;)
    nt *= PSF_NSNAP;
  diag_w = w = psf->size[0];
  diag_h = h = psf->size[1];
  m = w*h;
  QMALLOC(param, double, PSF_DIAGNPARAM);
  QMALLOC(fjac, double, PSF_DIAGNPARAM*m);
  QMALLOC(dloc, double, m);
  QMALLOC(dresi, double, m);
  dstep = 1.0/PSF_NSNAP;
  dstart = (1.0-dstep)/2.0;
  memset(dpos, 0, POLY_MAXDIM*sizeof(double));
  for (i=0; i<npc; i++)
     dpos[i] = -dstart;
/* For each snapshot of the PSF */ 
  for (n=0; n<nt; n++)
    {
    psf_build(psf, dpos);
    for (i=0; i<m; i++)
      dloc[i] = psf->loc[i];
    control.info = 0;
    control.nfev = 0;
    for (i=0; i<npc; i++)
      if (dpos[i]<dstart-0.01)
        {
        dpos[i] += dstep;
        break;
        }
      else
        dpos[i] = -dstart;
    }

/* This goes through the modified legacy interface */
/*
  lm_lmdif(m, PSF_DIAGNPARAM, param, dresi,
        control.ftol, control.xtol, control.gtol,
        control.maxcall*(PSF_DIAGNPARAM+1), control.epsilon, diag, 1,
        control.stepbound, &(control.info),
        &(control.nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
        psf_diagresi, psf_diagprintout, dloc);
*/
  free(param);
  free(fjac);
  free(dloc);
  free(dresi);

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
NOTES	Input arguments are there only for compatibility purposes (unused)
AUTHOR	E. Bertin (IAP)
VERSION	23/02/2007
 ***/
void	psf_diagprintout(int n_par, double* par, int m_dat, double* fvec,
		void *data, int iflag, int iter, int nfev)
  {
/*
printf("%d:     %g      %g      %g      %g      %g      %g      %g\n",
        iter, par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
*/
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
	pointer to the data structure (unused),
	pointer to the info structure (unused).
OUTPUT	-.
NOTES	Input arguments are there only for compatibility purposes (unused)
AUTHOR	E. Bertin (IAP)
VERSION	25/02/2007
 ***/
void	psf_diagresi(double *par, int m_dat, double *fvec, void *data,
		int *info)
  {
   double	*dloc,
		dx,dy,dy2,r2;
   int		x,y;

  dloc = (double *)data;
  for (y=0; y<diag_h; y++)
    {
    dy = y - par[2];
    dy2 = dy*dy;
    for (x=0; x<diag_w; x++)
      {
      dx = x-par[1];
      r2 = par[3]*dx*dx+par[4]*dy2+par[5]*dx*dy;
      *(fvec++) = *(dloc++)- par[0]*pow(1.0 + r2, -par[6]);
      }
    }

  return;
  }


