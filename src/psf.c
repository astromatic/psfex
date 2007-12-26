  /*
 				psf.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Stuff related to building the PSF.
*
*	Last modify:	26/12/2007
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
#include	"prefs.h"
#include	"sample.h"
#include	"poly.h"
#include	"psf.h"
#include	"vignet.h"
#include	ATLAS_LAPACK_H

static double	psf_laguerre(double x, int p, int q);

/****** psf_clean *************************************************************
PROTO	double	psf_clean(psfstruct *psf, setstruct *set)
PURPOSE	Filter out PSF candidates
INPUT	Pointer to the PSF,
	Pointer to the sample set.
OUTPUT	Reduced chi2.
NOTES	-
AUTHOR	E. Bertin (IAP)
VERSION	12/11/2007
 ***/
double	psf_clean(psfstruct *psf, setstruct *set)
  {
#define	EPS	(1e-4)  /* a small number */
   samplestruct	*sample;
   double	chi2,chimean,chivar,chisig,chisig1,chival, locut,hicut;
   float	*chi, *chit,*chit2,
		chimed, chi2max;
   int		i, n, nsample;

/* First compute residuals for each sample (chi^2) */
  NFPRINTF(OUTPUT,"Computing residuals...");
  psf_makeresi(psf, set, prefs.recenter_flag, 0.2);

/* Store the chi's (sqrt(chi2) pdf close to Gaussian) */
  NFPRINTF(OUTPUT,"Computing Chi2 statistics...");
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
    locut = chimed - 4.0*chisig;
    hicut = chimed + 4.0*chisig;
    }

  free(chi);
/*
  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "<Chi2/dof> = %.3f\n",chivar/(nsample-(nsample>1?1:0)));
*/
  chi2 = chivar/(nsample-(nsample>1?1:0));

/* Clip outliers */
  NFPRINTF(OUTPUT,"Filtering PSF-candidates...");
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
VERSION	12/11/2007
 ***/
double	psf_chi2(psfstruct *psf, setstruct *set)
  {
   samplestruct	*sample;
   double	chi2;
   int		n, nsample;

/* First compute residuals for each sample (chi^2) */
  NFPRINTF(OUTPUT,"Computing residuals...");
  psf_makeresi(psf, set, prefs.recenter_flag, 0.0);

/* Store the chi's (sqrt(chi2) pdf close to gaussian) */
  NFPRINTF(OUTPUT,"Computing Chi2 statistics...");
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
PROTO	psfstruct *psf_init(int *dim, int ndim)
PURPOSE	Allocate and initialize a PSF structure.
INPUT	1D array of degrees of the polynom,
	Array of char pointers to the context names,
	Number of dimensions.
OUTPUT  psfstruct pointer.
NOTES   The maximum degrees and number of dimensions allowed are set in poly.h.
AUTHOR  E. Bertin (IAP)
VERSION 01/03/2007
 ***/
psfstruct	*psf_init(char **names, int *group, int ndim,
			int *dim, int ngroup,
			int wpsf, int hpsf, float psfstep, int nsample)
  {
   psfstruct	*psf;
   static char	str[MAXCHAR];
   char		**names2, **names2t;
   int		*group2, *dim2,
		d, ndim2,ngroup2, npix, nsnap;

/* Allocate memory for the PSF structure itself */
  QCALLOC(psf, psfstruct, 1);
  psf->dim = PSF_NMASKDIM;	/* This is constant */
  QMALLOC(psf->size, int, psf->dim);

/* The polynom */
  names2 = NULL;
  group2 = dim2 = NULL;
  if ((ndim2=ndim))
    {
    QMEMCPY(names, names2, char *, ndim);
    QMEMCPY(group, group2, int, ndim);
    }
  if ((ngroup2=ngroup))
    {
    QMEMCPY(dim, dim2, int, ngroup);
    }

  psf->poly = poly_init(group2, ndim2, dim2, ngroup2);

/*-- Compute the maximum advised number of degrees of freedom */
  if (ngroup2)
    while (psf->poly->ncoeff>nsample
	|| (int)(psf->poly->ncoeff/(psfstep*psfstep*PSF_FREEDFACTOR)+0.499)
	   >nsample)
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
  psf->npix = psf->size[0] = wpsf;
  psf->npix *= (psf->size[1] = hpsf);
  psf->npix *= (psf->size[2] = psf->poly->ncoeff);
  QCALLOC(psf->comp, float, psf->npix);
  npix = psf->size[0]*psf->size[1];
  QCALLOC(psf->loc, float, npix);
  QCALLOC(psf->resi, float, npix);

/* Context arrays */
  nsnap = 1;
  if (ndim2)
    {
    QMALLOC(psf->contextoffset, double, ndim2);
    QMALLOC(psf->contextscale, double, ndim2);
    QMALLOC(psf->contextname, char *, ndim2);
    for (names2t=names2, d=0; d<ndim2; d++)
      {
      nsnap *= prefs.context_nsnap;
      QMALLOC(psf->contextname[d], char, 80);
      strcpy(psf->contextname[d], *(names2t++));
      }
    }

/* Allocate an array of Moffat function fits */
  QMALLOC(psf->moffat, moffatstruct, nsnap);

/* Free temporary arrays */
  if (ndim)
    {
    free(names2);
    free(group2);
    free(dim2);
    }

 return psf;
  }


/****** psf_end ***************************************************************
PROTO   void psf_end(psfstruct *psf)
PURPOSE Free a PSF structure and everything it contains.
INPUT   psfstruct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/12/2007
 ***/
void	psf_end(psfstruct *psf)
  {
   int	d, ndim;

  ndim = psf->poly->ndim;
  for (d=0; d<ndim; d++)
    free(psf->contextname[d]);
  free(psf->contextname);
  poly_end(psf->poly);
  free(psf->pixmask);
  free(psf->basis);
  free(psf->comp);
  free(psf->loc);
  free(psf->resi);
  free(psf->size);
  free(psf->moffat);
  free(psf->homo_kernel);
  free(psf);

  return;
  }


/****** psf_make **************************************************************
PROTO	void	psf_make(psfstruct *psf, setstruct *set)
PURPOSE	Make the PSF.
INPUT	Pointer to the PSF,
	Pointer to the sample set.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2007
 ***/
void	psf_make(psfstruct *psf, setstruct *set)
  {
   polystruct	*poly;
   samplestruct	*sample;
   double	*pstack,*wstack, *basis, *pix,*wpix, *coeff, *pos, *post;
   float	*comp;
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
  QMALLOC(pstack, double, nsample);
  QMALLOC(wstack, double, nsample);
  QMALLOC(pos, double, poly->ndim?(nsample*poly->ndim):1);
  QMALLOC(basis, double, poly->ncoeff*nsample);

  for (sample=set->sample, post=pos, n=nsample; n--; sample++)
    {
    update_retina(set, sample, psf->pixstep);
    for (i=0; i<poly->ndim; i++)
      *(post++) = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    }

/* Make a polynomial fit to each pixel */
  for (i=0; i<npix; i++)
    {
/*-- Stack ith pixel from each PSF candidate */
    for (sample=set->sample, pix=pstack,wpix=wstack, n=nsample; n--; sample++)
      {
      *(pix++) = (double)*(sample->retina+i);
      *(wpix++) = (double)*(sample->retiweight+i);
      }

/*-- Polynomial fitting */
    poly_fit(poly, i?NULL:pos, pstack, wstack, nsample, basis);

/*-- Store as a PSF component */
    for (coeff=poly->coeff, comp=psf->comp+i,  c=ncoeff; c--; comp+=npix)
      *comp = *(coeff++);
    }

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
VERSION 12/11/2007
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
    for (p=npix; p--;)
      *(pl++) +=  fac**(ppc++);
    }

  return;
  }


/****** psf_makeresi **********************************************************
PROTO	void	psf_makeresi(psfstruct *psf, setstruct *set, int centflag,
		float psf_extraccu)
PURPOSE	Compute PSF residuals.
INPUT	Pointer to the PSF,
	Pointer to the sample set,
	Re-centering flag (0=no),
	PSF accuracy parameter.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2007
 ***/
void	psf_makeresi(psfstruct *psf, setstruct *set, int centflag,
		float psf_extraccu)
  {
   samplestruct		*sample;
   static double	pos[MAXCONTEXT], amat[9], bmat[3];
   double		*dresi, *dresit, *amatt,
			*cvigx,*cvigxt, *cvigy,*cvigyt,
			nm1, chi2, dx,dy, ddx,ddy, dval,dvalx,dvaly,dwval,
			radmin2,radmax2, hcw,hch, yb, mx2,my2,mxy,
			xc,yc,rmax2,x,y, mse, xi2, xyi;
   float		*vigresi, *vig, *vigw, *fresi,*fresit,
			*cbasis,*cbasist, *cdata,*cdatat, *cvigw,*cvigwt,
			norm, fval, vigstep, psf_extraccu2, wval;
   int			i,j,n,ix,iy, ndim,npix,nsample, cw,ch,ncpix, okflag,
			accuflag, nchi2;

  accuflag = (psf_extraccu > 1.0/BIG);
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
        clapack_dpotrf(CblasRowMajor, CblasUpper, 3, amat, 3);
        clapack_dpotrs(CblasRowMajor, CblasUpper, 3, 1, amat, 3, bmat, 3);

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
    chi2 = 0.0;
    dresit = dresi;
    psf_extraccu2 = psf_extraccu*psf_extraccu*norm*norm;
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
    mse = 0.0;
    for (iy=set->vigsize[1]; iy--; y+=1.0)
      {
      x = -xc;
      for (ix=set->vigsize[0]; ix--; x+=1.0, vig++, vigresi++, dresit++)
        if ((wval=*(vigw++))>0.0)
          {
          if (accuflag)
            wval = 1.0/(1.0 / wval + psf_extraccu2**vigresi**vigresi);
          *vigresi = fval = (*vig-*vigresi*norm);
          if (x*x+y*y<rmax2)
            {
            mse += fval*fval;
            nchi2++;
            chi2 += (double)(wval*fval*fval);
            *dresit += fval;
            }
          }
      }

    sample->chi2 = (nchi2> 1)? chi2/(nchi2-1) : chi2;
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
PROTO	void	psf_refine(psfstruct *psf, setstruct *set)
PURPOSE	Refine PSF by solving a system to recover "aliased" components.
INPUT	Pointer to the PSF,
	Pointer to the sample set.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 20/11/2007
 ***/
void	psf_refine(psfstruct *psf, setstruct *set)
  {
   polystruct		*poly;
   samplestruct		*sample;
   double		pos[MAXCONTEXT];
   char			str[MAXCHAR];
   double		*desmat,*desmatt,*desmatt2, *desmat0,*desmat02,
			*bmat,*bmatt, *basis,*basist, *basist2,
			*sigvig,*sigvigt, *alphamat,*alphamatt,
			*betamat,*betamatt, *coeffmat,*coeffmatt,
			dx,dy, dval, norm, tikfac;
   float		*vig,*vigt,*vigt2, *wvig,
			*vecvig,*vecvigt, *ppix, *vec,
			vigstep;
   int			*desindex,*desindext,*desindext2,
			*desindex0,*desindex02;
   int			i,j,jo,k,l,c,n, npix,nvpix, ndata,ncoeff,nsample,npsf,
			ncontext, nunknown, matoffset, dindex;

/* Exit if no pixel is to be "refined" or if no sample is available */
  if (!set->nsample || !psf->basis)
    return;

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

  NFPRINTF(OUTPUT,"Processing samples...");
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

/* Go through each sample */
  for (sample=set->sample, n=0; n<nsample ; n++, sample++)
    {
    sprintf(str, "Processing sample #%d", n+1);
    NFPRINTF(OUTPUT, str);
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
    basis = poly->basis;
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
    tikfac= 0.008;
    tikfac = 1.0/(tikfac*tikfac);
    for (i=0; i<nunknown; i++)
      alphamat[i+nunknown*i] += tikfac;
    }

  NFPRINTF(OUTPUT,"Solving the system...");

  clapack_dpotrf(CblasRowMajor, CblasUpper, nunknown, alphamat, nunknown);
  clapack_dpotrs(CblasRowMajor, CblasUpper, nunknown, 1, alphamat, nunknown,
	betamat, nunknown);

  NFPRINTF(OUTPUT,"Updating the PSF...");
  if (!psf->pixmask)
    memset(psf->comp, 0, npix*ncoeff*sizeof(float));
  betamatt = betamat;
  for (j=0; j<npsf; j++)
    {
    ppix = psf->comp;
    for (c=ncoeff; c--;)
      {
      vec = &psf->basis[j*npix];
      dval = *(betamatt++);
      for (i=npix; i--;)
        *(ppix++) += dval**(vec++);
      }
    }

/* Free all */
  free(alphamat);
  free(betamat);

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
VERSION 20/11/2007
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
      NFPRINTF(OUTPUT,"Selecting pixels...");
      psforder = (float *)NULL;		/* To avoid gcc -Wall warnings */
      QMEMCPY(psf->comp, psforder, float, npix);
      for (psfordert=psforder, i=npix; i--; psfordert++)
        *psfordert = fabs(*psfordert);
      hmedian(psforder, npix);
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
		*(1+(int)(INTERPH*psf->pixstep))+1;
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
      NFPRINTF(OUTPUT, str);
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


/****** psf_save **************************************************************
PROTO   void	psf_save(psfstruct *psf, char *filename, int ext, int next)
PURPOSE Save the PSF data as a FITS file.
INPUT   Pointer to the PSF structure,
	Filename,
	Extension number,
	Number of extensions.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2007
 ***/
void	psf_save(psfstruct *psf, char *filename, int ext, int next)
  {
   static catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   char		*head, str[80];
   int		i, temp;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  if (!ext)
    {
    cat = new_cat(1);
    init_cat(cat);
    strcpy(cat->filename, filename);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);
    save_tab(cat, cat->tab);
    }
  tab = new_tab("PSF_DATA");

  head = tab->headbuf;
  addkeywordto_head(tab, "LOADED", "Number of loaded sources");
  fitswrite(head, "LOADED", &psf->samples_loaded, H_INT, T_LONG);
  addkeywordto_head(tab, "ACCEPTED", "Number of accepted sources");
  fitswrite(head, "ACCEPTED", &psf->samples_accepted, H_INT, T_LONG);
  addkeywordto_head(tab, "CHI2", "Final Chi2");
  fitswrite(head, "CHI2", &psf->chi2, H_FLOAT, T_DOUBLE);
  addkeywordto_head(tab, "POLNAXIS", "Number of context parameters");
  fitswrite(head, "POLNAXIS", &psf->poly->ndim, H_INT, T_LONG);
  for (i=0; i<psf->poly->ndim; i++)
    {
    sprintf(str, "POLGRP%1d", i+1);
    addkeywordto_head(tab, str, "Polynom group for this context parameter");
    temp = psf->poly->group[i]+1;
    fitswrite(head, str, &temp, H_INT, T_LONG);
    sprintf(str, "POLNAME%1d", i+1);
    addkeywordto_head(tab, str, "Name of this context parameter");
    fitswrite(head, str, psf->contextname[i], H_STRING, T_STRING);
    sprintf(str, "POLZERO%1d", i+1);
    addkeywordto_head(tab, str, "Offset value for this context parameter");
    fitswrite(head, str, &psf->contextoffset[i], H_EXPO, T_DOUBLE);
    sprintf(str, "POLSCAL%1d", i+1);
    addkeywordto_head(tab, str, "Scale value for this context parameter");
    fitswrite(head, str, &psf->contextscale[i], H_EXPO, T_DOUBLE);
    }

  addkeywordto_head(tab, "POLNGRP", "Number of context groups");
  fitswrite(head, "POLNGRP", &psf->poly->ngroup, H_INT, T_LONG);
  for (i=0; i<psf->poly->ngroup; i++)
    {
    sprintf(str, "POLDEG%1d", i+1);
    addkeywordto_head(tab, str, "Polynom degree for this context group");
    fitswrite(head, str, &psf->poly->degree[i], H_INT, T_LONG);
    }

/* Add and write important scalars as FITS keywords */
  /* -- FM -- : write fwhm too */
  addkeywordto_head(tab, "PSF_FWHM", "PSF FWHM");
  fitswrite(head, "PSF_FWHM", &psf->fwhm, H_FLOAT, T_FLOAT);
  addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data");
  fitswrite(head, "PSF_SAMP", &psf->pixstep, H_FLOAT, T_FLOAT);
  addkeywordto_head(tab, "PSFNAXIS", "Dimensionality of the PSF data");
  fitswrite(head, "PSFNAXIS", &psf->dim, H_INT, T_LONG);
  for (i=0; i<psf->dim; i++)
    {
    sprintf(str, "PSFAXIS%1d", i+1);
    addkeywordto_head(tab, str, "Number of element along this axis");
    fitswrite(head, str, &psf->size[i], H_INT, T_LONG);
    }

/* Create and fill the arrays */
  key = new_key("PSF_MASK");
  key->naxis = psf->dim;
  QMALLOC(key->naxisn, int, key->naxis);
  for (i=0; i<psf->dim; i++)
    key->naxisn[i] = psf->size[i];
  strcat(key->comment, "Tabulated PSF data");
  key->htype = H_FLOAT;
  key->ttype = T_FLOAT;
  key->nbytes = psf->npix*t_size[T_FLOAT];
  key->nobj = 1;
  key->ptr = psf->comp;
  add_key(key, tab, 0);

  save_tab(cat, tab);
/* But don't touch my arrays!! */
  blank_keys(tab);
  free_tab(tab);

  if (ext==next-1)
    free_cat(&cat , 1);

  return;
  }

