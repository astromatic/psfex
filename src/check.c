  /*
 				check.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Production of check-images for the PSF.
*
*	Last modify:	15/11/2007
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
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"check.h"
#include	"diagnostic.h"
#include	"poly.h"
#include	"psf.h"
#include	"sample.h"
#include	"vignet.h"


/****** psf_writecheck ********************************************************
PROTO	void	psf_writecheck(psfstruct *psf, setstruct *set, char *filename,
		checkenum checktype, int ext, int next, int cubeflag)
PURPOSE	Write a FITS image for check.
INPUT	Pointer to the PSF,
	Pointer to the sample set,
	Check-image filename,
	Check-image type,
	Extension number,
	Number of extensions,
	Datacube flag.
OUTPUT  -.
NOTES   Check-image is written as a datacube if cubeflag!=0.
AUTHOR  E. Bertin (IAP)
VERSION 15/11/2007
 ***/
void	psf_writecheck(psfstruct *psf, setstruct *set, char *filename,
		checkenum checktype, int ext, int next, int cubeflag)
  {
   static catstruct	*ccat[MAXCHECK];
   catstruct		*cat;
   tabstruct		*tab;
   samplestruct		*sample;
   char			str[82],
			*head;
   static double	dpos[POLY_MAXDIM], *dpost;
   double		dstep,dstart, dval1,dval2, scalefac;
   float		*pix,*pix0, *vig,*vig0, *fpix,*fpixsym,
			val;
   int			i,x,y, w,h,n, npc,nt,nr, nw,nh, step, ival1,ival2, npix;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  if (!ext)
    {
    cat = new_cat(1);
    init_cat(cat);
    strcpy(cat->filename, filename);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);
    if (next>1)
      {
      addkeywordto_head(cat->tab, "NEXTEND ", "Number of extensions");
      fitswrite(cat->tab->headbuf, "NEXTEND", &next, H_INT, T_LONG);
      save_tab(cat, cat->tab);
      }
    ccat[checktype] = cat;
    }
  else
    cat = ccat[checktype];

  sprintf(str, "chip%02d", ext+1);

  tab = new_tab(str);
  head = tab->headbuf;
  tab->bitpix =  BP_FLOAT;
  tab->bytepix = t_size[T_FLOAT];

  switch(checktype)
    {
    case PSF_BASIS:
/*----  View basis vectors as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = psf->size[0];
        tab->naxisn[1] = psf->size[1];
        tab->naxisn[2] = psf->nbasis;
        npix = tab->naxisn[0]*tab->naxisn[1]*tab->naxisn[2];
        tab->tabsize = tab->bytepix*npix;
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        fpix = psf->basis;
        for (i=npix; i--;)
          *(pix++) = *(fpix++);
        }
      else
        {
        nw = (int)sqrt((double)psf->nbasis);
        nw = ((nw-1)/10+1)*10;
        nh = (psf->nbasis-1)/nw + 1;
        w = psf->size[0];
        h = psf->dim>1? psf->size[1] : 1;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        step = (nw-1)*w;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        fpix = psf->basis;
        for (n=0; n<psf->nbasis; n++)
          {
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
      break;
    case PSF_CHI:
/*---- sqrt(chi2) map in PSF pixel-space */
      nw = 1;
      nh = 1;
      w = psf->size[0];
      h = psf->size[1];
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      pix = pix0;
      fpix = psf->resi;
      for (i=w*h; i--;)
        *(pix++) = *(fpix++);
      break;

    case PSF_PROTO:
/*---- PSF data for all components are arranged as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = psf->size[0];
        tab->naxisn[1] = psf->size[1];
        tab->naxisn[2] = psf->size[2];
        npix = tab->naxisn[0]*tab->naxisn[1]*tab->naxisn[2];
        tab->tabsize = tab->bytepix*npix;
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        fpix = psf->comp;
        for (i=npix; i--;)
          *(pix++) = *(fpix++);
        }
      else
        {
        npc = psf->size[2];
        nw = npc<10? npc:10;
        nh = (npc-1)/nw + 1;
        w = psf->size[0];
        h = psf->size[1];
        step = (nw-1)*w;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
/*------ Normalize the components in the image corners: pos=(0.5,0.5,..) */
        for (dpost=dpos, i=psf->poly->ndim; i--;)
          *(dpost++) = 0.5;
        poly_func(psf->poly, dpos);
        dpost = psf->poly->basis;
        fpix = psf->comp;
        for (n=0; n<npc; n++)
          {
          val = *(dpost++);
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++)*val;
          }
        }
      break;

    case PSF_RESIDUALS:
/*---- Residual vectors for all samples are arranged as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = set->vigsize[0];
        tab->naxisn[1] = set->vigsize[1];
        tab->naxisn[2] = set->nsample;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          fpix = (sample++)->vigresi;
          for (i=npix; i--;)
            *(pix++) = *(fpix++);
          }
        }
      else
        {
        nw = (int)sqrt((double)set->nsample);
        nw = ((nw-1)/10+1)*10;
        nh = (set->nsample-1)/nw + 1;
        w = set->vigsize[0];
        h = set->vigdim>1? set->vigsize[1] : 1;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        step = (nw-1)*w;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          fpix = (sample++)->vigresi;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
      break;

    case PSF_RAWDATA:
/*----  View original samples as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = set->vigsize[0];
        tab->naxisn[1] = set->vigsize[1];
        tab->naxisn[2] = set->nsample;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          fpix = (sample++)->vig;
          for (i=npix; i--;)
            *(pix++) = *(fpix++);
          }
        }
      else
        {
        nw = (int)sqrt((double)set->nsample);
        nw = ((nw-1)/10+1)*10;
        nh = (set->nsample-1)/nw + 1;
        w = set->vigsize[0];
        h = set->vigdim>1? set->vigsize[1] : 1;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        step = (nw-1)*w;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          fpix = (sample++)->vig;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
      break;

    case PSF_SAMPLES:
/*----  View all training samples as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = set->retisize[0];
        tab->naxisn[1] = set->retisize[1];
        tab->naxisn[2] = set->nsample;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          fpix = (sample++)->retina;
          for (i=npix; i--;)
            *(pix++) = *(fpix++);
          }
        }
      else
        {
        nw = (int)sqrt((double)set->nsample);
        nw = ((nw-1)/10+1)*10;
        nh = (set->nsample-1)/nw + 1;
        w = set->retisize[0];
        h = set->retidim>1? set->retisize[1] : 1;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        step = (nw-1)*w;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          fpix = (sample++)->retina;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
      break;

    case PSF_SNAPSHOTS:
/*----  View reconstructed PSFs as small vignets */
      npc = prefs.ncontext_group;
      nw = npc? prefs.context_nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= prefs.context_nsnap;
      nh = nt/nw;
      w = set->retisize[0];
      h = set->retidim>1? set->retisize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      dstep = 1.0/prefs.context_nsnap;
      dstart = (1.0-dstep)/2.0;
      memset(dpos, 0, POLY_MAXDIM*sizeof(double));
      for (i=0; i<npc; i++)
        dpos[i] = -dstart;
      for (n=0; n<nt; n++)
        {
        psf_build(psf, dpos);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        for (i=0; i<npc; i++)
          if (dpos[i]<dstart-0.01)
            {
            dpos[i] += dstep;
            break;
            }
          else
            dpos[i] = -dstart;
        }
      break;

    case PSF_SNAPSHOTS_IMRES:
/*----  View reconstructed PSFs as small vignets */
      npc = prefs.ncontext_group;
      nw = npc? prefs.context_nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= prefs.context_nsnap;
      nh = nt/nw;
      w = set->vigsize[0];
      h = set->vigdim>1? set->vigsize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QMALLOC(vig0, float, w*h);
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      dstep = 1.0/prefs.context_nsnap;
      dstart = (1.0-dstep)/2.0;
      memset(dpos, 0, POLY_MAXDIM*sizeof(double));
      for (i=0; i<npc; i++)
        dpos[i] = -dstart;
      for (n=0; n<nt; n++)
        {
        psf_build(psf, dpos);
        vignet_resample(psf->loc, psf->size[0], psf->size[1],
	vig0, set->vigsize[0], set->vigsize[1], 0.0, 0.0, 1.0/psf->pixstep, 1.0);
        vig = vig0;
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(vig++);
        for (i=0; i<npc; i++)
          if (dpos[i]<dstart-0.01)
            {
            dpos[i] += dstep;
            break;
            }
          else
            dpos[i] = -dstart;
        }
      free(vig0);
      break;

    case PSF_WEIGHTS:
/*----  View all training sample weights as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = set->retisize[0];
        tab->naxisn[1] = set->retisize[1];
        tab->naxisn[2] = set->nsample;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          fpix = (sample++)->retiweight;
          for (i=npix; i--;)
            *(pix++) = *(fpix++);
          }
        }
      else
        {
        nw = (int)sqrt((double)set->nsample);
        nw = ((nw-1)/10+1)*10;
        nh = (set->nsample-1)/nw + 1;
        w = set->retisize[0];
        h = set->retidim>1? set->retisize[1] : 1;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        step = (nw-1)*w;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          fpix = (sample++)->retiweight;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
      break;

    case PSF_MOFFAT:
/*----  View reconstructed PSFs as Moffat fits */
      npc = prefs.ncontext_group;
      nw = npc? prefs.context_nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= prefs.context_nsnap;
      nh = nt/nw;
      w = set->retisize[0];
      h = set->retidim>1? set->retisize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      for (nr=1, i=psf->poly->ndim; (i--)>0;)
        nr *= prefs.context_nsnap;	/* nr is the true number of Moffats */
      for (n=0; n<nt; n++)
        {
        psf_moffat(psf, &psf->moffat[n%nr]);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        }
      break;
 
   case PSF_SUBMOFFAT:
/*----  View the difference between reconstructed PSFs and Moffat fits */
      npc = prefs.ncontext_group;
      nw = npc? prefs.context_nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= prefs.context_nsnap;
      nh = nt/nw;
      w = set->retisize[0];
      h = set->retidim>1? set->retisize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      for (nr=1, i=psf->poly->ndim; (i--)>0;)
        nr *= prefs.context_nsnap;	/* nr is the true number of Moffats */
      for (n=0; n<nt; n++)
        {
        psf_build(psf, dpos);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        psf_moffat(psf, &psf->moffat[n%nr]);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) -= *(fpix++);
        }
      break;
   case PSF_SUBSYM:
/*----  View the difference between reconstructed PSFs and their symmetrical */
      npc = prefs.ncontext_group;
      nw = npc? prefs.context_nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= prefs.context_nsnap;
      nh = nt/nw;
      w = set->retisize[0];
      h = set->retidim>1? set->retisize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      for (n=0; n<nt; n++)
        {
        psf_build(psf, dpos);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        fpixsym = psf->loc + w*h;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++)-*(--fpixsym);
        }
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Yet unavailable CHECKIMAGE type",
	"");
    }

/* Put some (dummy) WCS in the output check-images */
  fitsread(set->head, "NAXIS1  ", &ival1, H_INT, T_LONG);
  fitsread(set->head, "NAXIS2  ", &ival2, H_INT, T_LONG);
  if (fitsread(set->head, "CTYPE1  ", str, H_STRING, T_STRING)==RETURN_OK)
    {
    addkeywordto_head(tab,"CTYPE1  ","WCS axis type");
    fitswrite(tab->headbuf, "CTYPE1  ", str,  H_STRING, T_STRING);
    if (fitsread(set->head, "CTYPE2  ", str, H_STRING, T_STRING)==RETURN_OK)
      {
      addkeywordto_head(tab,"CTYPE2  ","WCS axis type");
      fitswrite(tab->headbuf, "CTYPE2  ", str,  H_STRING, T_STRING);
      }
    if (cubeflag)
      {
      addkeywordto_head(tab,"CTYPE3  ","WCS axis type");
      fitswrite(tab->headbuf, "CTYPE3  ", " ",  H_STRING, T_STRING);
      }
    }
  if (fitsread(set->head, "CRVAL1  ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    addkeywordto_head(tab, "CRVAL1  ",
	"WCS coordinates of the reference pixel");
    fitswrite(tab->headbuf, "CRVAL1  ", &dval1,  H_EXPO, T_DOUBLE);
    if (fitsread(set->head, "CRVAL2  ", &dval2, H_EXPO, T_DOUBLE)==RETURN_OK)
      {
      addkeywordto_head(tab, "CRVAL2  ",
	"WCS coordinates of the reference pixel");
      fitswrite(tab->headbuf, "CRVAL2  ", &dval2,  H_EXPO, T_DOUBLE);
      }
    if (cubeflag)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CRVAL3  ",
	"WCS coordinates of the reference pixel");
      fitswrite(tab->headbuf, "CRVAL3  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    }
  scalefac = ((dval1=ival1/(double)tab->naxisn[0])
	< (dval2=ival2/(double)tab->naxisn[1])) ? dval1 : dval2;
  if (fitsread(set->head, "CD1_1   ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    dval1 *= scalefac;
    addkeywordto_head(tab, "CD1_1   ", "WCS transformation matrix");
    fitswrite(tab->headbuf, "CD1_1   ", &dval1,  H_EXPO, T_DOUBLE);
    if (fitsread(set->head, "CD1_2   ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
      {
      dval1 *= scalefac;
      addkeywordto_head(tab, "CD1_2   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD1_2   ", &dval1,  H_EXPO, T_DOUBLE);
      }
    if (fitsread(set->head, "CD2_1   ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
      {
      dval1 *= scalefac;
      addkeywordto_head(tab, "CD2_1   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD2_1   ", &dval1,  H_EXPO, T_DOUBLE);
      }
    if (fitsread(set->head, "CD2_2   ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
      {
      dval1 *= scalefac;
      addkeywordto_head(tab, "CD2_2   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD2_2   ", &dval1,  H_EXPO, T_DOUBLE);
      }
    if (cubeflag)
      {
      dval1 = 0.0;
      addkeywordto_head(tab, "CD1_3   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD1_3   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD2_3   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD2_3   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD3_1   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD3_1   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD3_2   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD3_2   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 1.0;
      addkeywordto_head(tab, "CD3_3   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD3_3   ", &dval1,  H_EXPO, T_DOUBLE);
      }
    }
  if (fitsread(set->head, "CDELT1  ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    dval1 *= scalefac;
    addkeywordto_head(tab, "CDELT1  ", "WCS pixel scale");
    fitswrite(tab->headbuf, "CDELT1  ", &dval1,  H_EXPO, T_DOUBLE);
    if (fitsread(set->head, "CDELT2  ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
      {
      dval1 *= scalefac;
      addkeywordto_head(tab, "CDELT2  ", "WCS pixel scale");
      fitswrite(tab->headbuf, "CDELT2  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    if (cubeflag)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CDELT3  ","WCS pixel scale");
      fitswrite(tab->headbuf, "CDELT3  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    }
  if (fitsread(set->head, "CRPIX1  ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    addkeywordto_head(tab, "CRPIX1  ",
	"pixel coordinates of the reference pixel");
    dval1 = (dval1 - (ival1+1)/2)/scalefac + (tab->naxisn[0]+1)/2.0;
    fitswrite(tab->headbuf, "CRPIX1  ", &dval1,  H_EXPO, T_DOUBLE);
    if (fitsread(set->head, "CRPIX2  ", &dval2, H_EXPO, T_DOUBLE)==RETURN_OK)
      {
      addkeywordto_head(tab, "CRPIX2  ",
	"pixel coordinates of the reference pixel");
      dval2 = (dval2 - (ival2+1)/2)/scalefac + (tab->naxisn[1]+1)/2.0;
      fitswrite(tab->headbuf, "CRPIX2  ", &dval2,  H_EXPO, T_DOUBLE);
      }
    if (cubeflag)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CRPIX3  ",
	"pixel coordinates of the reference pixel");
      fitswrite(tab->headbuf, "CRPIX3  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    }

  if (next == 1)
    prim_head(tab);
  fitswrite(head, "XTENSION", "IMAGE   ", H_STRING, T_STRING);
/* save table */
  save_tab(cat, tab);
  free_tab(tab);
  if (ext==next-1)
    free_cat(&cat, 1);

  return;
  }

