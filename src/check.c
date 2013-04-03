/*
*				check.c
*
* Produce check-images of the PSF
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
#include	"check.h"
#include	"diagnostic.h"
#include	"field.h"
#include	"wcs/poly.h"
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"vignet.h"


/****** check_write ********************************************************
PROTO	void	check_write(fieldstruct *field, char *checkname,
		checkenum checktype, int ext, int next, int cubeflag)
PURPOSE	Write a FITS image for check.
INPUT	Pointer to the field,
	Pointer to the set,
	Check-image filename,
	Check-image type,
	Extension number,
	Number of extensions,
	Datacube flag.
OUTPUT  -.
NOTES   Check-image is written as a datacube if cubeflag!=0.
AUTHOR  E. Bertin (IAP)
VERSION 19/07/2012
 ***/
void	check_write(fieldstruct *field, setstruct *set, char *checkname,
		checkenum checktype, int ext, int next, int cubeflag)
  {
   psfstruct		*psf;
   catstruct		*cat;
   tabstruct		*tab;
   samplestruct		**gridsample,
			*sample,*osample;
   char			filename[MAXCHAR], str[82],
			*head, *pstr,*pstr2;
   static double	dpos[POLY_MAXDIM], *dpost;
   double		dstep,dstart, dval1,dval2, scalefac, dstepx,dstepy;
   float		*pix,*pix0, *vig,*vig0, *fpix,*fpixsym,
			val;
   int			i,j,l,x,y, w,h,n, npc,nt, nw,nh,np, npos,npos2,
			step,step2, ipos, inpos, ival1,ival2, npix;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  if (!ext)
    {
    cat = new_cat(1);
    init_cat(cat);
    strcpy(cat->filename, checkname);
    if (!(pstr = strrchr(cat->filename, '.')))
      pstr = cat->filename+strlen(cat->filename);
    strcpy(filename, field->rcatname);
    if (!(pstr2 = strrchr(filename, '.')))
      pstr2 = filename+strlen(filename);
    *pstr2 = '\0';
    sprintf(pstr, "_%s.fits", filename);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", cat->filename);
    if (next>1)
      {
      addkeywordto_head(cat->tab, "NEXTEND ", "Number of extensions");
      fitswrite(cat->tab->headbuf, "NEXTEND", &next, H_INT, T_LONG);
      save_tab(cat, cat->tab);
      }
    field->ccat[checktype] = cat;
    }
  else
    cat = field->ccat[checktype];

  sprintf(str, "chip%02d", ext+1);

  psf = field->psf[ext];
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
/*---- sqrt(chi2) map in PSF pixel-space 
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
*/
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = set->vigsize[0];
        tab->naxisn[1] = set->vigsize[1];
        tab->naxisn[2] = set->nsample? set->nsample : 1;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          fpix = (sample++)->vigchi;
          for (i=npix; i--;)
            *(pix++) = *(fpix++);
          }
        }
      else
        {
        nw = (int)sqrt((double)set->nsample);
        nw = nw? ((nw-1)/10+1)*10 : 1;
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
          fpix = (sample++)->vigchi;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
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
        tab->naxisn[2] = set->nsample? set->nsample : 1;
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
        nw = nw? ((nw-1)/10+1)*10 : 1;
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

    case PSF_RESIDUALS_GRID:
/*----  View the brightest samples arranged as a grid */
      npc = 2;
      npos = psf->nsnap;
      for (nt=1, i=npc; (i--)>0;)
        nt *= psf->nsnap;
      QCALLOC(gridsample, samplestruct *, nt);
      dstepx = (float)prefs.context_nsnap / field->wcs[ext]->naxisn[0];
      dstepy = (float)prefs.context_nsnap / field->wcs[ext]->naxisn[1];
      sample = set->sample;
      for (n=set->nsample; n--; sample++)
        {
        ipos = (int)(dstepx * (sample->x+0.5)) ;
        if (ipos<0)
          ipos = 0;
        else if (ipos>=npos)
          ipos = npos-1; 
        inpos = ipos;
        ipos = (int)(dstepy * (sample->y+0.5));
        if (ipos<0)
          ipos = 0;
        else if (ipos>=npos)
          ipos = npos-1; 
        inpos += ipos*npos;
        osample = gridsample[inpos];
        if (!osample || osample->norm < sample->norm)
          gridsample[inpos] = sample;
        }

      nw = npc? psf->nsnap : 1;
      w = set->vigsize[0];
      h = set->vigsize[1];
      if (cubeflag)
        {
        nh = npc>2? psf->nsnap : nt/nw;
        np = npc>2? nt/(nw*nh) : 1;
        tab->naxis = 4;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = w;
        tab->naxisn[1] = h;
        tab->naxisn[2] = nw;
        tab->naxisn[3] = nh;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2]*tab->naxisn[3];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        for (n=0; n<nt; n++)
          {
          if ((sample = gridsample[n]))
            memcpy(pix, sample->vigresi, npix*sizeof(float));
          pix += npix;
          }
        }
      else
        {
        nh = nt/nw;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        step = (nw-1)*w;
        for (n=0; n<nt; n++)
          {
          if ((sample = gridsample[n]))
            {
            pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
            fpix = sample->vigresi;
            for (y=h; y--; pix += step)
              for (x=w; x--;)
                *(pix++) = *(fpix++);
            }
          }
        }
      free(gridsample);
      break;

    case PSF_SAMPLES:
/*----  View original samples as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = set->vigsize[0];
        tab->naxisn[1] = set->vigsize[1];
        tab->naxisn[2] = set->nsample? set->nsample : 1;
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
        nw = nw? ((nw-1)/10+1)*10 : 1;
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

    case PSF_SAMPLES_GRID:
/*----  View the brightest samples arranged as a grid */
      npc = 2;
      npos = psf->nsnap;
      for (nt=1, i=npc; (i--)>0;)
        nt *= psf->nsnap;
      QCALLOC(gridsample, samplestruct *, nt);
      dstepx = (float)prefs.context_nsnap / field->wcs[ext]->naxisn[0];
      dstepy = (float)prefs.context_nsnap / field->wcs[ext]->naxisn[1];
      sample = set->sample;
      for (n=set->nsample; n--; sample++)
        {
        ipos = (int)(dstepx * (sample->x+0.5)) ;
        if (ipos<0)
          ipos = 0;
        else if (ipos>=npos)
          ipos = npos-1; 
        inpos = ipos;
        ipos = (int)(dstepy * (sample->y+0.5));
        if (ipos<0)
          ipos = 0;
        else if (ipos>=npos)
          ipos = npos-1; 
        inpos += ipos*npos;
        osample = gridsample[inpos];
        if (!osample || osample->norm < sample->norm)
          gridsample[inpos] = sample;
        }

      nw = npc? psf->nsnap : 1;
      w = set->vigsize[0];
      h = set->vigsize[1];
      if (cubeflag)
        {
        nh = npc>2? psf->nsnap : nt/nw;
        np = npc>2? nt/(nw*nh) : 1;
        tab->naxis = 4;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = w;
        tab->naxisn[1] = h;
        tab->naxisn[2] = nw;
        tab->naxisn[3] = nh;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2]*tab->naxisn[3];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        for (n=0; n<nt; n++)
          {
          if ((sample = gridsample[n]))
            memcpy(pix, sample->vig, npix*sizeof(float));
          pix += npix;
          }
        }
      else
        {
        nh = nt/nw;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        step = (nw-1)*w;
        for (n=0; n<nt; n++)
          {
          if ((sample = gridsample[n]))
            {
            pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
            fpix = sample->vig;
            for (y=h; y--; pix += step)
              for (x=w; x--;)
                *(pix++) = *(fpix++);
            }
          }
        }
      free(gridsample);
      break;

    case PSF_WEIGHTS:
/*----  View original weights as small vignets */
      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = set->vigsize[0];
        tab->naxisn[1] = set->vigsize[1];
        tab->naxisn[2] = set->nsample? set->nsample : 1;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        sample = set->sample;
        for (n=0; n<set->nsample; n++)
          {
          fpix = (sample++)->vigweight;
          for (i=npix; i--;)
            *(pix++) = *(fpix++);
          }
        }
      else
        {
        nw = (int)sqrt((double)set->nsample);
        nw = nw? ((nw-1)/10+1)*10 : 1;
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
          fpix = (sample++)->vigweight;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
      break;

    case PSF_SNAPSHOTS:
/*----  View reconstructed PSFs as small vignets */
      npc = psf->poly->ndim;
      nw = npc? psf->nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= psf->nsnap;
      w = psf->size[0];
      h = psf->dim>1? psf->size[1] : 1;
      memset(dpos, 0, POLY_MAXDIM*sizeof(double));
      dstep = 1.0/prefs.context_nsnap;
      dstart = (1.0-dstep)/2.0;
      for (i=0; i<npc; i++)
        dpos[i] = -dstart;
      if (cubeflag)
        {
        nh = npc>2? psf->nsnap : nt/nw;
        np = npc>2? nt/(nw*nh) : 1;
        tab->naxis = 4;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = w;
        tab->naxisn[1] = h;
        tab->naxisn[2] = nw;
        tab->naxisn[3] = nh;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2]*tab->naxisn[3];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        for (n=0; n<nt; n++)
          {
          psf_build(psf, dpos);
          fpix = psf->loc;
          for (i=npix; i--;)
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
        }
      else
        {
        nh = nt/nw;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        step = (nw-1)*w;
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
        }
      break;

    case PSF_SNAPSHOTS_IMRES:
/*----  View reconstructed PSFs as small vignets */
      npc = psf->poly->ndim;
      nw = npc? psf->nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= psf->nsnap;
      w = set->vigsize[0];
      h = set->vigdim>1? set->vigsize[1] : 1;
      QMALLOC(vig0, float, w*h);
      memset(dpos, 0, POLY_MAXDIM*sizeof(double));
      dstep = 1.0/prefs.context_nsnap;
      dstart = (1.0-dstep)/2.0;
      for (i=0; i<npc; i++)
        dpos[i] = -dstart;
      if (cubeflag)
        {
        nh = npc>2? psf->nsnap : nt/nw;
        tab->naxis = 4;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = w;
        tab->naxisn[1] = h;
        tab->naxisn[2] = nw;
        tab->naxisn[3] = nh;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2]*tab->naxisn[3];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        for (n=0; n<nt; n++)
          {
          psf_build(psf, dpos);
          vignet_resample(psf->loc, psf->size[0], psf->size[1],
		vig0, set->vigsize[0], set->vigsize[1], 0.0, 0.0,
		1.0/psf->pixstep, 1.0);
          vig = vig0;
          for (i=npix; i--;)
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
        }
      else
        {
        nh = nt/nw;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        step = (nw-1)*w;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        for (n=0; n<nt; n++)
          {
          psf_build(psf, dpos);
          vignet_resample(psf->loc, psf->size[0], psf->size[1],
		vig0, set->vigsize[0], set->vigsize[1], 0.0, 0.0,
		1.0/psf->pixstep, 1.0);
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
        }
      free(vig0);
      break;

    case PSF_MOFFAT:
/*----  View reconstructed PSFs as Moffat fits */
      npc = psf->poly->ndim;
      nw = npc? psf->nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= psf->nsnap;
      nh = nt/nw;
      w = psf->size[0];
      h = psf->dim>1? psf->size[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      dstep = 1.0/psf->nsnap;
      dstart = (1.0-dstep)/2.0;
      memset(dpos, 0, POLY_MAXDIM*sizeof(double));
      for (i=0; i<npc; i++)
        dpos[i] = -dstart;
      for (n=0; n<nt; n++)
        {
        psf_moffat(psf, &psf->pfmoffat[n]);
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
 
   case PSF_SUBMOFFAT:
/*----  View the difference between reconstructed PSFs and Moffat fits */
      npc = psf->poly->ndim;
      nw = npc? psf->nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= psf->nsnap;
      nh = nt/nw;
      w = psf->size[0];
      h = psf->dim>1? psf->size[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      dstep = 1.0/psf->nsnap;
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
        psf_moffat(psf, &psf->pfmoffat[n]);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) -= *(fpix++);
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
   case PSF_SUBSYM:
/*----  View the difference between reconstructed PSFs and their symmetrical */
      npc = psf->poly->ndim;
      nw = npc? psf->nsnap : 1;
      for (nt=1, i=npc; (i--)>0;)
        nt *= psf->nsnap;
      nh = nt/nw;
      w = psf->size[0];
      h = psf->dim>1? psf->size[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      dstep = 1.0/psf->nsnap;
      dstart = (1.0-dstep)/2.0;
      memset(dpos, 0, POLY_MAXDIM*sizeof(double));
      for (i=0; i<npc; i++)
        dpos[i] = -dstart;
      for (n=0; n<nt; n++)
        {
        psf_build(psf, dpos);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        fpixsym = psf->loc + w*h;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++)-*(--fpixsym);
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
    case PSF_GREAT:
      {
       FILE	*listfile;
       char	listfilename[MAXCHARL], str[MAXCHARL], str2[80],
		*ptr;
       double	*list, invcscalex,invcscaley, scalex,scaley, coffx,coffy;
       int	ispoon,spoonsize,size, nlist;

/*---- Reconstruction of the PSF at positions given in an external ASCII file */
      NFPRINTF(OUTPUT, "Reading list of input positions...");
/*--- Create a file name with a ".dat" extension */
      strcpy(listfilename, field->catname);
      if (!(pstr = strrchr(listfilename, '.')))
      pstr = listfilename+strlen(listfilename);
      sprintf(pstr, "%s", ".dat");
      if (!(listfile = fopen(listfilename, "r")))
        error(EXIT_FAILURE, "*Error*: position list not found: ", listfilename);
      if (psf->cx>=0)
        {
        scalex = psf->contextscale[psf->cx];
        invcscalex = 1.0/scalex;
        coffx = psf->contextoffset[psf->cx];
        }
      else
        {
        invcscalex = scalex = 1.0;
        coffx = 0.0;
        }
      if (psf->cy>=0)
        {
        scaley = psf->contextscale[psf->cy];
        invcscaley = 1.0/scaley;
        coffy = psf->contextoffset[psf->cy];
        }
      else
        {
        invcscaley = scaley = 1.0;
        coffy = 0.0;
        }

      list = NULL;
      size = 0;
      for (i=0; fgets(str, MAXCHARL, listfile);)
        {
/*------ Examine current input line (discard empty and comment lines) */
        if (!*str || strchr("#\t\n",*str))
          continue;

        if (!i)
          {
/*-------- Allocate memory for the filtered list */
          ispoon = 1000;
          spoonsize = ispoon*2*sizeof(double);
          QMALLOC(list, double, size = spoonsize);
          }
        else  if (!(i%ispoon))
          QREALLOC(list, double, size += spoonsize);

         if (!(i%1000))
          {
          sprintf(str2, "Reading input list... (%d objects)", i);
          NFPRINTF(OUTPUT, str2);
          }

        list[2*i] = (strtof(str, &ptr)-coffx)*invcscalex;
        list[2*i+1] = (strtof(ptr, NULL)-coffy)*invcscaley;
        i++;
        }
      fclose(listfile);
      nlist = i;
      if (!nlist)
        warning("No valid positions found in ", listfilename);

      if (cubeflag)
        {
        tab->naxis = 3;
        QREALLOC(tab->naxisn, int, tab->naxis);
        tab->naxisn[0] = 48;
        tab->naxisn[1] = 48;
        tab->naxisn[2] = nlist? nlist : 1;
        npix = tab->naxisn[0]*tab->naxisn[1];
        tab->tabsize = tab->bytepix*npix*tab->naxisn[2];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        pix = pix0;
        for (n=0; n<nlist; n++)
          {
          psf_build(psf, &list[2*n]);
          list[2*n] = list[2*n]*scalex + coffx;
          list[2*n+1] = list[2*n+1]*scaley + coffy;
          vignet_resample(psf->loc, psf->size[0], psf->size[1],
		pix, 48, 48,
		(floor(list[2*n]+0.49999) - list[2*n] + 0.5) / psf->pixstep,
		(floor(list[2*n+1]+0.49999) - list[2*n+1] + 0.5) / psf->pixstep,
		1.0/psf->pixstep, 1.0);
          pix += 48*48;
          }
        }
      else
        {
        nw = (int)sqrt((double)nlist);
        nw = nw? ((nw-1)/10+1)*10 : 1;
        nh = (nlist-1)/nw + 1;
        w = 48;
        h = 48;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        step = (nw-1)*w;
        QCALLOC(pix0, float, tab->tabsize);
        QMALLOC(vig0, float, w*h);
        tab->bodybuf = (char *)pix0; 
        for (n=0; n<nlist; n++)
          {
          psf_build(psf, &list[2*n]);
          list[2*n] = list[2*n]*scalex + coffx;
          list[2*n+1] = list[2*n+1]*scaley + coffy;
          vignet_resample(psf->loc, psf->size[0], psf->size[1],
		vig0, 48, 48,
		(floor(list[2*n]+0.49999) - list[2*n] + 0.5) / psf->pixstep,
		(floor(list[2*n+1]+0.49999) - list[2*n+1] + 0.5) / psf->pixstep,
		1.0/psf->pixstep, 1.0);
          vig = vig0;
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(vig++);
          }
        free(vig0);
        }
      free(list);
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
    if (tab->naxis>2)
      {
      addkeywordto_head(tab,"CTYPE3  ","WCS axis type");
      fitswrite(tab->headbuf, "CTYPE3  ", " ",  H_STRING, T_STRING);
      }
    if (tab->naxis>3)
      {
      addkeywordto_head(tab,"CTYPE4  ","WCS axis type");
      fitswrite(tab->headbuf, "CTYPE4  ", " ",  H_STRING, T_STRING);
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
    if (tab->naxis>2)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CRVAL3  ",
	"WCS coordinates of the reference pixel");
      fitswrite(tab->headbuf, "CRVAL3  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    if (tab->naxis>3)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CRVAL4  ",
	"WCS coordinates of the reference pixel");
      fitswrite(tab->headbuf, "CRVAL4  ", &dval1,  H_EXPO, T_DOUBLE);
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
    if (tab->naxis>2)
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
    if (tab->naxis>3)
      {
      dval1 = 0.0;
      addkeywordto_head(tab, "CD1_4   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD1_4   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD2_4   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD2_4   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD3_4   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD3_4   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD4_1   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD4_1   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD4_2   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD4_2   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 0.0;
      addkeywordto_head(tab, "CD4_3   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD4_3   ", &dval1,  H_EXPO, T_DOUBLE);
      dval1 = 1.0;
      addkeywordto_head(tab, "CD4_4   ", "WCS transformation matrix");
      fitswrite(tab->headbuf, "CD4_4   ", &dval1,  H_EXPO, T_DOUBLE);
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
    if (tab->naxis>2)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CDELT3  ","WCS pixel scale");
      fitswrite(tab->headbuf, "CDELT3  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    if (tab->naxis>3)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CDELT4  ","WCS pixel scale");
      fitswrite(tab->headbuf, "CDELT4  ", &dval1,  H_EXPO, T_DOUBLE);
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
    if (tab->naxis>2)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CRPIX3  ",
	"pixel coordinates of the reference pixel");
      fitswrite(tab->headbuf, "CRPIX3  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    if (tab->naxis>3)
      {
      dval1 = 1.0;
      addkeywordto_head(tab,"CRPIX4  ",
	"pixel coordinates of the reference pixel");
      fitswrite(tab->headbuf, "CRPIX4  ", &dval1,  H_EXPO, T_DOUBLE);
      }
    }

  if (fitsread(set->head, "EPOCH   ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    addkeywordto_head(tab, "EPOCH   ", "Epoch of observation");
    fitswrite(tab->headbuf, "EPOCH   ", &dval1,  H_EXPO, T_DOUBLE);
    }
  if (fitsread(set->head, "EQUINOX ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    addkeywordto_head(tab, "EQUINOX ", "Equinox");
    fitswrite(tab->headbuf, "EQUINOX ", &dval1,  H_EXPO, T_DOUBLE);
    }
  if (fitsread(set->head, "RADECSYS", str, H_STRING, T_STRING)==RETURN_OK)
    {
    addkeywordto_head(tab, "RADECSYS", "Coordinate system");
    fitswrite(tab->headbuf, "RADECSYS", str,  H_STRING, T_STRING);
    }
  if (fitsread(set->head, "LONGPOLE", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    addkeywordto_head(tab, "LONGPOLE", "Longitude of pole");
    fitswrite(tab->headbuf, "LONGPOLE", &dval1,  H_EXPO, T_DOUBLE);
    }
  if (fitsread(set->head, "LATPOLE ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    addkeywordto_head(tab, "LATPOLE ", "Latitude of pole");
    fitswrite(tab->headbuf, "LATPOLE ", &dval1,  H_EXPO, T_DOUBLE);
    }
  if (fitsread(set->head, "EPOCH   ", &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
    {
    addkeywordto_head(tab, "EPOCH   ", "Epoch of observation");
    fitswrite(tab->headbuf, "EPOCH   ", &dval1,  H_EXPO, T_DOUBLE);
    }
  if (fitsfind(set->head, "PV?_????") != RETURN_ERROR)
    {
    for (l=0; l<2; l++)
      for (j=0; j<100; j++)
        {
        sprintf(str, "PV%d_%d ", l+1, j);
        if (fitsread(set->head, str, &dval1, H_EXPO, T_DOUBLE)==RETURN_OK)
          {
          addkeywordto_head(tab, str, "Distortion parameter");
          fitswrite(tab->headbuf, str, &dval1,  H_EXPO, T_DOUBLE);
          }
        }
    }

/*-- Add and write important scalars as FITS keywords */
  addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data");
  val = 0.0;
  fitswrite(tab->headbuf, "PSF_SAMP",
	psf->samples_accepted? &psf->pixstep : &val, H_FLOAT, T_FLOAT);

  if (next == 1)
    prim_head(tab);
  fitswrite(tab->headbuf, "XTENSION", "IMAGE   ", H_STRING, T_STRING);

/* save table */
  save_tab(cat, tab);
  free_tab(tab);
  if (ext==next-1)
    free_cat(&cat, 1);

  return;
  }

