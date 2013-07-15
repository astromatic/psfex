/*
*				utils.c
*
* Miscellaneous functions.
*
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
*	Last modified:		14/07/2013, Bastille Day
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<assert.h>
#include	<ctype.h>
#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include        "define.h"
#include        "globals.h"
#include        "prefs.h"

/*****************************************************************************/
/*
 * Copied here from sample.c as the rest of that file's all about fits I/O
 */
#include "sample.h"

/****** malloc_samples *******************************************************
PROTO   void malloc_samples(setstruct *set, int nsample)
PURPOSE Allocate memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
void	malloc_samples(setstruct *set, int nsample)

  {
   samplestruct	*sample;
   int		n;

  QMALLOC(set->sample, samplestruct *, nsample);
  set->nsamplemax = nsample;
  if (set != set->samples_owner)
    return;

  QMALLOC(set->s_sample, samplestruct, nsample);
  for (n=0; n!=nsample; ++n)
    {
    sample = set->sample[n] = &set->s_sample[n];
    QMALLOC(sample->vig, float, set->nvig);
    QMALLOC(sample->vigresi, float, set->nvig);
    QMALLOC(sample->vigweight, float, set->nvig);
    QMALLOC(sample->vigchi, float, set->nvig);
    if (set->ncontext)
      QMALLOC(sample->context, double, set->ncontext);
    }

  return;
  }


/****** realloc_samples ******************************************************
PROTO   void realloc_samples(setstruct *set, int nsample)
PURPOSE Re-allocate memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
void	realloc_samples(setstruct *set, int nsample)

  {
   samplestruct	*sample;
   int		n;

  if (set != set->samples_owner)
    {
    assert(nsample <= set->nsamplemax);
    return;
    }

/* If we want to reallocate 0 samples, better free the whole thing! */
  if (!nsample)
    free_samples(set);

/* Two cases: either more samples are required, or the opposite! */
  if (nsample>set->nsamplemax)
    {
    QREALLOC(set->s_sample, samplestruct, nsample);
    free(set->sample);
    QMALLOC(set->sample, samplestruct *, nsample);
    for (n = set->nsamplemax; n<nsample; n++)
      {
      sample = &set->s_sample[n];
      QMALLOC(sample->vig, float, set->nvig);
      QMALLOC(sample->vigresi, float, set->nvig);
      QMALLOC(sample->vigchi, float, set->nvig);
      QMALLOC(sample->vigweight, float, set->nvig);
      if (set->ncontext)
        QMALLOC(sample->context, double, set->ncontext);
      }
    }
  else if (nsample<set->nsamplemax)
    {
    for (n = nsample; n<set->nsamplemax; n++)
      {
      sample = &set->s_sample[n];
      free(sample->vig);
      free(sample->vigresi);
      free(sample->vigchi);
      free(sample->vigweight);
      if (set->ncontext)
        free(sample->context);
      }
    QREALLOC(set->s_sample, samplestruct, nsample);
    free(set->sample);
    QMALLOC(set->sample, samplestruct *, nsample);
    }

  set->nsamplemax = nsample;

  for (n = 0; n<nsample; n++)
     set->sample[n] = &set->s_sample[n];

  return;
  }


/****** free_samples *********************************************************
PROTO   void free_samples(setstruct *set, int nsample)
PURPOSE free memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
void	free_samples(setstruct *set)

  {
   samplestruct	*sample;
   int		n;

   if (set->samples_owner == set)
     {
     for (n = 0; n<set->nsamplemax; ++n )
       {
       sample = set->sample[n];
       free(sample->vig);
       free(sample->vigresi);
       free(sample->vigweight);
       free(sample->vigchi);
       if (set->ncontext)
	 free(sample->context);
       }
     free(set->s_sample);
     }

  free(set->sample);
  set->sample = NULL;
  set->nsample = set->nsamplemax = 0;

  return;
  }


/****** remove_sample ********************************************************
PROTO   samplestruct *remove_sample(setstruct *set, int isample)
PURPOSE Remove an element from a set of samples.
INPUT   set structure pointer,
        sample number.
OUTPUT  The new pointer for the element that replaced the removed one.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 01/03/99
*/
samplestruct	*remove_sample(setstruct *set, int isample)

  {
   static samplestruct	exsample;
   samplestruct		*sample;
   int			nsample;

/* If we want to reallocate 0 samples, better free the whole thing! */
  nsample = set->nsample-1;
  if (nsample>0)
    {
    sample = set->sample[isample];
    exsample = *set->sample[nsample];
    *set->sample[nsample] = *sample;
    *sample = exsample;
    }
   else
     nsample=0;
  realloc_samples(set, nsample);
  set->nsample = nsample;

  return set->sample[isample];
  }


/****** init_set ************************************************************
PROTO   setstruct *init_set()
PURPOSE Allocate and initialize a set structure.
INPUT   -.
OUTPUT  -.
NOTES   See prefs.h.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
setstruct	*init_set(contextstruct *context)

  {
   setstruct	*set;
   int		i;

  QCALLOC(set, setstruct, 1);
  set->nsample = set->nsamplemax = 0;
  set->samples_owner = set;		/* we'll own the samples, so we'll need to free them */
  set->vigdim = 2;
  QMALLOC(set->vigsize, int, set->vigdim);
  set->ncontext = context->ncontext;
  if (set->ncontext)
    {
    QMALLOC(set->contextoffset, double, set->ncontext);
    QMALLOC(set->contextscale, double, set->ncontext);
    QMALLOC(set->contextname, char *, set->ncontext);
    for (i=0; i<set->ncontext; i++)
      QMALLOC(set->contextname[i], char, 80);
    }

  return set;
  }


/****** end_set *************************************************************
PROTO   void end_set(setstruct *set)
PURPOSE free memory allocated by a complete set structure.
INPUT   set structure pointer,
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
void	end_set(setstruct *set)

  {
   int	i;

  free_samples(set);
  free(set->vigsize);
  if (set->ncontext)
    {
    for (i=0; i<set->ncontext; i++)
      free(set->contextname[i]);
    free(set->contextname);
    free(set->contextoffset);
    free(set->contextscale);
    }
  free(set->head);
  free(set);

  return;
  }


/****** make_weights *********************************************************
PROTO   void make_weights(setstruct *set, samplestruct *sample)
PURPOSE Produce a weight-map for each sample vignet.
INPUT   set structure pointer,
        sample structure pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP,Leiden observatory & ESO)
VERSION 13/08/2007
*/
void make_weights(setstruct *set, samplestruct *sample)

  {
   float	*vig, *vigweight,
		backnoise2, gain, noise2, profaccu2, pix;
   int		i;

/* Produce a weight-map */
  profaccu2 = prefs.prof_accuracy*prefs.prof_accuracy;
  gain = sample->gain;
  backnoise2 = sample->backnoise2;
  for (vig=sample->vig, vigweight=sample->vigweight, i=set->nvig; i--;)
    {
    if (*vig <= -BIG)
      *(vig++) = *(vigweight++) = 0.0;
    else
      {
      pix = *(vig++);
      noise2 = backnoise2 + profaccu2*pix*pix;
      if (pix>0.0 && gain>0.0)
        noise2 += pix/gain;
      *(vigweight++) = 1.0/noise2;      
      }
    }

  return;
  }


/****** recenter_sample ******************************************************
PROTO   void recenter_samples(samplestruct sample,
		setstruct *set, float fluxrad)
PURPOSE Recenter sample image using windowed barycenter.
INPUT   sample structure pointer,
	set structure pointer,
	flux radius.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 16/11/2009
*/
void	recenter_sample(samplestruct *sample, setstruct *set, float fluxrad)

  {
   double	tv, dxpos, dypos;
   float	*ima,*imat,*weight,*weightt,
		pix, var, locpix, locarea, sig,twosig2, raper,raper2,
		offsetx,offsety, mx,my, mx2ph,my2ph,
		rintlim,rintlim2,rextlim2, scalex,scaley,scale2,
		dx,dy, dx1,dy2, r2;
   int		i, x,y,x2,y2, w,h, sx,sy, xmin,xmax,ymin,ymax, pos;

  sig = fluxrad*2.0/2.35; /* From half-FWHM to sigma */
  twosig2 = 2.0*sig*sig;

  w = set->vigsize[0];
  h = set->vigsize[1];
/* Integration radius */
  raper = RECENTER_NSIG*sig;
  raper2 = raper*raper;
/* Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
  rintlim = raper - 0.75;
  rintlim2 = (rintlim>0.0)? rintlim*rintlim: 0.0;
/* External radius of the oversampled annulus (>r+sqrt(2)/2) */
  rextlim2 = (raper + 0.75)*(raper + 0.75);
  scaley = scalex = 1.0/RECENTER_OVERSAMP;
  scale2 = scalex*scaley;
  offsetx = 0.5*(scalex-1.0);
  offsety = 0.5*(scaley-1.0);
/* Use isophotal centroid as a first guess */
  mx = sample->dx + (float)(w/2);
  my = sample->dy + (float)(h/2);

  for (i=0; i<RECENTER_NITERMAX; i++)
    {
    xmin = (int)(mx-raper+0.499999);
    xmax = (int)(mx+raper+1.499999);
    ymin = (int)(my-raper+0.499999);
    ymax = (int)(my+raper+1.499999);
    mx2ph = mx*2.0 + 0.49999;
    my2ph = my*2.0 + 0.49999;

    if (xmin < 0)
      xmin = 0;
    if (xmax > w)
      xmax = w;
    if (ymin < 0)
      ymin = 0;
    if (ymax > h)
      ymax = h;
    tv = 0.0;
    dxpos = dypos = 0.0;
    ima = sample->vig;
    weight = sample->vigweight;
    for (y=ymin; y<ymax; y++)
      {
      imat= ima + (pos = y*w + xmin);
      weightt = weight + pos;
      dy = y - my;
      for (x=xmin; x<xmax; x++, imat++, weightt++)
        {
        dx = x - mx;
        if ((r2=dx*dx+dy*dy)<rextlim2)
          {
          if (RECENTER_OVERSAMP>1 && r2> rintlim2)
            {
            dx += offsetx;
            dy += offsety;
            locarea = 0.0;
            for (sy=RECENTER_OVERSAMP; sy--; dy+=scaley)
              {
              dx1 = dx;
              dy2 = dy*dy;
              for (sx=RECENTER_OVERSAMP; sx--; dx1+=scalex)
                if (dx1*dx1+dy2<raper2)
                  locarea += scale2;
              }
            }
          else
            locarea = 1.0;
          locarea *= expf(-r2/twosig2);
/*-------- Here begin tests for pixel and/or weight overflows. Things are a */
/*-------- bit intricated to have it running as fast as possible in the most */
/*-------- common cases */
          pix = *imat;
          if ((var=*weightt)==0.0)
            {
            if ((x2=(int)(mx2ph-x))>=0 && x2<w
		&& (y2=(int)(my2ph-y))>=0 && y2<h
		&& (var=*(weight + (pos = y2*w + x2)))>0.0)
              {
              var = 1.0/var;
              pix = *(ima+pos);
              }
            else
              pix = var = 0.0;
            }
          dx = x - mx;
          dy = y - my;
          locpix = locarea*pix;
          tv += locpix;
          dxpos += locpix*dx;
          dypos += locpix*dy;
          }
        }
      }

    if (tv>0.0)
      {
      mx += (dxpos /= tv)*RECENTER_GRADFAC;
      my += (dypos /= tv)*RECENTER_GRADFAC;
      }
    else
      break;

/*-- Stop here if position does not change */
    if (dxpos*dxpos+dypos*dypos < RECENTER_STEPMIN*RECENTER_STEPMIN)
      break;
    }

  sample->dx = mx - (float)(w/2);
  sample->dy = my - (float)(h/2);

  return;
  }
