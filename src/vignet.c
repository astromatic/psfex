/*
*				vignet.c
*
* Manipulate image rasters.
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
*	Last modified:		20/11/2013
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
#include	"vignet.h"


/****** vignet_resample ******************************************************
PROTO	int	vignet_resample(float *pix1, int w1, int h1,
		float *pix2, int w2, int h2, double dx, double dy, float step2,
		float stepi)
PURPOSE	Scale and shift a small image through sinc interpolation, with
	adjustable spatial wavelength cut-off. Image parts which lie outside
	boundaries are set to 0.

INPUT	Input raster,
	input raster width,
	input raster height,
	output raster,
	output raster width,
	output raster height,
	shift in x,
	shift in y,
	output pixel scale.	
OUTPUT	RETURN_ERROR if the images do not overlap, RETURN_OK otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	13/09/2010
 ***/
int	vignet_resample(float *pix1, int w1, int h1,
		float *pix2, int w2, int h2, double dx, double dy, float step2,
		float stepi)
  {
   static float	*statpix2;
   double	*mask,*maskt, mx1,mx2,my1,my2, xs1,ys1, x1,y1, x,y, dxm,dym,
		val, dstepi, norm;
   float	*pix12, *pixin,*pixin0, *pixout,*pixout0;
   int		i,j,k,n,t, *start,*startt, *nmask,*nmaskt,
		ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
		ix,iy, ix1,iy1, interpw, interph;

  if (stepi <= 0.0)
    stepi = 1.0;
  dstepi = 1.0/stepi;
  mx1 = (double)(w1/2);		/* Im1 center x-coord*/
  mx2 = (double)(w2/2);		/* Im2 center x-coord*/
  xs1 = mx1 + dx - mx2*step2;	/* Im1 start x-coord */

  if ((int)xs1 >= w1)
    return RETURN_ERROR;
  ixs2 = 0;			/* Int part of Im2 start x-coord */
  if (xs1<0.0)
    {
    dix2 = (int)(1-xs1/step2);
/*-- Simply leave here if the images do not overlap in x */
    if (dix2 >= w2)
      return RETURN_ERROR;
    ixs2 += dix2;
    xs1 += dix2*step2;
    }
  nx2 = (int)((w1-1-xs1)/step2+1);/* nb of interpolated Im2 pixels along x */
  if (nx2>(ix2=w2-ixs2))
    nx2 = ix2;
  if (nx2<=0)
    return RETURN_ERROR;
  my1 = (double)(h1/2);		/* Im1 center y-coord */
  my2 = (double)(h2/2);		/* Im2 center y-coord */
  ys1 = my1 + dy - my2*step2;	/* Im1 start y-coord */
  if ((int)ys1 >= h1)
    return RETURN_ERROR;
  iys2 = 0;			/* Int part of Im2 start y-coord */
  if (ys1<0.0)
    {
    diy2 = (int)(1-ys1/step2);
/*-- Simply leave here if the images do not overlap in y */
    if (diy2 >= h2)
      return RETURN_ERROR;
    iys2 += diy2;
    ys1 += diy2*step2;
    }
  ny2 = (int)((h1-1-ys1)/step2+1);/* nb of interpolated Im2 pixels along y */
  if (ny2>(iy2=h2-iys2))
    ny2 = iy2;
  if (ny2<=0)
    return RETURN_ERROR;

/* Set the yrange for the x-resampling with some margin for interpolation */
  iys1a = (int)ys1;		/* Int part of Im1 start y-coord with margin */
  hmh = (int)((INTERPW/2)/dstepi) + 2;	/* Interpolant start */
  interph = 2*hmh;
  hmw = (int)((INTERPW/2)/dstepi) + 2;
  interpw =  2*hmw;
  if (iys1a<0 || ((iys1a -= hmh)< 0))
    iys1a = 0;
  ny1 = (int)(ys1+ny2*step2)+interpw-hmh;	/* Interpolated Im1 y size */
  if (ny1>h1)					/* with margin */
    ny1 = h1;
/* Express everything relative to the effective Im1 start (with margin) */
  ny1 -= iys1a;
  ys1 -= (double)iys1a;

/* Allocate interpolant stuff for the x direction */
  QMALLOC(mask, double, nx2*interpw);	/* Interpolation masks */
  QMALLOC(nmask, int, nx2);		/* Interpolation mask sizes */
  QMALLOC(start, int, nx2);		/* Int part of Im1 conv starts */
/* Compute the local interpolant and data starting points in x */
  x1 = xs1;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=nx2; j--; x1+=step2)
    {
    ix = (ix1=(int)x1) - hmw;
    dxm = (ix1 - x1 - hmw)*dstepi;/* starting point in the interp. func */
    if (ix < 0)
      {
      n = interpw+ix;
      dxm -= (double)ix*dstepi;
      ix = 0;
      }
    else
      n = interpw;
    if (n>(t=w1-ix))
      n=t;
    *(startt++) = ix;
    *(nmaskt++) = n;
    norm = 0.0;
    for (x=dxm, i=n; i--; x+=dstepi)
      norm +=( *(maskt++) = INTERPF(x));
    norm = norm>0.0? 1.0/norm : dstepi;
    maskt -= n;
    for (i=n; i--;)
      *(maskt++) *= norm;
    }

  QCALLOC(pix12, float, nx2*ny1);	/* Intermediary frame-buffer */

/* Make the interpolation in x (this includes transposition) */
  pixin0 = pix1+iys1a*w1;
  pixout0 = pix12;
  for (k=ny1; k--; pixin0+=w1, pixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    pixout = pixout0;
    for (j=nx2; j--; pixout+=ny1)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)*(double)*(pixin++);
      *pixout = (float)val;
      }
    }

/* Reallocate interpolant stuff for the y direction */
  QREALLOC(mask, double, ny2*interph);	/* Interpolation masks */
  QREALLOC(nmask, int, ny2);		/* Interpolation mask sizes */
  QREALLOC(start, int, ny2);		/* Int part of Im1 conv starts */

/* Compute the local interpolant and data starting points in y */
  y1 = ys1;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=ny2; j--; y1+=step2)
    {
    iy = (iy1=(int)y1) - hmh;
    dym = (iy1 - y1 - hmh)*dstepi;/* starting point in the interp. func */
    if (iy < 0)
      {
      n = interph+iy;
      dym -= (double)iy*dstepi;
      iy = 0;
      }
    else
      n = interph;
    if (n>(t=ny1-iy))
      n=t;
    *(startt++) = iy;
    *(nmaskt++) = n;
    norm = 0.0;
    for (y=dym, i=n; i--; y+=dstepi)
      norm += (*(maskt++) = INTERPF(y));
    norm = norm>0.0? 1.0/norm : dstepi;
    maskt -= n;
    for (i=n; i--;)
      *(maskt++) *= norm;
    }

/* Initialize destination buffer to zero if pix2 != NULL */
  if (!pix2)
    pix2 = statpix2;
  else
    {
    memset(pix2, 0, (size_t)(w2*h2)*sizeof(float));
    statpix2 = pix2;
    }

/* Make the interpolation in y  and transpose once again */
  pixin0 = pix12;
  pixout0 = pix2+ixs2+iys2*w2;
  for (k=nx2; k--; pixin0+=ny1, pixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    pixout = pixout0;
    for (j=ny2; j--; pixout+=w2)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)*(double)*(pixin++);
      *pixout = (float)val;
      }
    }

/* Free memory */
  free(pix12);
  free(mask);
  free(nmask);
  free(start);

  return RETURN_OK;
  }


/******************************** vignet_copy ********************************/
/*
Copy a small part of the image. Image parts which lie outside boundaries are
set to 0.
*/
int     vignet_copy(float *pix1, int w1, int h1,
		float *pix2, int w2, int h2, int idx, int idy, vigopenum vigop)
  {
   int          x,y, xmin,ymin, nx,ny, off1,off2;

  if (vigop==VIGNET_CPY)
/*-- First put the pix2 background to zero */
    memset(pix2, 0, (size_t)(w2*h2)*sizeof(float));

/* Set the image boundaries */
  ymin = h2/2+idy-h1/2;
  if ((ny=h2-ymin)>h1)
    ny = h1;
  else if (ny<=0)
    return RETURN_ERROR;
  if (ymin<0)
    {
    pix1 -= ymin*w1;
    ny += ymin;
    }
  else
    pix2 += ymin*w2;

  xmin = w2/2+idx-w1/2;
  if ((nx=w2-xmin)>w1)
    nx = w1;
  else if (nx<=0)
    return RETURN_ERROR;
  if (xmin<0)
    {
    pix1 -= xmin;
    nx += xmin;
    }
  else
    pix2 += xmin;

/* Offsets */
  off1 = w1-nx;
  off2 = w2-nx;
/* Copy the right pixels to the destination */
  switch(vigop)
    {
    case VIGNET_CPY:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) = *(pix1++);
      break;

    case VIGNET_ADD:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) += *(pix1++);
      break;

    case VIGNET_SUB:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) -= *(pix1++);
      break;

    case VIGNET_MUL:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--;)
          *(pix2++) *= *(pix1++);
      break;

    case VIGNET_DIV:
      for (y=ny; y--; pix1+=off1, pix2+=off2)
        for (x=nx; x--; pix2++)
          if (*pix1)
            *pix2 /= *(pix1++);
          else
            *pix2 = (*pix2>0.0)? BIG:-BIG;

    default:
      error(EXIT_FAILURE, "*Internal Error*: unknown operation in ",
			"vignet_copy()");
      
      break;
    }

  return RETURN_OK;
  }


 /**************************** vignet_aperflux******************************/
/*
Compute the total flux within a circular aperture.
*/
float	vignet_aperflux(float *ima, float *var, int w, int h,
			float dxc, float dyc, float aper,
			float gain, float backnoise, float *fluxvar)

  {
   float		*imat,*vart,
			r2, raper,raper2, rintlim,rintlim2,rextlim2, mx,my,
			dx,dx1,dy,dy2, pix,pvar,invbacknoise2, invgain, area,
			offsetx,offsety,scalex,scaley,scale2, locarea, vthresh;
   double		tv, sigtv;
   int			x,y, x2,y2, xmin,xmax,ymin,ymax, sx,sy,
			fymin,fymax, pflag,corrflag;
   long			pos;

  pvar = backnoise*backnoise;
  invbacknoise2 = pvar>0.0? 1.0/pvar : 0.0;
  invgain = gain>0.0? 1.0/gain : 0.0;
/* Integration radius */
  raper = aper/2.0;
  raper2 = raper*raper;
/* Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
  rintlim = raper - 0.75;
  rintlim2 = (rintlim>0.0)? rintlim*rintlim: 0.0;
/* External radius of the oversampled annulus (>r+sqrt(2)/2) */
  rextlim2 = (raper + 0.75)*(raper + 0.75);
  tv = sigtv = area = 0.0;
  scaley = scalex = 1.0/APER_OVERSAMP;
  scale2 = scalex*scaley;
  offsetx = 0.5*(scalex-1.0);
  offsety = 0.5*(scaley-1.0);
  vthresh = BIG/2.0;
  mx = dxc + (float)(w/2);
  my = dyc + (float)(h/2);

  xmin = (int)(mx-raper+0.499999);
  xmax = (int)(mx+raper+1.499999);
  ymin = (int)(my-raper+0.499999);
  ymax = (int)(my+raper+1.499999);
  if (xmin < 0 || xmax > w || ymin < 0 || ymax > h)
    return 0.0;

  for (y=ymin; y<ymax; y++)
    {
    imat = ima + (pos = y*w + xmin);
    vart = var + pos;
    for (x=xmin; x<xmax; x++, imat++, vart++)
      {
      dx = x - mx;
      dy = y - my;
      if ((r2=dx*dx+dy*dy) < rextlim2)
        {
        if (r2> rintlim2)
          {
          dx += offsetx;
          dy += offsety;
          locarea = 0.0;
          for (sy=APER_OVERSAMP; sy--; dy+=scaley)
            {
            dx1 = dx;
            dy2 = dy*dy;
            for (sx=APER_OVERSAMP; sx--; dx1+=scalex)
              if (dx1*dx1+dy2<raper2)
                locarea += scale2;
            }
          }
        else
          locarea = 1.0;
        area += locarea;
/*------ Here begin tests for pixel and/or weight overflows. Things are a */
/*------ bit intricated to have it running as fast as possible in the most */
/*------ common cases */
        if ((pix=*imat)<=-BIG || (var && (pvar=*vart)>=vthresh))
          {
          if ((x2=(int)(2*mx+0.49999-x))>=0 && x2<w
		&& (y2=(int)(2*my+0.49999-y))>=0 && y2<h
		&& (pix=*(imat + (pos = y2*w + x2)))>-BIG)
            {
            if (var)
              {
              pvar = *(var + pos);
              if (pvar>=vthresh)
                pix = pvar = 0.0;
              }
            }
          else
            {
            pix = 0.0;
            if (var)
              pvar = 0.0;
            }
          }
        tv += locarea*pix;
        sigtv += locarea*pvar;
        if (pix>0.0 && gain>0.0)
          sigtv += pix*invgain*pvar*invbacknoise2;
        }
      }
    }

  if (tv>0.0)
    sigtv += tv*invgain;

  if (fluxvar)
    *fluxvar = sqrt(sigtv);

  return tv;
  }



