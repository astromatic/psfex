/*
*				cplot.c
*
* Produce check plots.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2008-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		16/04/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<dlfcn.h>
#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	PLPLOT_H
#include	PLPLOTP_H

#include	"define.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"cplot.h"
#include	"fitswcs.h"
#include	"prefs.h"
#include	"field.h"
#include	"psf.h"

struct {cplotdevenum device; char *devname; char *extension;}
		cplot_device[] = {{CPLOT_NULL, "null", ""},
		{CPLOT_XWIN, "xwin",":0"},
		{CPLOT_TK, "tk",""},
		{CPLOT_XTERM, "xterm",""},
		{CPLOT_PLMETA, "plmeta",".plm"},
		{CPLOT_PS, "ps", ".ps"},
		{CPLOT_PSC, "psc", ".ps"},
		{CPLOT_XFIG, "xfig", ".fig"},
		{CPLOT_LJIIP, "ljiip", ".lj"},
		{CPLOT_LJHPGL, "lj_hpgl",".hpg"},
		{CPLOT_IMP, "imp",".imp"},
		{CPLOT_PBM, "pbm",".pbm"},
		{CPLOT_PNG, "png",".png"},
		{CPLOT_JPEG, "jpeg",".jpg"},
		{CPLOT_PSTEX, "pstex", ".ps"},
                {CPLOT_AQT, "aqt", ""},
                {CPLOT_PDF, "pdf", ".pdf"},
                {CPLOT_SVG, "svg", ".svg"},
		{CPLOT_NULL,"",""}};

int	plotnum[CPLOT_NTYPES];
int	plotdev[CPLOT_NTYPES];
char	plotfilename[MAXCHAR];
int	plotaaflag;

/****** cplot_check ***********************************************************
PROTO	int cplot_check(cplotenum cplottype)
PURPOSE	Check that the specified check-plot type has been requested by user,
	by returning its index in CPLOT_TYPE keyword list, or RETURN_ERROR
	(-1) otherwise.
INPUT	Check-plot type.
OUTPUT	Index in CPLOT_TYPE keyword list, or RETURN_ERROR (-1) otherwise.
[5~[5~NOTES	Uses the global preferences.
5~AUTHOR	E. Bertin (IAP)
VERSION	18/10/2002
 ***/
int	 cplot_check(cplotenum cplottype)

  {
   int		i;

  for (i=0; i<prefs.ncplot_type; i++)
    if (cplottype == prefs.cplot_type[i])
      return i;

  return RETURN_ERROR;
  }


/****** cplot_init ***********************************************************
PROTO	int cplot_init(char *name, int nx, int ny, cplotenum cplottype)
PURPOSE	Initialize a check plot.
INPUT	Input name,
	Number of plots along the x axis,
	number of plots along the y axis,
	plot type,
	device number.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_init(char *name, int nx, int ny, cplotenum cplottype)
  {
   char		**devmenu, **devname;
   char		str[MAXCHAR],
		*pstr;
   int		j, num, cval, dev, argc, ndev, gdflag;

/* Check that plot was requested */
  cval = cplot_check(cplottype);
/* return otherwise */
  if (cval == RETURN_ERROR)
    return RETURN_ERROR;
  dev = plotdev[cval]++;
/* Run convert to antialias the check-plot image (quick and dirty) */
  if (prefs.cplot_antialiasflag && dev>0
	&& (prefs.cplot_device[dev-1] == CPLOT_PNG
	|| prefs.cplot_device[dev-1] == CPLOT_JPEG))
    {
    sprintf(str, "convert -geometry \"%dx%d\" -antialias %s %s",
	prefs.cplot_res[0]? prefs.cplot_res[0] : CPLOT_DEFRESX,
	prefs.cplot_res[1]? prefs.cplot_res[1] : CPLOT_DEFRESY,
	plotfilename, plotfilename);
    system(str);
    }

  if (dev>=prefs.ncplot_device)
    return RETURN_ERROR;
  num = plotnum[cval]+1;

/* Provide the right extension to the output filename */
  strcpy(str, name);
  if (!(pstr = strrchr(str, '.')))
    pstr = str+strlen(str);
  *pstr = '\0';
  for (j=0; *(cplot_device[j].devname); j++)
    if (prefs.cplot_device[dev]==cplot_device[j].device)
      break;

  if (cplot_device[j].device != CPLOT_XWIN)
    {
    sprintf(plotfilename, "%s_%s%s", prefs.cplot_name[cval], str,
			cplot_device[j].extension);
    plsfnam(plotfilename);		/* file name */
    }
  plssub(nx, ny);
  plsdev(cplot_device[j].devname);

  plotaaflag = 0;
  if (cplot_device[j].device == CPLOT_PNG
	|| cplot_device[j].device == CPLOT_JPEG)
    {
/*-- Hack to check whether the PNG/JPEG/GIF driver is based on GD or not */
    ndev = 100;	/* max number of plplot devices */
    QMALLOC(devmenu, char *, ndev);
    QMALLOC(devname, char *, ndev);
    plgFileDevs((const char ***)&devmenu, (const char ***)&devname, &ndev);
    gdflag = 0;
    for (dev=0; dev<ndev; dev++)
      if (!strncmp(devmenu[dev],"PNG", 3))
        {
        gdflag = 1;
        break;
        }
    free(devmenu);
    free(devname);

/*-- gd driver is 8 bits by default */ 
    if (gdflag)
      plsetopt("-drvopt","24bit");

/*-- Set custom resolutions */
    if (gdflag && prefs.cplot_antialiasflag)
      {
/*---- Oversample for antialiasing (only GD driver as others are already) */
      sprintf(str, "%dx%d",
	(prefs.cplot_res[0]? prefs.cplot_res[0] : CPLOT_DEFRESX)*CPLOT_AAFAC,
	(prefs.cplot_res[1]? prefs.cplot_res[1] : CPLOT_DEFRESY)*CPLOT_AAFAC);
      plsetopt("-geometry", str);
      plotaaflag = 1;
      }
    else if (prefs.cplot_res[0])
      {
/*---- No oversampling */
      sprintf(str, "%dx%d", prefs.cplot_res[0], prefs.cplot_res[1]);
      plsetopt("-geometry", str);
      }
    }
  else
    {
/*-- Small hack to reset driver options */
    argc = 0;
    plparseopts(&argc, NULL, PL_PARSE_NOPROGRAM);
    }

  plfontld(1);
  plscolbg(255,255,255);	/* Force the background colour to white */
  plscol0(15, 0,0,0);		/* Force the foreground colour to black */
  plscol0(3, 0,200,0);		/* Force the green colour to darken */
  plscol0(7, 128,128,128);	/* Force the midground colour to grey */
  plscol0(8, 64,0,0);		/* Force the brown colour to darken */
  plinit();

  return RETURN_OK;
  }


/****** cplot_end ************************************************************
PROTO	int cplot_end(void)
PURPOSE	Terminate a CPlot (and save additional plots if required).
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2004
 ***/
int	cplot_end(cplotenum cplottype)
  {
   int		cval;

/* Check that plot was requested */
  cval = cplot_check(cplottype);
/* return otherwise */
  if (cval == RETURN_ERROR)
    return RETURN_ERROR;

  plotdev[cval] = 0;
  plotnum[cval]++;

  return RETURN_OK;
  }


/****** cplot_degtosexal ******************************************************
PROTO	char	*cplot_degtosexal(char *str, double alpha, double step)
PURPOSE	Convert degrees to hh mm ss.xx coordinates in PLPLOT string format.
INPUT	Pointer to the output character string,
	alpha coordinate in degrees,
	step (precision) in degrees.
OUTPUT	Pointer to the output character string.
NOTES	At least 30 bytes must have been allocated for str.
AUTHOR	E. Bertin (IAP)
VERSION	20/07/2004
 ***/
char	*cplot_degtosexal(char *str, double alpha, double step)
  {
   int		hh, mm;
   double	dh, dm, ss;

  if (alpha<0.0 || alpha >360.0)
    alpha = fmod(alpha+360.0, 360.0);
  alpha += 1e-8;
  hh = (int)(dh = alpha/15.0);
  mm = (int)(dm = 60.0*(dh - hh));
  ss = 60.0*(dm - mm);
  if (step*DEG<0.999*15.0*ARCSEC)
    sprintf(str,"%02d#uh#d%02d#um#d%05.2f#us#d", hh, mm, ss);
  else if (step*DEG<0.999*15.0*ARCMIN)
    sprintf(str,"%02d#uh#d%02d#um#d%02.0f#us#d", hh, mm, ss);
  else if (step*DEG<0.999*15.0*DEG)
    sprintf(str,"%02d#uh#d%02d#um#d", hh, mm);
  else
    sprintf(str,"%02d#uh#d", hh);

  return str;
  }


/****** cplot_degtosexde *****************************************************
PROTO	char	*cplot_degtosexde(char *str, double delta, double step)
PURPOSE	Convert degrees to dd mm ss.xx coordinates in PLPLOT string format.
INPUT	Pointer to the output character string,
	delta coordinate in degrees,
	step (precision) in degrees.
OUTPUT	Pointer to the output character string.
NOTES	At least 30 bytes must have been allocated for str.
AUTHOR	E. Bertin (IAP)
VERSION	26/10/2005
 ***/
char	*cplot_degtosexde(char *str, double delta, double step)
  {
   char		sign;
   double	ddm, ds;
   int		dd, dm;

  sign = delta<1e-8?'-':'+';
  if (delta<1e-8)
    delta = -delta;
  delta += 1e-8;
  if (delta<-90.0)
    delta = -90.0;
  else if (delta>90.0)
    delta = 90.0;
  dd = (int)delta;
  dm = (int)(ddm = (60.0*(delta - dd)));
  ds = 60.0*(ddm - dm);
  if (step*DEG<0.999*ARCSEC)
    sprintf(str,"%c%02d#(2218)%02d#(2216)%05.2f#(2216)#(2216)", sign,dd,dm,ds);
  else if (step*DEG<0.999*ARCMIN)
    sprintf(str,"%c%02d#(2218)%02d#(2216)%2.0f#(2216)#(2216)", sign,dd,dm,ds);
  else if (step*DEG<0.999*DEG)
    sprintf(str,"%c%02d#(2218)%02d#(2216)", sign,dd,dm);
  else
    sprintf(str,"%c%02d#(2218)", sign,dd);

  return str;
  }


/****** cplot_drawbounds *****************************************************
PROTO	int cplot_drawbounds(wcsstruct *wcsin, wcsstruct *wcsout)
PURPOSE	Draw the projected image boundaries in a given projection.
INPUT	Pointer to the image WCS structure,
	pointer to the plot WCS structure.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	29/10/2008
 ***/
int	cplot_drawbounds(wcsstruct *wcsin, wcsstruct *wcsout)
  {
   PLFLT	x[5],y[5];
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS];
   int		i, lngin,latin, lngout,latout;

/* Initialize the input coordinates to an "average" value */
  for (i=0; i<wcsin->naxis; i++)
    rawpos2[i] = wcsin->naxisn[i]/2.0;

  lngin = wcsin->lng;
  latin = wcsin->lat;
  if (lngin==-1)
    lngin = 0;
  if (latin==-1)
    latin = 1;
  lngout = wcsout->lng;
  latout = wcsout->lat;
  if (lngout==-1)
    lngout = 0;
  if (latout==-1)
    latout = 1;

/* 1st corner */
  rawpos2[lngin] = 0.0;
  rawpos2[latin] = 0.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[lngout] = wcspos2[lngin];
  wcspos[latout] = wcspos2[latin];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[4] = x[0] = rawpos[lngout];
  y[4] = y[0] = rawpos[latout];
/* 2nd corner */
  rawpos2[lngin] = wcsin->naxisn[lngin]-1.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[lngout] = wcspos2[lngin];
  wcspos[latout] = wcspos2[latin];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[1] = rawpos[lngout];
  y[1] = rawpos[latout];
/* 3rd corner */
  rawpos2[latin] = wcsin->naxisn[latin]-1.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[lngout] = wcspos2[lngin];
  wcspos[latout] = wcspos2[latin];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[2] = rawpos[lngout];
  y[2] = rawpos[latout];
/* Last corner */
  rawpos2[lngin] = 0.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[lngout] = wcspos2[lngin];
  wcspos[latout] = wcspos2[latin];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[3] = rawpos[lngout];
  y[3] = rawpos[latout];
/* Draw */
  plline(5, x,y);

  return RETURN_OK;
  }


/****** cplot_drawloccoordgrid ************************************************
PROTO	int cplot_drawloccoordgrid(wcsstruct *wcs, double xmin, double xmax,
				double ymin, double ymax)
PURPOSE	Draw an atlas-like grid of angular coordinates with respect to
	projection center.
INPUT	Pointer to the WCS projection structure,
	left projected coordinate limit,
	right projected coordinate limit,
	bottom projected coordinate limit,
	top projected coordinate limit.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	26/10/2005
 ***/
int	cplot_drawloccoordgrid(wcsstruct *wcs, double xmin, double xmax,
				double ymin, double ymax)
  {
  wcsstruct	*wcs2;
  char		str[32];
  PLFLT		*x,*y;
  double	rawpos[NAXIS], wcspos[NAXIS],
		dx,dy, alphabeg,deltabeg, alphaend,deltaend,
		alphastep,deltastep, xd,yd,xdo,ydo, xm,ym, xmd,ymd, xmu,ymu;
  int		i, lng,lat;

  lng = wcs->lng;
  lat = wcs->lat;
/* Exit if we are not dealing with angular coordinates */
  if (lng==lat)
    return RETURN_ERROR;

/* A small trick to compute the min and max longitudes and latitudes */
  wcs2 = copy_wcs(wcs);
  dx = xmax - xmin;
  dy = ymax - ymin;
  wcs2->naxisn[lng] = (int)(dx+1.0);
  wcs2->naxisn[lat] = (int)(dy+1.0);
  wcs2->crpix[lng] -= xmin;
  wcs2->crpix[lat] -= ymin;
  wcs2->crval[lng] = 0.0;
  wcs2->crval[lat] = 0.0;
  init_wcs(wcs2);
  range_wcs(wcs2);
  alphabeg = fmod(wcs2->wcsmin[lng], 360.0);
  alphaend = fmod(wcs2->wcsmax[lng], 360.0);
  while (alphaend<alphabeg)
    alphaend += 360.0;
  deltabeg = wcs2->wcsmin[lat];
  deltaend = wcs2->wcsmax[lat];

  alphastep = alphaend - alphabeg;
/* Quantize at the 15 degrees level */
  if (alphastep*DEG>45.0*DEG)
     alphastep = 15.0*DEG/DEG;
/* Quantize at the 5 degrees level */
  else if (alphastep*DEG>15.0*DEG)
    alphastep = 5.0*DEG/DEG;
/* Quantize at the 1 degree level */
  else if (alphastep*DEG>5.0*DEG)
    alphastep = 1.0*DEG/DEG;
/* Quantize at the 20 arcmin level */
  else if (alphastep*DEG>1.0*DEG)
    alphastep = 20.0*ARCMIN/DEG;
/* Quantize at the 5 arcmin level */
  else if (alphastep*DEG>20.0*ARCMIN)
    alphastep = 5.0*ARCMIN/DEG;
/* Quantize at the 1 arcmin level */
  else if (alphastep*DEG>5.0*ARCMIN)
    alphastep = 1.0*ARCMIN/DEG;
/* Quantize at the 20 arcsec level */
  else if (alphastep*DEG>1.0*ARCMIN)
    alphastep = 20.0*ARCSEC/DEG;
/* Quantize at the 5 arcsec level */
  else if (alphastep*DEG>20.0*ARCSEC)
    alphastep = 5.0*ARCSEC/DEG;
/* Quantize at the 1 arcsec level */
  else if (alphastep*DEG>5.0*ARCSEC)
    alphastep = 1.0*ARCSEC/DEG;
/* Quantize at the 0.25 arcsec level */
  else if (alphastep*DEG>1.0*ARCSEC)
    alphastep = 0.25*ARCSEC/DEG;
  else
    alphastep = 0.1*ARCSEC/DEG;
  alphabeg -= fmod(alphabeg, alphastep) + alphastep;
  alphabeg = fmod (alphabeg, 360.0);
  alphaend += alphastep;
  while (alphaend<alphabeg)
    alphaend += 360.0;

  deltastep = deltaend - deltabeg;
/* Quantize at the 15 degrees level */
  if (deltastep*DEG>45.0*DEG)
    deltastep = 15.0*DEG/DEG;
/* Quantize at the 5 degrees level */
  else if (deltastep*DEG>15.0*DEG)
    deltastep = 5.0*DEG/DEG;
/* Quantize at the 1 degree level */
  else if (deltastep*DEG>5.0*DEG)
    deltastep = 1.0*DEG/DEG;
/* Quantize at the 20 arcmin level */
  else if (deltastep*DEG>1.0*DEG)
    deltastep = 20.0*ARCMIN/DEG;
/* Quantize at the 5 arcmin level */
  else if (deltastep*DEG>20.0*ARCMIN)
    deltastep = 5.0*ARCMIN/DEG;
/* Quantize at the 1 arcmin level */
  else if (deltastep*DEG>5.0*ARCMIN)
    deltastep = 1.0*ARCMIN/DEG;
/* Quantize at the 20 arcsec level */
  else if (deltastep*DEG>1.0*ARCMIN)
    deltastep = 20.0*ARCSEC/DEG;
/* Quantize at the 5 arcsec level */
  else if (deltastep*DEG>20.0*ARCSEC)
    deltastep = 5.0*ARCSEC/DEG;
/* Quantize at the 1 arcsec level */
  else if (deltastep*DEG>5.0*ARCSEC)
    deltastep = 1.0*ARCSEC/DEG;
/* Quantize at the 0.25 arcsec level */
  else if (deltastep*DEG>1.0*ARCSEC)
    deltastep = 0.25*ARCSEC/DEG;
  else
    deltastep = 0.1*ARCSEC/DEG;
  deltabeg -= fmod(deltabeg, deltastep) + deltastep;
  deltaend += deltastep;
  if (deltabeg<-90.0)
    deltabeg = -90.0;
  if (deltaend>90.0)
    deltabeg = 90.0;

  QMALLOC(x, PLFLT, CPLOT_NPOINTDEF);
  QMALLOC(y, PLFLT, CPLOT_NPOINTDEF);

/* Draw meridians */
  plschr(0.0, 0.33);
  plwid(0);
  pllsty(2);
  xmd = xmu = xdo = -0.5;
  ymd = ymu = ydo = -0.5;
  for (wcspos[0] = alphabeg; wcspos[0]<=alphaend; wcspos[0] += alphastep)
    {
    i=0;
    for (wcspos[1]=deltabeg; wcspos[1]<deltaend && i<CPLOT_NPOINTDEF;
	wcspos[1]+=deltastep, i++)
      {
      wcs_to_raw(wcs2, wcspos, rawpos);
      x[i] = (PLFLT)(xd = rawpos[0]);
      y[i] = (PLFLT)(yd = rawpos[1]);
      if (i>0)
        {
        if ((yd-ymin)*(ydo-ymin) < 0.0 
		&& (xm = (xd-xmin)/dx) > 0.0 && xm < 1.0
		&& fabs(xm-xmd) > 0.1)
          {
          plmtex("b", 2.0, (PLFLT)xm, 0.5, cplot_degtosexde(str,wcspos[0],
		alphastep));
          xmd = xm;
          }
        if ((yd-ymax)*(ydo-ymax) < 0.0
		&& (xm = (xd-xmin)/dx) > 0.0 && xm < 1.0
		&& fabs(xm-xmu) > 0.1)
          {
          plmtex("t", 1.5, (PLFLT)xm, 0.5, cplot_degtosexde(str,wcspos[0],
		alphastep));
          xmu = xm;
          }
        }
      xdo = xd;
      ydo = yd;
      }
    plline(i,x,y);
    }

/* Draw parallels */
  for (wcspos[1] = deltabeg; wcspos[1]<=deltaend; wcspos[1] += deltastep)
    {
    i=0;
    for (wcspos[0]=alphabeg; wcspos[0]<alphaend && i<CPLOT_NPOINTDEF;
	wcspos[0]+=alphastep, i++)
      {
      wcs_to_raw(wcs2, wcspos, rawpos);
      x[i] = (PLFLT)(xd = rawpos[0]);
      y[i] = (PLFLT)(yd = rawpos[1]);
      if (i>0)
        {
        if ((xd-xmin)*(xdo-xmin) < 0.0
		&& (ym = (yd-ymin)/dy) > 0.0 && ym < 1.0
		&& fabs(ym-ymd) > 0.1)
          {
          plmtex("lv", 1.0, (PLFLT)ym, 1.0, cplot_degtosexde(str,wcspos[1],
		deltastep));
          ymd = ym;
          }
        if ((xd-xmax)*(xdo-xmax) < 0.0
		&& (ym = (yd-ymin)/dy) > 0.0 && ym < 1.0
		&& fabs(ym-ymu) > 0.1)
          {
          plmtex("rv", 1.0, (PLFLT)ym, 0.0, cplot_degtosexde(str,wcspos[1],
		deltastep));
          ymu = ym;
          }
	}
      xdo = xd;
      ydo = yd;
      }
    plline(i,x,y);
    }

  free(x);
  free(y);
  end_wcs(wcs2);

  return RETURN_OK;

  }


/****** distort_map *******************************************************   
PROTO   void distort_map(PLFLT x,PLFLT y,PLFLT *tx,PLFLT *ty, void *pltr_data)
PURPOSE Astrometric mapping function for shade plots in cplot_distort().
INPUT   TBW
OUTPUT  -.
NOTES   see http://plplot.sourceforge.net/docbook-manual/
            plplot-html-5.5.3/contour-plots.html#contour-plots-c
AUTHOR  E. Bertin (IAP)
VERSION 03/11/2009
***/
static void distort_map(PLFLT x,PLFLT y,PLFLT *tx,PLFLT *ty, void *pltr_data)
{
 wcsstruct	**wcs;
 double	 	rawpos[NAXIS], wcspos[NAXIS], wcspos2[NAXIS];
 int		i, lng0,lat0, lng1,lat1, naxis, nsnap;


  wcs = (wcsstruct **)pltr_data;
  if ((naxis=wcs[0]->naxis) > 2)
    for (i=0; i<naxis; i++)
      rawpos[i]= wcs[0]->naxisn[i]/2.0 + 0.5;
  lng0 = wcs[0]->lng;
  lat0 = wcs[0]->lat;
  if (lng0==-1)
    lng0 = 0;
  if (lat0==-1)
    lat0 = 1;
  lng1 = wcs[1]->lng;
  lat1 = wcs[1]->lat;
  if (lng1==-1)
    lng1 = 0;
  if (lat1==-1)
    lat1 = 1;
  nsnap = prefs.context_nsnap > 1? prefs.context_nsnap - 1 : 1;
  rawpos[lng0] = x*wcs[0]->naxisn[lng0]/nsnap + 0.5;
  rawpos[lat0] = y*wcs[0]->naxisn[lat0]/nsnap + 0.5;
  raw_to_wcs(wcs[0], rawpos, wcspos);
  wcspos2[lng1] = wcspos[lng0];
  wcspos2[lat1] = wcspos[lat0];
  wcs_to_raw(wcs[1], wcspos2, rawpos);
  *tx = rawpos[lng1];
  *ty = rawpos[lat1];

  return;
  }


/****** cplot_fwhm ***********************************************************
PROTO	int cplot_fwhm(fieldstruct *field)
PURPOSE	Plot a map of the PSF FWHM in the instrument field.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_fwhm(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   psfstruct	*psf;
   PLFLT	**fwhm,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		afwhm,fwhmmin,fwhmmax, mfwhm,dfwhm;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		i,j, e, n,n2,ncx,ncy,nt, nfwhm, naxis, nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_FWHM) == RETURN_ERROR)
    {
    cplot_end(CPLOT_FWHM);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: FWHM map", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  fwhmmin = BIG;
  fwhmmax = -BIG;

/* First pass through the data to find min and max FWHMs */
  flag = 0;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (!psf->samples_accepted)
      continue;
    flag = 1;
/*-- Compute total number of snapshots */
    for (nt=1, n=psf->poly->ndim; (n--)>0;)
      nt *= psf->nsnap;
    for (n=0; n<nt; n++)
      {
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      if (psf->cx >= 0)
        raw[0] = (psf->moffat[n].context[psf->cx]-psf->contextoffset[psf->cx])
		/ psf->contextscale[psf->cx];
      if (psf->cy >= 0)
        raw[1] = (psf->moffat[n].context[psf->cy]-psf->contextoffset[psf->cy])
		/ psf->contextscale[psf->cy];
      afwhm = sqrt(psf->moffat[n].fwhm_min*psf->moffat[n].fwhm_max
		* wcs_scale(wcs, raw));
      if (afwhm<fwhmmin)
        fwhmmin = afwhm;
      if (afwhm>fwhmmax)
        fwhmmax = afwhm;
      }
    }

/* Lower bound to variability in FWHM is 1e-4 */
  if (!flag)
    fwhmmin = fwhmmax = ARCSEC/DEG;
  if ((mfwhm=(fwhmmin+fwhmmax)/2.0) < 1.0e-10*ARCSEC/DEG)
    mfwhm = 1.0e-10*ARCSEC/DEG;
  if ((dfwhm=fwhmmax-fwhmmin) < 1.0e-4*mfwhm)
    dfwhm = 1.0e-4*mfwhm;
  fwhmmin = mfwhm - dfwhm/2.0;
  fwhmmax = mfwhm + dfwhm/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = fwhmmin + (i-0.5) * dfwhm / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real 2D FWHM mapping */
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (psf->samples_accepted)
      {
      ncx = ncy = nt = 1;
      nsnap2 = psf->nsnap>1? psf->nsnap : 2;
      for (n=0; n<psf->poly->ndim; n++)
        {
        nt *= nsnap2;
        if (psf->cx>=0 && n<psf->cx)
          ncx *= nsnap2;
        if (psf->cy>=0 && n<psf->cy)
          ncy *= nsnap2;
        }
      plAlloc2dGrid(&fwhm, nsnap2, nsnap2);
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      xstep = wcs->naxisn[0] / (nsnap2-1);
      ystep = wcs->naxisn[1] / (nsnap2-1);
      raw[1] = 0.5;
      for (j=0; j<nsnap2; j++)
        {
        raw[0] = 0.5;
        for (i=0; i<nsnap2; i++)
          {
          dval = 0.0;
          nfwhm = 0;
/*-------- We average all PSF FWHMs at a given X and Y set of coordinates */
          for (n=0; n<nt; n++)
            if ((psf->cx<0 || (n/ncx)%nsnap2 == i)
		&& (psf->cy<0 || (n/ncy)%nsnap2 == j))
              {
              n2 = psf->nsnap>1? n : 0;
              dval += sqrt(psf->moffat[n2].fwhm_min*psf->moffat[n2].fwhm_max
			* wcs_scale(wcs, raw));
              nfwhm++;
              }
          fwhm[i][j] = dval / nfwhm ;
          raw[0] += xstep;
          }
        raw[1] += ystep;
        }

      wcsptr[0] = wcs;
      wcsptr[1] = wcsout;
      plshades(fwhm, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
      plFree2dGrid(fwhm, nsnap2, nsnap2);
      }
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&fwhm, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    fwhm[0][j] = fwhm[1][j] = fwhmmin + j * dfwhm/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  if (wcs->lng == wcs->lat)
    {
    plwind(0.0,1.0,fwhmmin,fwhmmax);
    plshades(fwhm, 2, CPLOT_NSHADES, NULL, 0.0, 1.0,
	   fwhmmin,fwhmmax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
    }
  else
    {
    plwind(0.0,1.0,fwhmmin*DEG/ARCSEC,fwhmmax*DEG/ARCSEC);
    plshades(fwhm, 2, CPLOT_NSHADES, NULL, 0.0, 1.0,
	   fwhmmin*DEG/ARCSEC,fwhmmax*DEG/ARCSEC, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
    }
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  sprintf(str, "%s", (wcs->lng == wcs->lat)?
	  "pixels"
	: (mfwhm < 0.09*ARCSEC/DEG?
	  "mas" : (mfwhm < ARCMIN/DEG?
		   "arcsec" : (mfwhm < 1.0? "arcmin": "deg"))));
  plschr(0.0, 0.6);
  plmtex("l", 5.0, 0.5, 0.5, str);
  plmtex("b", 2.0, 0.5, 0.5, "PSF FWHM");

/* Draw right colour scale */
  mfwhm /= 100.0;			/* convert to percentage */
  fwhmmin = fwhmmin/mfwhm - 100.0;
  fwhmmax = fwhmmax/mfwhm - 100.0;
  plwind(0.0,1.0,fwhmmin,fwhmmax);
  plschr(0.0, 0.5);
  plbox("", 0.0, 0, "cmstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("r", 5.0, 0.5, 0.5, "%");

  plFree2dGrid(fwhm, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_fwhm(field);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_ellipticity *****************************************************
PROTO	int cplot_ellipticity(fieldstruct *field)
PURPOSE	Plot a map of the PSF ellipticity in the instrument field.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_ellipticity(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   psfstruct	*psf;
   PLFLT	**ellip,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		aellip,ellipmin,ellipmax, mellip,dellip;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		i,j, e, n,n2, ncx,ncy,nt, nellip, naxis, nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_ELLIPTICITY) == RETURN_ERROR)
    {
    cplot_end(CPLOT_ELLIPTICITY);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: ellipticity map", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  ellipmin = BIG;
  ellipmax = -BIG;

/* First pass through the data to find min and max ellipticities */
  flag = 0;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (!psf->samples_accepted)
      continue;
    flag = 1;
/*-- Compute total number of snapshots */
    for (nt=1, n=psf->poly->ndim; (n--)>0;)
      nt *= psf->nsnap;
    for (n=0; n<nt; n++)
      {
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      if (psf->cx >= 0)
        raw[0] = (psf->moffat[n].context[psf->cx]-psf->contextoffset[psf->cx])
		/ psf->contextscale[psf->cx];
      if (psf->cy >= 0)
        raw[1] = (psf->moffat[n].context[psf->cy]-psf->contextoffset[psf->cy])
		/ psf->contextscale[psf->cy];
      aellip = (psf->moffat[n].fwhm_max-psf->moffat[n].fwhm_min)
		/ (psf->moffat[n].fwhm_max + psf->moffat[n].fwhm_min);
      if (aellip<ellipmin)
        ellipmin = aellip;
      if (aellip>ellipmax)
        ellipmax = aellip;
      }
    }

/* Lower bound to variability in ellipticity is 1e-4 */
  if (!flag)
    ellipmin = ellipmax = 0.0;
  if ((mellip=(ellipmin+ellipmax)/2.0) < 1.0e-10)
    mellip = 1.0e-10;
  if ((dellip=ellipmax-ellipmin) < 1.0e-4*mellip)
    dellip = 1.0e-4*mellip;
  ellipmin = mellip - dellip/2.0;
  ellipmax = mellip + dellip/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = ellipmin + (i-0.5) * dellip / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real 2D ellipticity mapping */
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (psf->samples_accepted)
      {
      ncx = ncy = nt = 1;
      nsnap2 = psf->nsnap>1? psf->nsnap : 2;
      for (n=0; n<psf->poly->ndim; n++)
        {
        nt *= nsnap2;
        if (psf->cx>=0 && n<psf->cx)
          ncx *= nsnap2;
        if (psf->cy>=0 && n<psf->cy)
          ncy *= nsnap2;
        }
      plAlloc2dGrid(&ellip, nsnap2, nsnap2);
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      xstep = wcs->naxisn[0] / (nsnap2-1);
      ystep = wcs->naxisn[1] / (nsnap2-1);
      raw[1] = 0.5;
      for (j=0; j<nsnap2; j++)
        {
        raw[0] = 0.5;
        for (i=0; i<nsnap2; i++)
          {
          dval = 0.0;
          nellip = 0;
/*-------- We average all PSF ellips at a given X and Y set of coordinates */
          for (n=0; n<nt; n++)
            if ((psf->cx<0 || (n/ncx)%nsnap2 == i)
		&& (psf->cy<0 || (n/ncy)%nsnap2 == j))
              {
              n2 = psf->nsnap>1? n : 0;
              dval += (psf->moffat[n2].fwhm_max-psf->moffat[n2].fwhm_min)
		/ (psf->moffat[n2].fwhm_max + psf->moffat[n2].fwhm_min);
              nellip++;
              }
          ellip[i][j] = dval / nellip ;
          raw[0] += xstep;
          }
        raw[1] += ystep;
        }

      wcsptr[0] = wcs;
      wcsptr[1] = wcsout;
      plshades(ellip, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
      plFree2dGrid(ellip, nsnap2, nsnap2);
      }
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&ellip, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    ellip[0][j] = ellip[1][j] = ellipmin + j * dellip/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,ellipmin,ellipmax);
  plshades(ellip, 2, CPLOT_NSHADES, NULL, 0.0, 1.0,
	   ellipmin,ellipmax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  sprintf(str, "(a-b)/(a+b)");
  plschr(0.0, 0.5);
  plmtex("l", 5.0, 0.5, 0.5, str);
  plmtex("b", 2.0, 0.5, 0.5, "PSF ellipticity");

/* Draw right colour scale */
  mellip /= 100.0;			/* convert to percentage */
  ellipmin = ellipmin/mellip - 100.0;
  ellipmax = ellipmax/mellip - 100.0;
  plwind(0.0,1.0,ellipmin,ellipmax);
  plschr(0.0, 0.5);
  plbox("", 0.0, 0, "cmstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("r", 5.0, 0.5, 0.5, "%");

  plFree2dGrid(ellip, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_ellipticity(field);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_moffatresi *****************************************************
PROTO	int cplot_moffatresi(fieldstruct *field)
PURPOSE	Plot a map of Moffat fit residuals in the instrument field.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_moffatresi(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   psfstruct	*psf;
   PLFLT	**resi,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		aresi,resimin,resimax, mresi,dresi;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		i,j, e, n,n2,ncx,ncy,nt, nresi, naxis, nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_MOFFATRESI) == RETURN_ERROR)
    {
    cplot_end(CPLOT_MOFFATRESI);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: map of Moffat fit residuals", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  resimin = BIG;
  resimax = -BIG;

/* First pass through the data to find min and max residuals */
  flag = 0;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (!psf->samples_accepted)
      continue;
    flag = 1;
/*-- Compute total number of snapshots */
    for (nt=1, n=psf->poly->ndim; (n--)>0;)
      nt *= psf->nsnap;
    for (n=0; n<nt; n++)
      {
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      if (psf->cx >= 0)
        raw[0] = (psf->pfmoffat[n].context[psf->cx]-psf->contextoffset[psf->cx])
		/ psf->contextscale[psf->cx];
      if (psf->cy >= 0)
        raw[1] = (psf->pfmoffat[n].context[psf->cy]-psf->contextoffset[psf->cy])
		/ psf->contextscale[psf->cy];
      aresi = psf->pfmoffat[n].residuals;
      if (aresi<resimin)
        resimin = aresi;
      if (aresi>resimax)
        resimax = aresi;
      }
    }

/* Lower bound to variability in residuals is 1e-4 */
  if (!flag)
    resimin = resimax = 0.0;
  if ((mresi=(resimin+resimax)/2.0) < 1.0e-10)
    mresi = 1.0e-10;
  if ((dresi=resimax-resimin) < 1.0e-4*mresi)
    dresi = 1.0e-4*mresi;
  resimin = mresi - dresi/2.0;
  resimax = mresi + dresi/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = resimin + (i-0.5) * dresi / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real mapping of residuals */
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (psf->samples_accepted)
      {
      ncx = ncy = nt = 1;
      nsnap2 = psf->nsnap>1? psf->nsnap : 2;
      for (n=0; n<psf->poly->ndim; n++)
        {
        nt *= nsnap2;
        if (psf->cx>=0 && n<psf->cx)
          ncx *= nsnap2;
        if (psf->cy>=0 && n<psf->cy)
          ncy *= nsnap2;
        }
      plAlloc2dGrid(&resi, nsnap2, nsnap2);
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      xstep = wcs->naxisn[0] / (nsnap2-1);
      ystep = wcs->naxisn[1] / (nsnap2-1);
      raw[1] = 0.5;
      for (j=0; j<nsnap2; j++)
        {
        raw[0] = 0.5;
        for (i=0; i<nsnap2; i++)
          {
          dval = 0.0;
          nresi = 0;
/*-------- We average all PSF residuals at a given X and Y set of coordinates */
          for (n=0; n<nt; n++)
            if ((psf->cx<0 || (n/ncx)%nsnap2 == i)
		&& (psf->cy<0 || (n/ncy)%nsnap2 == j))
              {
              n2 = psf->nsnap>1? n : 0;
              dval += psf->pfmoffat[n2].residuals;
              nresi++;
              }
          resi[i][j] = dval / nresi ;
          raw[0] += xstep;
          }
        raw[1] += ystep;
        }

      wcsptr[0] = wcs;
      wcsptr[1] = wcsout;
      plshades(resi, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
      plFree2dGrid(resi, nsnap2, nsnap2);
      }
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&resi, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    resi[0][j] = resi[1][j] = resimin + j * dresi/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,resimin,resimax);
  plshades(resi, 2, CPLOT_NSHADES, NULL, 0.0, 1.0,
	   resimin,resimax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("b", 2.0, 0.5, 0.5, "Fit residuals");

/* Draw right colour scale */
  mresi /= 100.0;			/* convert to percentage */
  resimin = resimin/mresi - 100.0;
  resimax = resimax/mresi - 100.0;
  plwind(0.0,1.0,resimin,resimax);
  plschr(0.0, 0.5);
  plbox("", 0.0, 0, "cmstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("r", 5.0, 0.5, 0.5, "%");
  plwind(0.0,1.0,resimin,resimax);
  plschr(0.0, 0.5);
  plbox("", 0.0, 0, "cmstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("r", 5.0, 0.5, 0.5, "%");

  plFree2dGrid(resi, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_moffatresi(field);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_asymresi ******************************************************
PROTO	int cplot_asymresi(fieldstruct *field)
PURPOSE	Plot a map of asymmetry residuals in the instrument field.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_asymresi(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   psfstruct	*psf;
   PLFLT	**resi,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		aresi,resimin,resimax, mresi,dresi;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		i,j, e, n,n2,ncx,ncy,nt, nresi, naxis, nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_ASYMRESI) == RETURN_ERROR)
    {
    cplot_end(CPLOT_ASYMRESI);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: PSF asymmetry map", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  resimin = BIG;
  resimax = -BIG;

/* First pass through the data to find min and max residuals */
  flag = 0;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (!psf->samples_accepted)
      continue;
    flag = 1;
/*-- Compute total number of snapshots */
    for (nt=1, n=psf->poly->ndim; (n--)>0;)
      nt *= psf->nsnap;
    for (n=0; n<nt; n++)
      {
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      if (psf->cx >= 0)
        raw[0] = (psf->moffat[n].context[psf->cx]-psf->contextoffset[psf->cx])
		/ psf->contextscale[psf->cx];
      if (psf->cy >= 0)
        raw[1] = (psf->moffat[n].context[psf->cy]-psf->contextoffset[psf->cy])
		/ psf->contextscale[psf->cy];
      aresi = psf->moffat[n].symresiduals;
      if (aresi<resimin)
        resimin = aresi;
      if (aresi>resimax)
        resimax = aresi;
      }
    }

/* Lower bound to variability in residuals is 1e-4 */
  if (!flag)
    resimin = resimax = 0.0;
  if ((mresi=(resimin+resimax)/2.0) < 1.0e-10)
    mresi = 1.0e-10;
  if ((dresi=resimax-resimin) < 1.0e-4*mresi)
    dresi = 1.0e-4*mresi;
  resimin = mresi - dresi/2.0;
  resimax = mresi + dresi/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = resimin + (i-0.5) * dresi / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real mapping of residuals */
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    psf = field->psf[e];
    if (psf->samples_accepted)
      {
      ncx = ncy = nt = 1;
      nsnap2 = psf->nsnap>1? psf->nsnap : 2;
      for (n=0; n<psf->poly->ndim; n++)
        {
        nt *= nsnap2;
        if (psf->cx>=0 && n<psf->cx)
          ncx *= nsnap2;
        if (psf->cy>=0 && n<psf->cy)
          ncy *= nsnap2;
        }
      plAlloc2dGrid(&resi, nsnap2, nsnap2);
      for (i=0; i<naxis; i++)
        raw[i] = wcs->naxisn[i]/2.0 + 0.5;
      xstep = wcs->naxisn[0] / (nsnap2-1);
      ystep = wcs->naxisn[1] / (nsnap2-1);
      raw[1] = 0.5;
      for (j=0; j<nsnap2; j++)
        {
        raw[0] = 0.5;
        for (i=0; i<nsnap2; i++)
          {
          dval = 0.0;
          nresi = 0;
/*-------- We average all PSF residuals at a given X and Y set of coordinates */
          for (n=0; n<nt; n++)
            if ((psf->cx<0 || (n/ncx)%nsnap2 == i)
		&& (psf->cy<0 || (n/ncy)%nsnap2 == j))
              {
              n2 = psf->nsnap>1? n : 0;
              dval += psf->moffat[n2].symresiduals;
              nresi++;
              }
          resi[i][j] = dval / nresi ;
          raw[0] += xstep;
          }
        raw[1] += ystep;
        }

      wcsptr[0] = wcs;
      wcsptr[1] = wcsout;
      plshades(resi, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
      plFree2dGrid(resi, nsnap2, nsnap2);
      }
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&resi, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    resi[0][j] = resi[1][j] = resimin + j * dresi/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,resimin,resimax);
  plshades(resi, 2, CPLOT_NSHADES, NULL, 0.0, 1.0,
	   resimin,resimax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("b", 2.0, 0.5, 0.5, "Asymmetry");

/* Draw right colour scale */
  mresi /= 100.0;			/* convert to percentage */
  resimin = resimin/mresi - 100.0;
  resimax = resimax/mresi - 100.0;
  plwind(0.0,1.0,resimin,resimax);
  plschr(0.0, 0.5);
  plbox("", 0.0, 0, "cmstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("r", 5.0, 0.5, 0.5, "%");

  plFree2dGrid(resi, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_asymresi(field);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_counts *********************************************************
PROTO	int cplot_counts(fieldstruct *field)
PURPOSE	Plot an x,y map of the number of sources that has been loaded.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_counts(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   PLFLT	**count,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		cmin,cmax, mc,dc;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		**icount,
		i,j, e, n,nt, naxis, nsnap,nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_COUNTS) == RETURN_ERROR)
    {
    cplot_end(CPLOT_COUNTS);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: source count map", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  cmin = BIG;
  cmax = -BIG;

/* First pass through the data to find min and max number counts */
  flag = 0;
  icount = field->lcount;
  nsnap = prefs.context_nsnap;
  nt = nsnap*nsnap;
  for (e=0; e<field->next; e++)
    for (n=0; n<nt; n++)
      {
      dval = (double)icount[e][n];
      if (dval < cmin)
        cmin = dval;
      if (dval > cmax)
        cmax = dval;
      }

/* Lower bound to variability in counts is 1e-4 */
  if ((mc=(cmin+cmax)/2.0) <  1.0e-10)
    mc = 1.0e-10;
  if ((dc=cmax-cmin) < 1.0e-4*mc)
    dc = 1.0e-4*mc;
  cmin = mc - dc/2.0;
  cmax = mc + dc/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = cmin + (i-0.5) * dc / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real 2D count mapping */
  nsnap2 = nsnap>1? nsnap : 2;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    plAlloc2dGrid(&count, nsnap2, nsnap2);
    for (i=0; i<naxis; i++)
      raw[i] = wcs->naxisn[i]/2.0 + 0.5;
    xstep = wcs->naxisn[0] / (nsnap2-1);
    ystep = wcs->naxisn[1] / (nsnap2-1);
    if (nsnap > 1)
      for (j=0; j<nsnap2; j++)
        for (i=0; i<nsnap2; i++)
          count[i][j] = (PLFLT)icount[e][j*nsnap+i];
    else
      for (j=0; j<nsnap2; j++)
        for (i=0; i<nsnap2; i++)
          count[i][j] = (PLFLT)icount[e][0];
    wcsptr[0] = wcs;
    wcsptr[1] = wcsout;
    plshades(count, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
    plFree2dGrid(count, nsnap2, nsnap2);
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&count, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    count[0][j] = count[1][j] = cmin + j * dc/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,cmin,cmax);
  plshades(count, 2, CPLOT_NSHADES, NULL, 0.0, 1.0, cmin,cmax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  plschr(0.0, 0.6);
  plmtex("l", 5.0, 0.5, 0.5, "number");
  plmtex("b", 2.0, 0.5, 0.5, "source ##");

  plFree2dGrid(count, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_counts(field);		/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_countfrac *******************************************************
PROTO	int cplot_countfrac(fieldstruct *field)
PURPOSE	Plot an x,y map of the fraction of sources that has been accepted.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_countfrac(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   PLFLT	**count,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		cmin,cmax, mc,dc;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		**iacount,**ilcount,
		i,j, e, n,nt, naxis, nsnap,nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_COUNTFRAC) == RETURN_ERROR)
    {
    cplot_end(CPLOT_COUNTFRAC);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: source count fraction map", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  cmin = BIG;
  cmax = -BIG;

/* First pass through the data to find min and max number counts */
  flag = 0;
  ilcount = field->lcount;
  iacount = field->acount;
  nsnap = prefs.context_nsnap;
  nt = nsnap*nsnap;
  for (e=0; e<field->next; e++)
    for (n=0; n<nt; n++)
      {
      dval = ilcount[e][n] > 0? 100.0*iacount[e][n] / ilcount[e][n] : 0.0;
      if (dval < cmin)
        cmin = dval;
      if (dval > cmax)
        cmax = dval;
      }

/* Lower bound to variability in counts fraction is 1e-4 */
  if ((mc=(cmin+cmax)/2.0) <  1.0e-10)
    mc = 1.0e-10;
  if ((dc=cmax-cmin) < 1.0e-4*mc)
    dc = 1.0e-4*mc;
  cmin = mc - dc/2.0;
  cmax = mc + dc/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = cmin + (i-0.5) * dc / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real 2D count mapping */
  nsnap2 = nsnap>1? nsnap : 2;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    plAlloc2dGrid(&count, nsnap2, nsnap2);
    for (i=0; i<naxis; i++)
      raw[i] = wcs->naxisn[i]/2.0 + 0.5;
    xstep = wcs->naxisn[0] / (nsnap2-1);
    ystep = wcs->naxisn[1] / (nsnap2-1);
    for (j=0; j<nsnap2; j++)
      for (i=0; i<nsnap2; i++)
        {
        n = (nsnap>1)? j*nsnap + i : 0;
        count[i][j] = (PLFLT)(ilcount[e][n] > 0?
			100.0*iacount[e][n] / ilcount[e][n] : 0.0);
        }

    wcsptr[0] = wcs;
    wcsptr[1] = wcsout;
    plshades(count, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
    plFree2dGrid(count, nsnap2, nsnap2);
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&count, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    count[0][j] = count[1][j] = cmin + j * dc/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,cmin,cmax);
  plshades(count, 2, CPLOT_NSHADES, NULL, 0.0, 1.0, cmin,cmax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  plschr(0.0, 0.6);
  plmtex("l", 5.0, 0.5, 0.5, "%");
  plmtex("b", 2.0, 0.5, 0.5, "fraction");

  plFree2dGrid(count, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_countfrac(field);		/* Recursive stuff */

  return RETURN_OK;
  }

/****** cplot_modchi2 ********************************************************
PROTO	int cplot_modchi2(fieldstruct *field)
PURPOSE	Plot an x,y map of the average fit chi2/dof.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_modchi2(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   PLFLT	**count,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		cmin,cmax, mc,dc;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		**chi2,
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		**scount,
		i,j, e, n,nt, naxis, nsnap,nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_CHI2) == RETURN_ERROR)
    {
    cplot_end(CPLOT_CHI2);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: #gx#u2#d/d.o.f. map", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  cmin = BIG;
  cmax = -BIG;

/* First pass through the data to find min and max chi2 */
  flag = 0;
  scount = field->count;
  chi2 = field->modchi2;
  nsnap = prefs.context_nsnap;
  nt = nsnap*nsnap;
  for (e=0; e<field->next; e++)
    for (n=0; n<nt; n++)
      {
      dval = scount[e][n] > 0? chi2[e][n] / scount[e][n] : 0.0;
      if (dval < cmin)
        cmin = dval;
      if (dval > cmax)
        cmax = dval;
      }

/* Lower bound to variability in chi2 fraction is 1e-4 */
  if ((mc=(cmin+cmax)/2.0) <  1.0e-10)
    mc = 1.0e-10;
  if ((dc=cmax-cmin) < 1.0e-4*mc)
    dc = 1.0e-4*mc;
  cmin = mc - dc/2.0;
  cmax = mc + dc/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = cmin + (i-0.5) * dc / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real 2D chi2 mapping */
  nsnap2 = nsnap>1? nsnap : 2;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    plAlloc2dGrid(&count, nsnap2, nsnap2);
    for (i=0; i<naxis; i++)
      raw[i] = wcs->naxisn[i]/2.0 + 0.5;
    xstep = wcs->naxisn[0] / (nsnap2-1);
    ystep = wcs->naxisn[1] / (nsnap2-1);
    for (j=0; j<nsnap2; j++)
      for (i=0; i<nsnap2; i++)
        {
        n = (nsnap>1)? j*nsnap + i : 0;
        count[i][j] = (PLFLT)(scount[e][n] > 0? chi2[e][n]/scount[e][n] : 0.0);
        }

    wcsptr[0] = wcs;
    wcsptr[1] = wcsout;
    plshades(count, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
    plFree2dGrid(count, nsnap2, nsnap2);
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&count, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    count[0][j] = count[1][j] = cmin + j * dc/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,cmin,cmax);
  plshades(count, 2, CPLOT_NSHADES, NULL, 0.0, 1.0, cmin,cmax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  plschr(0.0, 0.6);
  plmtex("l", 5.0, 0.5, 0.5, "<#gx#u2#d/d.o.f.>");

  plFree2dGrid(count, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_modchi2(field);		/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_modresi *******************************************************
PROTO	int cplot_modresi(fieldstruct *field)
PURPOSE	Plot an x,y map of the average fit chi2/dof.
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/04/2012
 ***/
int	cplot_modresi(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsout;
   PLFLT	**count,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		cmin,cmax, mc,dc;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[80];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		**resi,
		xmin,ymin,xmax,ymax, xstep,ystep, dval;
   int		naxisn[NAXIS],
		**scount,
		i,j, e, n,nt, naxis, nsnap,nsnap2, flag;

  if (cplot_init(field->rcatname, 1,1, CPLOT_MODRESI) == RETURN_ERROR)
    {
    cplot_end(CPLOT_MODRESI);
    return RETURN_OK;
    }

  wcs = field->wcs[0];
  if (!wcs || wcs->naxis<2)
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcs->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==wcs->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }

  wcsout = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol0(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uField %.24s: map of residuals", field->rtcatname);
  plschr(0.0, 1.0);
  pllab("","", str);
  plwid(0);
  plcol0(7);
  cplot_drawloccoordgrid(wcsout, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol0(15);
  plscmap1n(256);

  cmin = BIG;
  cmax = -BIG;

/* First pass through the data to find min and max chi2 */
  flag = 0;
  scount = field->count;
  resi = field->modresi;
  nsnap = prefs.context_nsnap;
  nt = nsnap*nsnap;
  for (e=0; e<field->next; e++)
    for (n=0; n<nt; n++)
      {
      dval = scount[e][n] > 0? resi[e][n] / scount[e][n] : 0.0;
      if (dval < cmin)
        cmin = dval;
      if (dval > cmax)
        cmax = dval;
      }

/* Lower bound to variability in chi2 fraction is 1e-4 */
  if ((mc=(cmin+cmax)/2.0) <  1.0e-10)
    mc = 1.0e-10;
  if ((dc=cmax-cmin) < 1.0e-4*mc)
    dc = 1.0e-4*mc;
  cmin = mc - dc/2.0;
  cmax = mc + dc/2.0;

  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = cmin + (i-0.5) * dc / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

/* Now the real 2D chi2 mapping */
  nsnap2 = nsnap>1? nsnap : 2;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    plAlloc2dGrid(&count, nsnap2, nsnap2);
    for (i=0; i<naxis; i++)
      raw[i] = wcs->naxisn[i]/2.0 + 0.5;
    xstep = wcs->naxisn[0] / (nsnap2-1);
    ystep = wcs->naxisn[1] / (nsnap2-1);
    for (j=0; j<nsnap2; j++)
      for (i=0; i<nsnap2; i++)
        {
        n = (nsnap>1)? j*nsnap + i : 0;
        count[i][j] = (PLFLT)(scount[e][n] > 0? resi[e][n]/scount[e][n] : 0.0);
        }

    wcsptr[0] = wcs;
    wcsptr[1] = wcsout;
    plshades(count, nsnap2, nsnap2, NULL,
	     xstep/2.0+0.5, wcs->naxisn[0]-xstep/2.0+0.5,
             ystep/2.0+0.5, wcs->naxisn[1]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 0, distort_map, wcsptr);
    plFree2dGrid(count, nsnap2, nsnap2);
    plcol0(7);
    plwid(lwid);
    cplot_drawbounds(wcs, wcsout);
    }

/* Draw left colour scale */
  plAlloc2dGrid(&count, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    count[0][j] = count[1][j] = cmin + j * dc/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,cmin,cmax);
  plshades(count, 2, CPLOT_NSHADES, NULL, 0.0, 1.0, cmin,cmax, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol0(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  plschr(0.0, 0.6);
  plmtex("l", 5.0, 0.5, 0.5, "normalized residuals");

  plFree2dGrid(count, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcsout);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_modresi(field);		/* Recursive stuff */

  return RETURN_OK;
  }


