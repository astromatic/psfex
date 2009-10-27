  /*
 				field.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Handling of multiple PSFs.
*
*	Last modify:	27/10/2009
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
#include	"check.h"
#include	"fitswcs.h"
#include	"prefs.h"
#include	"psf.h"
#include	"field.h"

/****** field_init ************************************************************
PROTO	fieldstruct *field_init(char *catname, int next)
PURPOSE	Allocate and initialize a PSF MEF structure (groups of PSFs).
INPUT	Catalog filename.
OUTPUT  fieldstruct pointer.
NOTES   .
AUTHOR  E. Bertin (IAP)
VERSION 30/03/2009
 ***/
fieldstruct	*field_init(char *catname)
  {
   fieldstruct	*field;
   catstruct	*cat;
   tabstruct	*tab, *imatab;
   keystruct	*key;
   char		*pstr;
   int		e, next, next0, ntab, countsize;

  QCALLOC(field, fieldstruct, 1);
/* Compute the number of valid input extensions */
  if (!(cat = read_cat(catname)))
    error(EXIT_FAILURE, "*Error*: cannot open ", catname);
  tab = cat->tab;
  next0 = 0;
  for (ntab = 0 ; ntab<cat->ntab; ntab++, tab = tab->nexttab)
    {
/*--  Check for the next valid image extension */
    if ((tab->naxis != 2)
	|| (strncmp(tab->xtension, "BINTABLE", 8)
		&& strncmp(tab->xtension, "ASCTABLE", 8))
	|| (strncmp(tab->extname, "LDAC_OBJECTS", 8)
		&& strncmp(tab->extname, "OBJECTS", 8)))
      continue;
    next0++;
    }
  field->next = next0;
  QMALLOC(field->psf, psfstruct *, next0);
  strcpy(field->catname, catname);
/* A short, "relative" version of the filename */
  if (!(field->rcatname = strrchr(field->catname, '/')))
    field->rcatname = field->catname;
  else
    field->rcatname++;
  strcpy(field->rtcatname, field->rcatname);
  if ((pstr=strrchr(field->rtcatname, '.')))
    *pstr = '\0';

  if (!next0)
    {
    field_end(field);
    error(EXIT_FAILURE,"*Error*: No SExtractor FITS-LDAC catalog found in ",
        catname);
    }

  QMALLOC(field->wcs, wcsstruct *, next0);
/* Compute the number of valid input extensions */
  tab = cat->tab;
  next = 0;
  for (ntab = 0 ; ntab<cat->ntab; ntab++, tab = tab->nexttab)
/*--  Check for the next valid FITS extension */
    if ((!strcmp("LDAC_IMHEAD",tab->extname))
	&& (key=read_key(tab, "Field Header Card")))
      {
/*---- Create a new table from scratch to hold the image header */
      imatab = new_tab("Image header");
      free(imatab->headbuf);
      imatab->headnblock = 1 + (key->nbytes-1)/FBSIZE;
      QCALLOC(imatab->headbuf, char, imatab->headnblock*FBSIZE);
      memcpy(imatab->headbuf, key->ptr, key->nbytes);
      imatab->cat = cat;
      readbasic_head(imatab);
      field->wcs[next++] = read_wcs(imatab);
      if (!imatab->headbuf
	|| fitsread(imatab->headbuf, "OBJECT  ", field->ident,
	H_STRING,T_STRING)!= RETURN_OK)
        strcpy(field->ident, "no ident");
      free_tab(imatab);
      }
    else if ((!strcmp("LDAC_OBJECTS", tab->extname)
	||  !strcmp("OBJECTS", tab->extname)) && tab->naxis == 2)    
      field->ndet += tab->naxisn[1];

  free_cat(&cat, 1);

  field_locate(field);
  QCALLOC(field->ccat, catstruct *, MAXCHECK);
  countsize = prefs.context_nsnap*prefs.context_nsnap;
  QMALLOC(field->lcount, int *, next0);
  QMALLOC(field->acount, int *, next0);
  for (e=0; e<next0; e++)
    {
    QCALLOC(field->lcount[e], int, countsize);
    QCALLOC(field->acount[e], int, countsize);
    }

  return field;
  }


/****** field_end *************************************************************
PROTO	void field_end(fieldstruct *field)
PURPOSE	Free a PSF MEF structure (groups of PSFs).
INPUT	Pointer to the fieldstruct.
OUTPUT  -.
NOTES   .
AUTHOR  E. Bertin (IAP)
VERSION 30/03/2009
 ***/
void	field_end(fieldstruct *field)
  {
   int	ext;

  for (ext=0; ext<field->next; ext++)
    {
    psf_end(field->psf[ext]);
    end_wcs(field->wcs[ext]);
    free(field->lcount[ext]);
    free(field->acount[ext]);
    }
  free(field->psf);
  free(field->wcs);
  free(field->ccat);
  free(field->lcount);
  free(field->acount);
  free(field);

  return;
  }


/****** field_locate *********************************************************
PROTO   void field_locate(fieldstruct *field)
PURPOSE Compute field position, scale and footprint.
INPUT   Pointer to field structure.
OUTPUT  A pointer to the created field structure.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 14/03/2008
*/
void	field_locate(fieldstruct *field)
  {
   wcsstruct		*wcs;
   double		*scale[NAXIS],*scalet[NAXIS],
			*wcsmean,
			cosalpha,sinalpha, sindelta, dist, maxradius;
   int			i, e, lat,lng, naxis;

/* Some initializations */
  cosalpha = sinalpha = sindelta = 0.0;
  wcs = field->wcs[0];
  naxis = wcs->naxis;
  wcsmean = field->meanwcspos;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(scale[i], double, field->next);
    scalet[i] = scale[i];
    wcsmean[i] = 0.0;
    }

/* Go through each extension */
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    lng = wcs->lng;
    lat = wcs->lat;
/*-- Locate set */
    if (lat != lng)
      {
      cosalpha += cos(wcs->wcsscalepos[lng]*DEG);
      sinalpha += sin(wcs->wcsscalepos[lng]*DEG);
      sindelta += sin(wcs->wcsscalepos[lat]*DEG);
      }
    for (i=0; i<naxis; i++)
      {
      if (lat==lng || (i!=lng && i!=lat))
        wcsmean[i] += wcs->wcsscalepos[i];
      *(scalet[i]++) = wcs->wcsscale[i];
      }
    }

/* Now make the stats on each axis */
  lng = field->wcs[0]->lng;
  lat = field->wcs[0]->lat;
  for (i=0; i<naxis; i++)
    {
    if (lat!=lng && (i==lng))
      {
      wcsmean[i] = atan2(sinalpha/field->next,cosalpha/field->next)/DEG;
      wcsmean[i] = fmod(wcsmean[i]+360.0, 360.0);
      }
    else if (lat!=lng && (i==lat))
      wcsmean[i] = asin(sindelta/field->next)/DEG;
    else
      wcsmean[i] /= field->next;
    field->meanwcsscale[i] = dhmedian(scale[i], field->next);
    }

/* Compute the radius of the field and mean airmass */
  maxradius = 0.0;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
/*-- The distance is the distance to the center + the diagonal of the image */
    dist = wcs_dist(wcs, wcs->wcsscalepos, field->meanwcspos)
		+ wcs->wcsmaxradius;
    if (dist>maxradius)
      maxradius = dist;
    }

  field->maxradius = maxradius;

/* Free memory */
  for (i=0; i<naxis; i++)
    free(scale[i]);

  return;
  }


/****** dhmedian **************************************************************
PROTO	double   dhmedian(double *ra, int n)
PURPOSE	Compute the median of an array of doubles, using the Heapsort
	algorithm (based on Num.Rec algo.).
INPUT	Pointer to the array,
	Number of array elements.
OUTPUT	Median of the array.
NOTES	Warning: the order of input data is modified!.
AUTHOR	E. Bertin (IAP)
VERSION	22/07/2002
 ***/
double   dhmedian(double *ra, int n)

  {
   int		l, j, ir, i;
   double	rra;


  if (n<2)
    return *ra;
  ra--;
  for (l = ((ir=n)>>1)+1;;)
    {
    if (l>1)
      rra = ra[--l];
    else
      {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1)
        {
        ra[1] = rra;
        return n&1? ra[n/2+1] : (ra[n/2]+ra[n/2+1])/2.0;
        }
      }
    for (j = (i=l)<<1; j <= ir;)
      {
      if (j < ir && ra[j] < ra[j+1])
        ++j;
      if (rra < ra[j])
        {
        ra[i] = ra[j];
        j += (i=j);
        }
      else
        j = ir + 1;
      }
    ra[i] = rra;
    }

/* (the 'return' is inside the loop!!) */
  }


/****** field_count ***********************************************************
PROTO	void field_count(fieldstruct **fields, setstruct *set, int counttype)
PURPOSE	Count the number of sources (samples) per image area.
INPUT	Pointer to an array of fieldstruct pointers,
	Pointer to the set to be counter,
	Sample type.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 30/03/2009
 ***/
void	field_count(fieldstruct **fields, setstruct *set, int counttype)
  {
   fieldstruct	*field;
   samplestruct	*sample;
   int		c,e,n,s, w,h, x,y, size;

  sample = set->sample;
  size = (double)prefs.context_nsnap;
  for (s=set->nsample; s--; sample++)
    {
    c = sample->catindex;
    e = sample->extindex;
    field = fields[c];
    w = field->wcs[e]->naxisn[0];
    h = field->wcs[e]->naxisn[1];
    x = (int)((sample->x-0.5)*size) / w;
    if (x<0)
      x = 0;
    else if (x>=size)
      x = size-1;
    y = (int)((sample->y-0.5)*size) / h;
    if (y<0)
      y = 0;
    else if (y>=size)
      y = size-1;
    n = y*size+x;
    if ((counttype & COUNT_LOADED))
      fields[c]->lcount[e][n]++;
    if ((counttype & COUNT_ACCEPTED))
      fields[c]->acount[e][n]++;
    }

  return;
  }


/****** field_psfsave *********************************************************
PROTO   void	field_psfsave(fieldstruct *field, char *filename)
PURPOSE Save PSF data as a Multi-extension FITS file.
INPUT   Pointer to the PSF structure,
	Filename,
	Extension number,
	Number of extensions.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 27/10/2009
 ***/
void	field_psfsave(fieldstruct *field, char *filename)
  {
   catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   psfstruct	*psf;
   float	zero = 0.0;
   char		*head,
		str[80];
   int		i, ext, temp;

  cat = new_cat(1);
  init_cat(cat);
  sprintf(cat->filename, filename);
  if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: cannot open for writing ", cat->filename);
/* Write primary HDU */
  save_tab(cat, cat->tab);

  for (ext=0; ext<field->next; ext++)
    {
    psf = field->psf[ext];
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

/*-- Add and write important scalars as FITS keywords */
    addkeywordto_head(tab, "PSF_FWHM", "PSF FWHM");
    fitswrite(head, "PSF_FWHM", psf->samples_accepted? &psf->fwhm : &zero,
	H_FLOAT, T_FLOAT);
    addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data");
    fitswrite(head, "PSF_SAMP", psf->samples_accepted? &psf->pixstep : &zero,
	H_FLOAT, T_FLOAT);
    addkeywordto_head(tab, "PSFNAXIS", "Dimensionality of the PSF data");
    fitswrite(head, "PSFNAXIS", &psf->dim, H_INT, T_LONG);
    for (i=0; i<psf->dim; i++)
      {
      sprintf(str, "PSFAXIS%1d", i+1);
      addkeywordto_head(tab, str, "Number of element along this axis");
      fitswrite(head, str, &psf->size[i], H_INT, T_LONG);
      }

/*-- Create and fill the arrays */
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
/*-- But don't touch my arrays!! */
    blank_keys(tab);
    free_tab(tab);
    }

  free_cat(&cat , 1);

  return;
  }


