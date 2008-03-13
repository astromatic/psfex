  /*
 				psfmef.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Handling of multiple PSFs.
*
*	Last modify:	13/03/2008
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
#include	"fitswcs.h"
#include	"prefs.h"
#include	"psf.h"
#include	"psfmef.h"

/****** psfmef_init ***********************************************************
PROTO	psfmefstruct *psfmef_init(char *catname, int next)
PURPOSE	Allocate and initialize a PSF MEF structure (groups of PSFs).
INPUT	Catalog filename.
OUTPUT  psfmefstruct pointer.
NOTES   .
AUTHOR  E. Bertin (IAP)
VERSION 13/03/2008
 ***/
psfmefstruct	*psfmef_init(char *catname)
  {
   psfmefstruct	*psfmef;
   catstruct	*cat;
   tabstruct	*tab, *imatab;
   keystruct	*key;
   int		next, ntab;

  QCALLOC(psfmef, psfmefstruct, 1);
/* Compute the number of valid input extensions */
  if (!(cat = read_cat(catname)))
    error(EXIT_FAILURE, "*Error*: cannot open ", catname);
  tab = cat->tab;
  next = 0;
  for (ntab = 0 ; ntab<cat->ntab; ntab++, tab = tab->nexttab)
    {
/*--  Check for the next valid image extension */
    if ((tab->naxis != 2)
	|| (strncmp(tab->xtension, "BINTABLE", 8)
		&& strncmp(tab->xtension, "ASCTABLE", 8))
	|| (strncmp(tab->extname, "LDAC_OBJECTS", 8)
		&& strncmp(tab->extname, "OBJECTS", 8)))
      continue;
    next++;
    }
  psfmef->next = next;
  QMALLOC(psfmef->psf, psfstruct *, next);
  strcpy(psfmef->catname, catname);
/* A short, "relative" version of the filename */
  if (!(psfmef->rcatname = strrchr(psfmef->catname, '/')))
    psfmef->rcatname = psfmef->catname;
  else
    psfmef->rcatname++;

  QMALLOC(psfmef->wcs, wcsstruct *, next);
/* Compute the number of valid input extensions */
  tab = cat->tab;
  next = 0;
  for (ntab = 0 ; ntab<cat->ntab; ntab++, tab = tab->nexttab)
    {
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
      psfmef->wcs[next++] = read_wcs(imatab);
      free_tab(imatab);
      }
      continue;
    }

  free_cat(&cat, 1);

  return psfmef;
  }


/****** psfmef_end ***********************************************************
PROTO	void psfmef_end(psfmefstruct *psfmef)
PURPOSE	Free a PSF MEF structure (groups of PSFs).
INPUT	Pointer to the psfmefstruct.
OUTPUT  -.
NOTES   .
AUTHOR  E. Bertin (IAP)
VERSION 13/03/2008
 ***/
void	psfmef_end(psfmefstruct *psfmef)
  {
   int	ext;

  for (ext=0; ext<psfmef->next; ext++)
    {
    psf_end(psfmef->psf[ext]);
    end_wcs(psfmef->wcs[ext]);
    }
  free(psfmef->psf);
  free(psfmef->wcs);
  free(psfmef);

  return;
  }


/****** psfmef_save ***********************************************************
PROTO   void	psfmef_save(psfmefstruct *psfmef, char *filename)
PURPOSE Save PSF data as a Multi-extension FITS file.
INPUT   Pointer to the PSF structure,
	Filename,
	Extension number,
	Number of extensions.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 11/03/2008
 ***/
void	psfmef_save(psfmefstruct *psfmef, char *filename)
  {
   catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   psfstruct	*psf;
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

  for (ext=0; ext<psfmef->next; ext++)
    {
    psf = psfmef->psf[ext];
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

