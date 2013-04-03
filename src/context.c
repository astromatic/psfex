/*
*				context.c
*
* Manage observation contexts.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		20/07/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "context.h"
#include "wcs/poly.h"
#include "psf.h"
#include "field.h"

/****** context_init *********************************************************
PROTO   contextstruct *context_init(char **names, int *group, int ndim,
			int *degree, int ngroup, int pcexflag)
PURPOSE Allocate and initialize a context structure.
INPUT   Pointer to an array of context names,
	Pointer to an array of group indices,
	Number of dependency parameters,
	Pointer to an array of group degrees,
	Number of groups,
	Principal Component exclusion flag.
OUTPUT  Pointer to an allocated context structure.
NOTES   See prefs.h.
AUTHOR  E. Bertin (IAP)
VERSION 13/02/2009
*/
contextstruct	*context_init(char **names, int *group, int ndim, int *degree,
	 int ngroup, int pcexflag)
  {
   contextstruct	*context;
   int			*groupflag,
			d,d2,g,g2, pcflag;

  QCALLOC(context, contextstruct, 1);
  QCALLOC(context->name, char *, ndim);
  QCALLOC(context->group, int, ndim);
  QCALLOC(context->degree, int, ndim);
  QCALLOC(context->pcflag, int, ndim);
  QCALLOC(groupflag, int, ndim);
/* Copy context names and group indices ...*/
/* ... and remove Principal Components if asked to */
  d2=0;
  for (d=0; d<ndim; d++)
    {
    if ((pcflag = !wstrncmp(names[d], "HIDDEN?", 80)))
      {
      context->npc++;
      if (!pcexflag)
        context->pcflag[d] = 1;
      }
    if (!pcexflag || !pcflag)
      {
      QMALLOC(context->name[d2], char, 80);
      strncpy(context->name[d2], names[d], 80);
      context->group[d2] = group[d];
      groupflag[group[d]-1]++;
      d2++;
      } 
    }

  context->ncontext = d2;

/* Reorganize groups to remove missing group indices */
  g2 = 0;
  for (g=0; g<ngroup; g++)
    {
    if (groupflag[g])
      {
      context->degree[g2] = degree[g];
      for (d2=0; d2<context->ncontext; d2++)
        if (context->group[d2]==g+1)
          context->group[d2]=g2+1;
      g2++;
      }
    }

  context->ngroup = g2;
  free(groupflag);

  return context;
  }


/****** context_apply ********************************************************
PROTO	void context_apply(contextstruct *context, psfstruct *psf,
		fieldstruct **fields, int ext, int catindex, int ncat)
PURPOSE	"Apply" hidden dependencies to a set of PSFs.
INPUT	Pointer to the full context,
	Pointer to the PSF which depends (or not) on hidden dependencies,
	Pointer to an array of fields where the final PSFs are to be copied.
	Current extension,
	Starting catalog index,
	Number of catalogs.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/11/2009
 ***/
void context_apply(contextstruct *context, psfstruct *psf,
		fieldstruct **fields, int ext, int catindex, int ncat)
  {
   psfstruct		*psf2;
   polystruct		*poly, *poly2;
   contextstruct	*context2;
   const double		prime[]={2.0,3.0,5.0,7.0,11.0,13.0,17.0,19.0,23.0,
				29.0,31.0,37.0,41.0,43.0,47.0,53.0,59.0,61.0};
   double		dpos[MAXCONTEXT],
			dval;
   float		*comp, *comp2, *comp2t, *bcoeff2;
   int			*polycopyflag, *polycopyflagt,
			c,c2, e, i, n,n2, p, npc, npix,comp2size,ncontext2,
			bcoeff2size;

  ncat += catindex;
  for (p=catindex; p<ncat; p++)
    if (ext==ALL_EXTENSIONS)
      for (e=0; e<fields[p]->next; e++)
        fields[p]->psf[e] = psf_copy(psf);
    else
      fields[p]->psf[ext] = psf_copy(psf);

/* No PC dependency: the PSF is simply duplicated */
  if (!context->npc)
    return;

  if (sizeof(prime)/sizeof(double) < context->ncontext)
    error(EXIT_FAILURE, "*Internal Error*: ",
		"not enough prime numbers in context_apply()");

/* Prepare component copying process */
/* Use the first PSF as template */
  poly = psf->poly;
  if (poly->ndim != context->ncontext)
    error(EXIT_FAILURE, "*Internal Error*: ",
		"poly->ndim != context->ncontext in context_apply()");

/* Create a new context with the PCs removed */
  context2 = context_init(context->name, context->group, context->ncontext,
			context->degree, context->ngroup, CONTEXT_REMOVEHIDDEN);
  ncontext2 = context2->ncontext;
  poly2 = poly_init(context2->group, context2->ncontext, context2->degree,
		context2->ngroup);

/* Identify the polynomial terms that have to be merged using prime numbers */
  c2 = 0;
  for (c=0; c<context->ncontext; c++)
    if (context->pcflag[c])
      dpos[c] = 1.0;
    else
      dpos[c] = prime[c2++];
  poly_func(poly, dpos);
  for (c=0; c<context2->ncontext; c++)
    dpos[c] = prime[c];
  poly_func(poly2, dpos);

/* The polycopyflag matrix flags terms that will be merged */
  QMALLOC(polycopyflag, int, poly2->ncoeff*poly->ncoeff);
  polycopyflagt = polycopyflag;
  for (n2=0; n2<poly2->ncoeff; n2++)
    {
    dval = poly2->basis[n2];
    for (n=0; n<poly->ncoeff; n++)
      *(polycopyflagt++) = (fabs(poly->basis[n] - dval) < 0.1);
    }

/* Merge PSF components for each PSF */
  npix = psf->size[0]*psf->size[1];
  comp2size = npix*poly2->ncoeff;
  bcoeff2 = NULL;			/* To avoid gcc -Wall warnings */
  bcoeff2size = psf->nbasis*poly2->ncoeff;
  npc = context->npc;
  for (p=catindex; p<ncat; p++)
    {
/*--- Update PC component values */
    c2 = 0;
    for (c=0; c<context->ncontext; c++)
      if (context->pcflag[c])
        dpos[c] = (context->pc[p*npc+c2++]-psf->contextoffset[c])
		/psf->contextscale[c];
      else
        dpos[c] = 1.0;
    poly_func(poly, dpos);
    QCALLOC(comp2, float, comp2size);
    if (psf->basiscoeff)
      QCALLOC(bcoeff2, float, bcoeff2size);
    polycopyflagt = polycopyflag;
    for (n2=0; n2<poly2->ncoeff; n2++)
      for (n=0; n<poly->ncoeff; n++)
        if (*(polycopyflagt++))
          {
          comp = psf->comp + n*npix;
          comp2t = comp2 + n2*npix;
          dval = poly->basis[n];
          for (i=npix; i--;)
            *(comp2t++) += dval**(comp++);
          if (psf->basiscoeff)
            bcoeff2[n2] += dval*psf->basiscoeff[n];
          }
/*-- Replace the new PSF components */
    if (ext==ALL_EXTENSIONS)
      for (e=0;e<fields[p]->next; e++)
        {
        fields[p]->psf[e] = psf2 = psf_inherit(context2, psf);
        free(psf2->comp);
        if (e)
          {
          QMEMCPY(comp2, psf2->comp, float, comp2size);
          if (psf->basiscoeff)
            QMEMCPY(bcoeff2, psf2->basiscoeff, float, bcoeff2size);
          }
        else
          {
          psf2->comp = comp2;
          if (psf->basiscoeff)
            psf2->basiscoeff = bcoeff2;
          }
        }
    else
      {
      fields[p]->psf[ext] = psf2 = psf_inherit(context2, psf);
      free(psf2->comp);
      psf2->comp = comp2;
      if (psf->basiscoeff)
        psf2->basiscoeff = bcoeff2;
      }
    }

  return;
  }


/****** context_end ***********************************************************
PROTO   void context_end(contextstruct *context)
PURPOSE free memory allocated by a complete context structure.
INPUT   Context structure pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 29/10/2008
*/
void	context_end(contextstruct *context)
  {
   int	d;

  for (d=0; d<context->ncontext; d++)
    free(context->name[d]);
  free(context->name);
  free(context->group);
  free(context->degree);
  free(context->pcflag);
  free(context->pc);
  free(context);

  return;
  }
